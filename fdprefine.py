#!/dls/science/groups/i23/pyenvs/tihana_conda/bin/python
import os
import re
from iotbx import reflection_file_reader 
import iotbx
import gemmi
import subprocess 
import copy 
import plotly.graph_objects as go
import plotly.io as pio

from multiprocessing import Pool
class refinefdoubleprime():
  def __init__(self):
    self.cpus = os.cpu_count() - 1
    self.pool = Pool(self.cpus)

    try:
      subprocess.run(["phenix.about"])
      print("Found Phenix installation")
    except:
      print("Cannot find Phenix installation. Try to run module load phenix")

    mtzIn = input("File location for MTZ: ")
    self.pdbIn = input("File location for PDB: ")
    self.projIn = input("Name of project: ")
    genMonomerLib = input("Do you have ligands in the PDB file? (y/n) ").lower()
    if genMonomerLib == "y":
      print("Generating monomer library, this should only take a few minutes...\n")
      subprocess.run(["phenix.ready_set", f"{self.pdbIn}"])
      pdbInBase, pdbInExt = self.pdbIn.rsplit('.', 1)
      pdbInBase = os.path.basename(pdbInBase)
      self.ligandIn = str(pdbInBase + ".ligands.cif")
    else:
      self.ligandIn = str(None)

    elementsToTry = input("Which elements to try, comma separated: ")
    self.elements = [x.strip() for x in elementsToTry.split(',')]

    mtzInfo = reflection_file_reader.any_reflection_file(mtzIn)
    millerArray = mtzInfo.as_miller_arrays()
    mtzInfo.file_content()
    mtzobj = iotbx.mtz.object(mtzIn)
    self.space_group = mtzobj.space_group_name()
    csym = mtzobj.crystals()[0].crystal_symmetry()
    unit_cell = csym.unit_cell()
    self.unit_cell_strip = (str(unit_cell)).strip('()').replace(',', '')
    self.mtz = gemmi.read_mtz_file(mtzIn)
    self.WV = (self.mtz.dataset(1).wavelength)
    self.scrapedData = []
    
  def wavelength_to_eV(self):
    h = 6.62607015e-34
    c = 2.99792458e8
    eV = 1.602176634e-19
    self.energyeV = (h * c) / (float(self.WV) * 1e-10) / eV

  def lookup_fprime(self, element):
    lookupPath = os.path.join("lookup", f"{element}.dat")
    closestEnergy = None
    self.closestValues = None
    with open(lookupPath, 'r') as lookupFile:
        for line in lookupFile:
            parts = line.split()
            currentEnergy, fPrime, fDoublePrime = map(float, parts)
            if closestEnergy is None or abs(self.energy - currentEnergy) < abs(self.energy - closestEnergy):
                closestEnergy = currentEnergy
                self.closestValues = (fPrime, fDoublePrime)

  def change_elem(self):
    pdbInBase, pdbInExt = self.pdbIn.rsplit('.', 1)
    pdbInBase = os.path.basename(pdbInBase)
    print(f"Running on {str(pdbInBase)}")
    pdbLinesWrite = []
    self.pdbList = []

    with open(self.pdbIn, 'r') as file:
      pdbLinesWrite = file.readlines()

    toChange = []
    for line in pdbLinesWrite:
      if line.startswith("HETATM"):
        print(line.strip())
        changeHETATM = input("Would you like to run refinement on this HETATM? ")
        if changeHETATM.lower() in ('y', 'yes'):
          toChange.append(line)

    for element in self.elements:
      toFDPRefine = []
      toWrite = copy.deepcopy(pdbLinesWrite)
      for i, line in enumerate(pdbLinesWrite):
        if line in toChange:
          elementSymbol = element.rjust(2)[:2] if len(element) == 2 else ' ' + element
          residueName = element.rjust(3)[:3]
          elementName = element.ljust(4)[:4]
          lineList = list(line)
          lineList[12:16] = elementName
          lineList[17:20] = residueName
          lineList[76:78] = elementSymbol
          toWrite[i] = ''.join(lineList)
          toFDPRefine.append([element, line[21:22].strip(), line[23:26].strip()])
          
      with open(f"{pdbInBase}_{element}.{pdbInExt}", 'w') as pdbOut:
        for line in toWrite:
          if not line.endswith('\n'):
            line += '\n'
          pdbOut.write(line)
        self.pdbList.append((f"{pdbInBase}_{element}.{pdbInExt}", element, toFDPRefine))
  
  def scrapeLastAnomalousGroupData(self, ele):
    log_file_path = (f'{self.projIn}_fdp_{ele}_1.log')
    with open(log_file_path, 'r') as file:
      content = file.read()

    last_macro_cycle_index = content.rfind('MACRO_CYCLE')

    content_after_macro_cycle = content[last_macro_cycle_index:]
      
    pattern = re.compile(
        r'Anomalous scatterer group:\s+'
        r'Selection: "chain (?P<chain>[A-Za-z]) and resid (?P<resid>\d+) and element (?P<element>[A-Za-z]+)"\s+'
        r'Number of selected scatterers: \d+\s+'
        r'f_prime: +[+-]?\d+\.\d+\s+'
        r'f_double_prime: (?P<f_double_prime>[+-]?\d+\.\d+)'
    )

    matches = pattern.findall(content_after_macro_cycle)
    data = [(m[0], int(m[1]), m[2], float(m[3]), float(self.closestValues[1]), float(abs(float(self.closestValues[1]) - (float(m[3]))))) for m in matches]

    self.scrapedData.append(data)

  def makeTable(self):
    header = ['Chain', 'ResidID', 'Element', 'Refined FDP', "Theoretical FDP", "Absolute Difference"]
    flattenedData = [item for sublist in self.scrapedData for item in sublist]
    tableData = list(map(list, zip(*flattenedData)))
    
    fig = go.Figure(data=[go.Table(
      header=dict(values=header,
                  fill_color='paleturquoise',
                  align='left'),
      cells=dict(values=tableData,
                fill_color='lavender',
                align='left'))
    ])
    pio.write_html(fig, f'{self.projIn}.html')


  # Create bpos eff file
  def runBPos(self, pdbIn, elementIn):
    with open(f'bposEffParam_{elementIn}.eff', 'w') as file:
      file.write(f'''refinement {{
      crystal_symmetry {{  
        unit_cell = {self.unit_cell_strip}
        space_group = {self.space_group}
      }}
      input {{  
        pdb {{  
          file_name = "{pdbIn}"  
        }}  
        xray_data {{  
          file_name = "{self.mtzIn}"
          labels = IMEAN,SIGIMEAN  
          r_free_flags {{  
            file_name = "{self.mtzIn}"
            label = FreeR_flag  
            test_flag_value = 0  
          }}  
        }}    
		    monomers {{
		      file_name = {self.ligandIn}
        }}
      }}  
      output {{  
        prefix = """{self.projIn}_bpos_{str(elementIn)}"""   
        job_title = """{self.projIn}"""  
        serial_format = "%d"
        write_def_file = False  
      }}  
      electron_density_maps {{  
        map_coefficients {{  
          map_type = 2mFo-DFc  
          mtz_label_amplitudes = 2FOFCWT  
          mtz_label_phases = PH2FOFCWT  
          fill_missing_f_obs = True  
        }}  
        map_coefficients {{  
          map_type = 2mFo-DFc  
          mtz_label_amplitudes = 2FOFCWT_no_fill  
          mtz_label_phases = PH2FOFCWT_no_fill  
        }}  
        map_coefficients {{  
          map_type = mFo-DFc  
          mtz_label_amplitudes = FOFCWT  
          mtz_label_phases = PHFOFCWT  
        }}  
        map_coefficients {{  
          map_type = anomalous  
          mtz_label_amplitudes = ANOM  
          mtz_label_phases = PHANOM  
        }}  
      }}  
      refine {{  
        strategy = *individual_sites individual_sites_real_space rigid_body \  
                  *individual_adp group_adp tls occupancies group_anomalous  
      }}  
      main {{  
        number_of_macro_cycles = 5  
        wavelength = {self.WV}
        nproc = {self.cpus}
      }}   
      }}''')

    subprocess.run(["phenix.refine", f"bposEffParam_{elementIn}.eff"])

  def runFdp(self, elementIn, toFDPRefine):
    anomalousScatterersStr = "anomalous_scatterers {\n"
    for elemID, chainID, residID in toFDPRefine:
        groupStr = f"""          group {{
              selection = chain {chainID} and resid {residID} and element {elemID}
              f_prime = {self.closestValues[0]}
              refine = f_prime *f_double_prime
            }}\n"""
        anomalousScatterersStr += groupStr
    anomalousScatterersStr += "        }"
    
    with open(f'fdpEffParam_{elementIn}.eff', 'w') as file:
      file.write(f'''refinement {{  
			crystal_symmetry {{
        unit_cell = {self.unit_cell_strip}
        space_group = {self.space_group}
      }}  
      input {{  
        pdb {{  
          file_name = "{self.projIn}_bpos_{elementIn}_1.pdb"  
        }}  
        xray_data {{  
          file_name = "{self.mtzIn}"
          labels = I(+),SIGI(+),I(-),SIGI(-)
          r_free_flags {{  
            file_name = "{self.mtzIn}"
            label = FreeR_flag  
            test_flag_value = 0  
          }}  
        }}
        monomers {{
		      file_name = {self.ligandIn}
        }}
      }}  
      output {{  
        prefix = """{self.projIn}_fdp_{elementIn}"""  
        job_title = """{self.projIn}_{elementIn}"""  
        serial_format = "%d"
        write_def_file = False  
      }}  
      electron_density_maps {{  
        map_coefficients {{  
          map_type = 2mFo-DFc  
          mtz_label_amplitudes = 2FOFCWT  
          mtz_label_phases = PH2FOFCWT  
          fill_missing_f_obs = True  
        }}  
        map_coefficients {{  
          map_type = 2mFo-DFc  
          mtz_label_amplitudes = 2FOFCWT_no_fill  
          mtz_label_phases = PH2FOFCWT_no_fill  
        }}  
        map_coefficients {{  
          map_type = mFo-DFc  
          mtz_label_amplitudes = FOFCWT  
          mtz_label_phases = PHFOFCWT  
        }}  
        map_coefficients {{  
          map_type = anomalous  
          mtz_label_amplitudes = ANOM  
          mtz_label_phases = PHANOM  
        }}  
      }}  
      refine {{  
        strategy = individual_sites individual_sites_real_space rigid_body \   
                  individual_adp group_adp tls occupancies *group_anomalous  
      	{anomalousScatterersStr}
      }}
      main {{  
        number_of_macro_cycles = 5  
        wavelength = {self.WV}
      }}  
      }}''')
    
    subprocess.run(["phenix.refine", f"fdpEffParam_{elementIn}.eff"])      

if __name__ == "__main__":  
  run = refinefdoubleprime()
  run.change_elem()
  run.wavelength_to_eV()
  for pdb, ele, tfdpr in run.pdbList:
    run.lookup_fprime()
    run.runBPos(pdb, ele)
    run.runFdp(ele, tfdpr)
    run.scrapeLastAnomalousGroupData(ele)
  run.makeTable()