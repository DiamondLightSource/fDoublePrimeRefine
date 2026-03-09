#!/dls/science/groups/i23/pyenvs/tihana_conda/bin/python
import os
import shutil
import re
from iotbx import reflection_file_reader
import iotbx
import gemmi
import subprocess
import copy
import plotly.graph_objects as go
import plotly.io as pio
from multiprocessing import Pool
from halo import Halo
import numpy as np

lookup_path = "/dls/science/groups/i23/scripts/chris/phenix_automation/lookup"
cwd = os.getcwd()


class refinefdoubleprime:
    def __init__(self):
        self.cpus = os.cpu_count() - 1

        ### should load phenix as part of the setup script, if it doesn't uncomment below to include a check
        # try:
        #   subprocess.run(["phenix.about"])
        #   print("Found Phenix installation")
        # except:
        #   print("Cannot find Phenix installation. Try to run module load phenix")
        self.toignore = ("HOH", "BOG", "IXX", "GOL", "PEG", "ATP")
        self.mtzIn = input("File location for MTZ: ")
        if os.path.exists(self.mtzIn):
            pass
        else:  
            shutil.copy(self.mtzIn, os.path.join(cwd, os.path.basename(self.mtzIn)))
        self.mtzIn = os.path.basename(self.mtzIn)
        self.pdbIn = input("File location for PDB: ")
        if os.path.exists(self.pdbIn):
            pass
        else:
            shutil.copy(self.pdbIn, os.path.join(os.path.basename(self.pdbIn)))
        self.pdbIn = os.path.basename(self.pdbIn)
        if os.path.exists(os.path.join(self.pdbIn + ".orig")):
            pass
        else:
            shutil.copy(self.pdbIn, self.pdbIn + ".orig")
        self.projIn = input("Name of project: ")

        with open(self.pdbIn, "r") as infile:
            lines = infile.readlines()

        cleanedLines = []
        skip = 0
        for i, line in enumerate(lines):
            if skip > 0:
                skip -= 1
                continue

            if line.startswith("ANISOU"):
                continue

            if "TLS DETAILS." in line:
                skip = 17
                continue

            if "TLS DETAILS" in line:
                skip = 2
                continue

            if "TLS GROUP :" in line:
                skip = 16
                continue

            cleanedLines.append(line)

        with open(self.pdbIn, "w") as outfile:
            outfile.writelines(cleanedLines)

        genMonomerLib = input("Do you have ligands in the PDB file? (y/n) ").lower()
        if genMonomerLib == "y":
            with Halo(
                text="Generating monomer library, this should only take a few minutes...",
                spinner="toggle",
            ):
                logFile = "monomerlib_output.log"
                with open(logFile, "a") as log:
                    subprocess.run(
                        ["phenix.ready_set", f"{self.pdbIn}"], stdout=log, stderr=log
                    )
                    pdbInBase, pdbInExt = self.pdbIn.rsplit(".", 1)
                    pdbInBase = os.path.basename(pdbInBase)
                    self.ligandIn = str(pdbInBase + ".ligands.cif")
                    if not os.path.exists(self.ligandIn):
                        self.ligandIn = str(None)
        else:
            self.ligandIn = str(None)
        print("")

        elementsToTry = input("Which elements to try, comma separated: ")
        self.elements = [x.strip() for x in elementsToTry.split(",")]

        mtzInfo = reflection_file_reader.any_reflection_file(self.mtzIn)
        mtzInfo.file_content()
        mtzobj = iotbx.mtz.object(self.mtzIn)
        self.space_group = mtzobj.space_group_name()
        csym = mtzobj.crystals()[0].crystal_symmetry()
        unit_cell = csym.unit_cell()
        self.unit_cell_strip = (str(unit_cell)).strip("()").replace(",", "")
        self.mtz = gemmi.read_mtz_file(self.mtzIn)
        self.WV = self.mtz.dataset(1).wavelength
        self.scrapedData = []
        self.mtz = None

    def wavelength_to_eV(self):
        h = 6.62607015e-34
        c = 2.99792458e8
        eV = 1.602176634e-19
        self.energyeV = (h * c) / (float(self.WV) * 1e-10) / eV

    def lookup_fprime(self, element):
        lookupPath = os.path.join(lookup_path, f"{element}.dat")
        closestEnergy = None
        closestValues = None
        with open(lookupPath, "r") as lookupFile:
            for line in lookupFile:
                parts = line.split()
                currentEnergy, fPrime, fDoublePrime = map(float, parts)
                if closestEnergy is None or abs(self.energyeV - currentEnergy) < abs(
                    self.energyeV - closestEnergy
                ):
                    closestEnergy = currentEnergy
                    closestValues = (fPrime, fDoublePrime)
        return closestValues

    def change_elem(self):
        pdbInBase, pdbInExt = self.pdbIn.rsplit(".", 1)
        pdbInBase = os.path.basename(pdbInBase)
        print(f"Running on {str(pdbInBase)}")
        pdbLinesWrite = []
        self.pdbList = []

        with open(self.pdbIn, "r") as file:
            pdbLinesWrite = file.readlines()

        toChange = []
        for line in pdbLinesWrite:
            if line.startswith("HETATM") and not any(
                value in line for value in self.toignore
            ):
                print(line.strip())
                changeHETATM = input(
                    "Would you like to run refinement on this HETATM? "
                )
                if changeHETATM.lower() in ("y", "yes"):
                    toChange.append(line)

        for element in self.elements:
            toFDPRefine = []
            toWrite = copy.deepcopy(pdbLinesWrite)
            for i, line in enumerate(pdbLinesWrite):
                if line in toChange:
                    elementSymbol = (
                        element.rjust(2)[:2] if len(element) == 2 else " " + element
                    )
                    residueName = element.rjust(3)[:3]
                    elementName = element.ljust(4)[:4]
                    lineList = list(line)
                    lineList[12:16] = elementName
                    lineList[17:20] = residueName
                    lineList[76:78] = elementSymbol
                    toWrite[i] = "".join(lineList)
                    toFDPRefine.append(
                        [element, line[21:22].strip(), line[22:26].strip()]
                    )

            with open(f"{pdbInBase}_{element}.{pdbInExt}", "w") as pdbOut:
                for line in toWrite:
                    if not line.endswith("\n"):
                        line += "\n"
                    pdbOut.write(line)
                self.pdbList.append(
                    (f"{pdbInBase}_{element}.{pdbInExt}", element, toFDPRefine)
                )

    def scrapeLastAnomalousGroupData(self, ele, closestValues):
        log_file_path = f"{self.projIn}_fdp_{ele}_1.log"
        if not os.path.exists(log_file_path):
            print(
                f"Warning: Log file '{log_file_path}' not found for element {ele}. "
                f"phenix.refine may have failed — check {ele}_output.log for details."
            )
            return
        with open(log_file_path, "r") as file:
            content = file.read()

        last_macro_cycle_index = content.rfind("MACRO_CYCLE")

        content_after_macro_cycle = content[last_macro_cycle_index:]

        pattern = re.compile(
            r"Anomalous scatterer group:\s+"
            r'Selection: "chain (?P<chain>[A-Za-z]) and resid (?P<resid>\d+) and element (?P<element>[A-Za-z]+)"\s+'
            r"Number of selected scatterers: \d+\s+"
            r"f_prime: +[+-]?\d+\.\d+\s+"
            r"f_double_prime: (?P<f_double_prime>[+-]?\d+\.\d+)"
        )

        matches = pattern.findall(content_after_macro_cycle)
        data = [
            (
                m[0],
                int(m[1]),
                m[2],
                np.round(float(m[3]), 4),
                np.round(float(closestValues[1]), 4),
                np.round(float(abs(float(closestValues[1]) - (float(m[3])))), 4),
            )
            for m in matches
        ]

        self.scrapedData.append(data)

    def makeTable(self):
        header = [
            "Chain",
            "ResidID",
            "Element",
            "Refined FDP",
            "Theoretical FDP",
            "Absolute Difference",
            "Percentage Difference"
        ]
        flattenedData = [item for sublist in self.scrapedData for item in sublist]
        newData = []
        for row in flattenedData:
            theoretical = row[4]
            abs_diff = row[5]
            if theoretical != 0:
                perc_diff = np.round(((abs_diff - theoretical) / theoretical) * 100, 4)
            else:
                perc_diff = None
            newData.append(row + (perc_diff,))
        tableData = list(map(list, zip(*newData)))
    
        fig = go.Figure(
            data=[
                go.Table(
                    header=dict(
                        values=header, fill_color="paleturquoise", align="left"
                    ),
                    cells=dict(values=tableData, fill_color="lavender", align="left"),
                )
            ]
        )
        pio.write_html(fig, f"{self.projIn}.html")
    
    def makeTable2(self):
        header = [
            "Chain",
            "ResidID",
            "Element",
            "Refined FDP",
            "Theoretical FDP",
            "Absolute Difference",
            "Percentage Difference"
        ]
    
        flattenedData = [item for sublist in self.scrapedData for item in sublist]
        newRows = []
        for row in flattenedData:
            theoretical = row[4]
            abs_diff = row[5]
            if theoretical != 0:
                perc_diff = np.round(((abs_diff - theoretical) / theoretical) * 100, 4)
            else:
                perc_diff = None
            newRows.append(row + (perc_diff,))
    
        rows_html = ""
        for row in newRows:
            row_html = "<tr>" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>"
            rows_html += row_html
    
        html_template = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>Sortable Table</title>
            <!-- DataTables CSS -->
            <link rel="stylesheet" href="https://cdn.datatables.net/1.11.3/css/jquery.dataTables.min.css">
            <!-- jQuery -->
            <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
            <!-- DataTables JS -->
            <script src="https://cdn.datatables.net/1.11.3/js/jquery.dataTables.min.js"></script>
            <script>
                $(document).ready(function() {{
                    $('#sortableTable').DataTable();
                }});
            </script>
        </head>
        <body>
            <table id="sortableTable" class="display">
                <thead>
                    <tr>
                        {"".join(f"<th>{col}</th>" for col in header)}
                    </tr>
                </thead>
                <tbody>
                    {rows_html}
                </tbody>
            </table>
        </body>
        </html>
        """
        
        with open(f"{self.projIn}.html", "w") as f:
            f.write(html_template)

    def runBPos(self, pdbIn, elementIn):
        prefix = f"{self.projIn}_bpos_{str(elementIn)}"
        with open(f"bposEffParam_{elementIn}.eff", "w") as file:
            file.write(
                f'''refinement {{
      crystal_symmetry {{
        unit_cell = {self.unit_cell_strip}
        space_group = {self.space_group}
      }}
      pdb_interpretation {{
        allow_polymer_cross_special_position = True
        restraints_library {{
          cdl = False
        }}
      }}
      output {{
        job_title = """{self.projIn}"""
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
        strategy = *individual_sites individual_sites_real_space rigid_body \\
                  *individual_adp group_adp tls occupancies group_anomalous
      }}
      main {{
        number_of_macro_cycles = 5
        wavelength = {self.WV}
        nproc = {self.cpus}
      }}
      }}'''
            )

        cmd = [
            "phenix.refine",
            pdbIn,
            self.mtzIn,
            f"bposEffParam_{elementIn}.eff",
            f"xray_data.labels=IMEAN,SIGIMEAN",
            f"xray_data.r_free_flags.label=FreeR_flag",
            f"xray_data.r_free_flags.test_flag_value=0",
            f"output.prefix={prefix}",
            f"output.serial_format=%d",
        ]
        if self.ligandIn != str(None):
            cmd.append(self.ligandIn)
        logFile = f"{elementIn}_output.log"
        with open(logFile, "a") as log:
            result = subprocess.run(cmd, stdout=log, stderr=log)
        if result.returncode != 0:
            raise RuntimeError(
                f"phenix.refine (bpos) failed for {elementIn} with return code {result.returncode}. "
                f"Check {logFile} for details."
            )

    def runFdp(self, elementIn, toFDPRefine, closestValues):
        anomalousScatterersStr = "anomalous_scatterers {\n"
        for elemID, chainID, residID in toFDPRefine:
            groupStr = f"""          group {{
              selection = chain {chainID} and resid {residID} and element {elemID}
              f_prime = {closestValues[0]}
              refine = f_prime *f_double_prime
            }}\n"""
            anomalousScatterersStr += groupStr
        anomalousScatterersStr += "        }"

        prefix = f"{self.projIn}_fdp_{elementIn}"
        bposPdb = f"{self.projIn}_bpos_{elementIn}_1.pdb"
        with open(f"fdpEffParam_{elementIn}.eff", "w") as file:
            file.write(
                f'''refinement {{
      crystal_symmetry {{
        unit_cell = {self.unit_cell_strip}
        space_group = {self.space_group}
      }}
      pdb_interpretation {{
        allow_polymer_cross_special_position = True
        restraints_library {{
          cdl = False
        }}
      }}
      output {{
        job_title = """{self.projIn}_{elementIn}"""
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
        strategy = individual_sites individual_sites_real_space rigid_body \\
                  individual_adp group_adp tls occupancies *group_anomalous
      	{anomalousScatterersStr}
      }}
      main {{
        number_of_macro_cycles = 5
        wavelength = {self.WV}
      }}
      }}'''
            )

        cmd = [
            "phenix.refine",
            bposPdb,
            self.mtzIn,
            f"fdpEffParam_{elementIn}.eff",
            f"xray_data.labels=I(+),SIGI(+),I(-),SIGI(-)",
            f"xray_data.r_free_flags.label=FreeR_flag",
            f"xray_data.r_free_flags.test_flag_value=0",
            f"output.prefix={prefix}",
            f"output.serial_format=%d",
        ]
        if self.ligandIn != str(None):
            cmd.append(self.ligandIn)
        logFile = f"{elementIn}_output.log"
        with open(logFile, "a") as log:
            result = subprocess.run(cmd, stdout=log, stderr=log)
        if result.returncode != 0:
            raise RuntimeError(
                f"phenix.refine (fdp) failed for {elementIn} with return code {result.returncode}. "
                f"Check {logFile} for details."
            )


def runParallel(args):
    run, pdb, ele, tfdpr = args
    try:
        closestValues = run.lookup_fprime(ele)
        run.runBPos(pdbIn=pdb, elementIn=ele)
        run.runFdp(elementIn=ele, toFDPRefine=tfdpr, closestValues=closestValues)
    except RuntimeError as e:
        print(f"\nError processing {pdb} ({ele}): {e}")
    except Exception as e:
        print(f"\nUnexpected error processing {pdb} ({ele}): {e}")


if __name__ == "__main__":
    run = refinefdoubleprime()
    run.change_elem()
    run.wavelength_to_eV()
    toRun = [(run, pdb, ele, tfdpr) for pdb, ele, tfdpr in run.pdbList]
    with Halo(text="Running phenix.refine in parallel", spinner="clock"):
        with Pool(processes=os.cpu_count() - 1) as pool:
            pool.map(runParallel, toRun)
    with Halo(text="Scraping data", spinner="hearts"):
        for _, ele, _ in run.pdbList:
            closestValues = run.lookup_fprime(ele)
            run.scrapeLastAnomalousGroupData(ele=ele, closestValues=closestValues)
        run.makeTable()
