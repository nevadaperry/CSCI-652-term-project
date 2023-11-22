import os
import MafProcessing
import VariantCalling

# region === Pathing Setup ===

# Simple class to validate and store program paths.
class Paths:
    def __init__(self):
        # Helper function to create path (if needed) and test for access.
        def createAndCheckAccess(pathString):
            os.makedirs(pathString, exist_ok=True)
            if(not os.access(pathString,os.R_OK | os.W_OK | os.X_OK)):
                raise PermissionError(f"Insufficient permissions for configured directory: {pathString}")

        # Root path for working program directory.
        self.base = os.getcwd()
        createAndCheckAccess(self.base)

        # Data path for all data used in the program.
        self.data = f"{self.base}\\data"
        createAndCheckAccess(self.base)

paths = Paths()

# endregion === Pathing Setup ===

#region === Data Setup ===

data = {"seqs" : {},
        "pw" : {}}


# Helper method to read genome structure files that include a unit (genome or subunit) name,
# start location, and end location.
# TODO proper pathing to where geneAnnotations files are located in file structure
def readStructureFile(fileName):
    units = {}
    with open(fileName, "r") as _file:
        for _line in _file:
            if(_line.strip() != ""):
                _parts = _line.strip().split()
                units[_parts[0]] = {"Location" : (int(_parts[1]), int(_parts[2]))}
    return units

# Store information about genome regions and subregions of the sarsCov2 genome.
geneAnnotations = readStructureFile(f"{paths.data}\\sarsCov2structure.txt")
for geneName in geneAnnotations.keys():
    if(os.path.exists(f"{paths.data}\\{geneName}geneSubUnits.txt")):
        geneAnnotations[geneName]["SubUnits"] = readStructureFile(f"{paths.data}\\{geneName}geneSubUnits.txt")
    else:
        geneAnnotations[geneName]["SubUnits"] = None


# Reading all seq files.
for filename in os.listdir(f"{paths.data}\\seqs"):
    with open(f"{paths.data}\\seqs\\{filename}", "r") as f:
        seq = ""
        for line in f:
            if(line.startswith(">")):
                # Skip header/comment lines
                continue
            else:
                seq += line.strip()
        data["seqs"][filename] = seq

# Reading all maf files.
for filename in os.listdir(f"{paths.data}\\pw"):
    if(filename.startswith("sars2")):
        data["pw"][filename.split(".sing.maf")[0]] = MafProcessing.processFullMafFile(f"{paths.data}\\pw\\{filename}")


# Perform variant calling on all pairwise alignments.
for value in data["pw"].values():
    VariantCalling.indexPotentialVariants(value)
    VariantCalling.classifyVariantTypes(value)
    VariantCalling.classifyVariantLocations(value, geneAnnotations)




#endregion === Data Setup ===
