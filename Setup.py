"""
Author: Alex Somheil
Student ID: Z1888439
Date: 11/18/2023

This file contains various functions for facilitating setup of all calculated data. It
also provides a method to print out the gap rate and substitution rate of found data.
"""

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

data = {"pw" : {}}


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

# Reading all maf files. useExtendedData option for parsing ALL sars-cov2 pw-alignments, instead of just the
# 5 chosen Variants of Concern.
useExtendedData = True
if(useExtendedData):
    for filename in os.listdir(f"{paths.data}\\pw-extended"):
        if(filename.startswith("sars2")):
            data["pw"][filename.split(".sing.maf")[0]] = MafProcessing.processFullMafFile(f"{paths.data}\\pw-extended\\{filename}")
else:
    for filename in os.listdir(f"{paths.data}\\pw-voc"):
        if(filename.startswith("sars2")):
            data["pw"][filename.split(".sing.maf")[0]] = MafProcessing.processFullMafFile(f"{paths.data}\\pw-voc\\{filename}")


# Reading all transmissibility files.
transmissibility = {}
def parseTransmissibilityTable(table):
    # Split the table into lines
    lines = table.strip().split('\n')

    # Initialize an empty dictionary to store the data
    tableDict = {}

    # Iterate over each line in the table
    for _line in lines:
        # Split each line into its components
        parts = _line.split()
        country = parts[0]
        r_value = float(parts[1])
        ci_lower = float(parts[2])
        ci_upper = float(parts[3])
        r_squared = float(parts[4])
        growth_rate = float(parts[5])
        growth_rate_ci_lower = float(parts[6])
        growth_rate_ci_upper = float(parts[7])

        # Add the data to the dictionary
        tableDict[country] = {
            'R': r_value,
            'CI Lower': ci_lower,
            'CI Upper': ci_upper,
            'R Squared': r_squared,
            'Growth Rate': growth_rate,
            'Growth Rate CI Lower': growth_rate_ci_lower,
            'Growth Rate CI Upper': growth_rate_ci_upper
        }
    return tableDict
for filename in os.listdir(f"{paths.data}\\transmissibility"):
    genomeName = filename.split("_")[0]
    transmissibility[genomeName] = {}

    with open(f"{paths.data}\\transmissibility\\{filename}", "r") as f:
        parsedTable = parseTransmissibilityTable(f.read())
    transmissibility[genomeName]["Raw"] = parsedTable
# Further summarizing transmissibility data.
for genome,tmisDict in transmissibility.items():

    totalR = 0.0
    totalGrowthRate = 0.0
    for values in tmisDict["Raw"].values():
        totalR += values["R"]
        totalGrowthRate += values["Growth Rate"]

    tmisDict["AvgR"] = round(totalR / len(tmisDict["Raw"]),4)
    tmisDict["AvgGrowthRate"] = round(totalGrowthRate / len(tmisDict["Raw"]),4)


# Perform variant calling on all pairwise alignments.
for value in data["pw"].values():
    VariantCalling.indexPotentialMutations(value)
    VariantCalling.classifyMutationTypes(value)
    VariantCalling.classifyMutationLocations(value, geneAnnotations)




#endregion === Data Setup ===

# Optional section used to print calculated Gap Rate and Substitution Rate.
printGapSubRateData = False
if(printGapSubRateData):

    print("==== GAP RATES ====\n")
    for key,value in data["pw"].items():
        print(f"{key}: {value['GapRate']}")

    print("\n==== SUBSTITUTION RATES ====\n")
    for key, value in data["pw"].items():
        print(f"{key}: {value['SubstitutionRate']}")
