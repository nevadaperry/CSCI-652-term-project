import tomli
import os
import MafProcessing

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
        "pw" : {},
        "pw-more" : {}}

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
    data["pw"][filename.split(".sing.maf")[0]] = MafProcessing.processFullMafFile(f"{paths.data}\\pw\\{filename}")
for filename in os.listdir(f"{paths.data}\\pw-more"):
    data["pw-more"][filename.split(".sing.maf")[0]] = MafProcessing.processFullMafFile(f"{paths.data}\\pw-more\\{filename}")

#endregion === Data Setup ===
