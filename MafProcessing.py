import re

"""
Author: Alex Somheil
Student ID: Z1888439
Date: 11/16/2023

This file contains various functions for parsing .maf files containing single pairwise alignments
between two genome alignment sequences into neatly formatted, adaptable Python dictionaries. It
also helps to calculate some basic information about the sequences such as gap rate, substitution rate,
and more.
"""




# This method converts a maf file into an organized, parsed dictionary for
# use with processing data.
def readMafFile(mafFilePath):
    mafFile = open(mafFilePath, 'r')

    comparedGenomes = []
    blockList = []
    currentBlock = {}
    for line in mafFile:
        line = line.strip()
        # Skip if a blank line
        if(line == ""):
            continue
        # Skip if operator line
        elif(line.startswith("##")):
            continue
        # Signifies the beginning of a new block.
        elif(line.startswith("a")):
            if(currentBlock != {}):
                blockList.append(currentBlock)
            currentBlock = {"Score": line.split("score=")[1]}
        elif(line.startswith("s")):
            sRowContents = re.split(r"\s+",line)
            genomeName = sRowContents[1].split(".")[0]
            if(genomeName not in comparedGenomes):
                comparedGenomes.append(genomeName)
            currentBlock[f"{genomeName}-Chromosome"] = sRowContents[1]
            currentBlock[f"{genomeName}-StartPosition"] = int(sRowContents[2])
            currentBlock[f"{genomeName}-SequenceSize"] = int(sRowContents[3])
            currentBlock[f"{genomeName}-StrandInfo"] = sRowContents[4]
            currentBlock[f"{genomeName}-ChromosomeSize"] = int(sRowContents[5])
            currentBlock[f"{genomeName}-AlignmentSequence"] = sRowContents[6]

    if (currentBlock != {}):
        blockList.append(currentBlock)
    return {"Genomes" : tuple(comparedGenomes), "Blocks" : blockList}

# This method simply classifies whether a base pair is a match, a transition, a transversion, of
# nothing.
def classifySingleBasePair(base1, base2):
    if base1 == base2:
        return "Match"
    elif base1 + base2 in ("AG", "GA", "CT", "TC"):
        return "Transition"
    elif base1 in "ACGT" and base2 in "ACGT":
        return "Transversion"
    else:
        return "None"

# This method processes the raw count of matches, transitions, transversions, and nones for a full
# single maf block based on the genomes specified to compare.
def classifyMafBlock(genomes : tuple,mafBlockDict : dict):
    # First get the total length of all target alignments.
    alignmentLength = len(mafBlockDict[f"{genomes[0]}-AlignmentSequence"])

    # Set classified values of the mafBlockDict to 0
    mafBlockDict["Match"] = 0
    mafBlockDict["Transition"] = 0
    mafBlockDict["Transversion"] = 0
    mafBlockDict["Gap"] = 0
    mafBlockDict["GapLengths"] = {}
    mafBlockDict["None"] = 0
    mafBlockDict["A"] = {"A" : 0, "G" : 0, "C" : 0, "T" : 0}
    mafBlockDict["G"] = {"A" : 0, "G" : 0, "C" : 0, "T" : 0}
    mafBlockDict["C"] = {"A" : 0, "G" : 0, "C" : 0, "T" : 0}
    mafBlockDict["T"] = {"A" : 0, "G" : 0, "C" : 0, "T" : 0}

    # Now, iterate through both alignments, comparing each column base.
    gap1 = 0
    gap2 = 0
    for i in range(alignmentLength):
        genome1Base = mafBlockDict[f"{genomes[0]}-AlignmentSequence"][i]
        genome2Base = mafBlockDict[f"{genomes[1]}-AlignmentSequence"][i]
        mafBlockDict[classifySingleBasePair(genome1Base,genome2Base)] += 1
        if(genome1Base in "ACGT"):
            if(gap1 > 0):
                mafBlockDict["Gap"] += 1
                if(mafBlockDict["GapLengths"].get(gap1) is None):
                    mafBlockDict["GapLengths"][gap1] = 0
                mafBlockDict["GapLengths"][gap1] += 1
                gap1 = 0
            regularBase1 = True
        elif(genome1Base == "-"):
            gap1 += 1
            regularBase1 = False
        else:
            if(gap1 > 0):
                mafBlockDict["Gap"] += 1
                if(mafBlockDict["GapLengths"].get(gap1) is None):
                    mafBlockDict["GapLengths"][gap1] = 0
                mafBlockDict["GapLengths"][gap1] += 1
                gap1 = 0
            regularBase1 = False

        if(genome2Base in "ACGT"):
            if(gap2 > 0):
                mafBlockDict["Gap"] += 1
                if(mafBlockDict["GapLengths"].get(gap2) is None):
                    mafBlockDict["GapLengths"][gap2] = 0
                mafBlockDict["GapLengths"][gap2] += 1
                gap2 = 0
            regularBase2 = True
        elif(genome2Base == "-"):
            gap2 += 1
            regularBase2 = False
        else:
            if (gap2 > 0):
                mafBlockDict["Gap"] += 1
                if(mafBlockDict["GapLengths"].get(gap2) is None):
                    mafBlockDict["GapLengths"][gap2] = 0
                mafBlockDict["GapLengths"][gap2] += 1
                gap2 = 0
            regularBase2 = False

        if(regularBase1 and regularBase2):
            mafBlockDict[genome1Base][genome2Base] += 1

# This method simply calculate the total values of an entire mafFile, given a classifiedMafDict.
def calculateTotalValues(mafFileDict : dict):
    returnDict = {"SubstitutionRate" : 0,
                  "TiTv" : 0,
                  "Gap" : 0,
                  "GapRate" : 0,
                  "GapLengths" : {},
                  "AlignmentCounts" : {"A" : {"A" : 0, "G" : 0, "C" : 0, "T" : 0},
                                       "G" : {"A" : 0, "G" : 0, "C" : 0, "T" : 0},
                                       "C" : {"A" : 0, "G" : 0, "C" : 0, "T" : 0},
                                       "T" : {"A" : 0, "G" : 0, "C" : 0, "T" : 0}}}

    totalMatches = 0
    totalTransitions = 0
    totalTransversions = 0
    totalGaps = 0
    for thisBlock in mafFileDict["Blocks"]:
        totalMatches += thisBlock["Match"]
        totalTransitions += thisBlock["Transition"]
        totalTransversions += thisBlock["Transversion"]
        totalGaps += thisBlock["Gap"]
        for key,value in thisBlock["GapLengths"].items():
            if(returnDict["GapLengths"].get(key) is None):
                returnDict["GapLengths"][key] = 0
            returnDict["GapLengths"][key] += value

        regularBases = ["A","G","C","T"]
        for base1 in regularBases:
            for base2 in regularBases:
                returnDict["AlignmentCounts"][base1][base2] += thisBlock[base1][base2]


    returnDict["Gap"] = totalGaps
    returnDict["SubstitutionRate"] = (totalTransitions + totalTransversions) / (totalMatches + totalTransitions + totalTransversions)
    returnDict["TiTv"] = totalTransitions / totalTransversions
    returnDict["GapRate"] = totalGaps / (totalMatches + totalTransitions + totalTransversions + totalGaps)

    return returnDict

# This method processed a full maf file, converting it from a raw file into a mafDataDict and returning
# the totalValues result dict.
def processFullMafFile(mafFilePath):
    # Read and process maf data into dictionary
    mafData = readMafFile(mafFilePath)

    # Classify the bases for each column for each block in the maf data.
    for block in mafData["Blocks"]:
        classifyMafBlock(mafData["Genomes"],block)


    finalMafData = mafData
    # noinspection PyTypeChecker
    finalMafData.update(calculateTotalValues(mafData))
    return finalMafData
