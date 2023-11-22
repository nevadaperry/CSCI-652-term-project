"""
Author: Alex Somheil
Student ID: Z1888439
Date: 11/18/2023

This file contains various functions for identifying and classifying variants between two
genomes that have been pairwise aligned and processed by MafProcessing.py.
"""




# This method accepts a full processed pairwise alignment dictionary, and attempts to find positions
# at which variations occur.
def indexPotentialVariants(pwDict : dict):

    genome1 = pwDict["Genomes"][0]
    genome2 = pwDict["Genomes"][1]

    foundVariants = []
    for block in pwDict["Blocks"]:
        blockSize = len(block[f"{genome1}-AlignmentSequence"])

        genome1_startIndex = block[f"{genome1}-StartPosition"]
        genome2_startIndex = block[f"{genome2}-StartPosition"]

        lastStartIndex = 0
        readingVariant = False
        runningGenome1SubSequence = ""
        runningGenome1GapCount = 0
        runningGenome2SubSequence = ""
        runningGenome2GapCount = 0
        isGapVariant = False
        for index in range(blockSize):
            genome1_Value = block[f"{genome1}-AlignmentSequence"][index]
            genome2_Value = block[f"{genome2}-AlignmentSequence"][index]

            # We make sure to count gaps to get the correct position of the variants.
            if(genome1_Value == "-"):
                runningGenome1GapCount += 1
            if(genome2_Value == "-"):
                runningGenome2GapCount += 1

            if(genome1_Value == genome2_Value):
                if(readingVariant):
                    readingVariant = False
                    foundVariants.append({"Length" : len(runningGenome1SubSequence),
                                          f"{genome1}-Location" : (lastStartIndex + genome1_startIndex - runningGenome1GapCount,index + genome1_startIndex - runningGenome1GapCount), f"{genome2}-Location" : (lastStartIndex + genome2_startIndex - runningGenome2GapCount, index + genome2_startIndex - runningGenome2GapCount),
                                          f"{genome1}-SubSequence" : runningGenome1SubSequence,f"{genome2}-SubSequence" : runningGenome2SubSequence})
                    runningGenome1SubSequence = ""
                    runningGenome2SubSequence = ""
            else:
                # We have to test to make sure this is still the same type of mutation and that it shouldn't
                # be treated as two separate variants.
                if(genome1_Value.isalnum() and genome2_Value.isalnum()):
                    # This special case means a polymorphism is right next to a gap, but we still must treat them
                    # separately.
                    if(readingVariant and isGapVariant):
                        foundVariants.append({"Length": len(runningGenome1SubSequence),
                                              f"{genome1}-Location": (lastStartIndex + genome1_startIndex - runningGenome1GapCount,index + genome1_startIndex - runningGenome1GapCount),f"{genome2}-Location": (lastStartIndex + genome2_startIndex - runningGenome2GapCount,index + genome2_startIndex - runningGenome2GapCount),
                                              f"{genome1}-SubSequence": runningGenome1SubSequence,f"{genome2}-SubSequence": runningGenome2SubSequence})
                        runningGenome1SubSequence = ""
                        runningGenome2SubSequence = ""
                        lastStartIndex = index
                    isGapVariant = False
                else:
                    # This special case means a polymorphism is right next to a gap, but we still must treat them
                    # separately.
                    if(readingVariant and not isGapVariant):
                        foundVariants.append({"Length": len(runningGenome1SubSequence),
                                              f"{genome1}-Location": (lastStartIndex + genome1_startIndex - runningGenome1GapCount,index + genome1_startIndex - runningGenome1GapCount),f"{genome2}-Location": (lastStartIndex + genome2_startIndex - runningGenome2GapCount,index + genome2_startIndex - runningGenome2GapCount),
                                              f"{genome1}-SubSequence": runningGenome1SubSequence,f"{genome2}-SubSequence": runningGenome2SubSequence})
                        runningGenome1SubSequence = ""
                        runningGenome2SubSequence = ""
                        lastStartIndex = index
                    isGapVariant = True

                if(not readingVariant):
                    readingVariant = True
                    lastStartIndex = index
                runningGenome1SubSequence += genome1_Value
                runningGenome2SubSequence += genome2_Value
        if (readingVariant):
            foundVariants.append(
                {f"{genome1}-Location": (lastStartIndex + genome1_startIndex - runningGenome1GapCount, (blockSize - 1) + genome1_startIndex - runningGenome1GapCount),
                 f"{genome2}-Location": (lastStartIndex + genome2_startIndex - runningGenome2GapCount, (blockSize - 1) + genome2_startIndex - runningGenome2GapCount), })

    pwDict["Variants"] = foundVariants

# This method accepts a single processed pairwise alignment dictionary, and attempts to classify what TYPE of variant
# it is.
def classifyVariants(pwDict : dict):
    genome1Name = pwDict["Genomes"][0]
    genome2Name = pwDict["Genomes"][1]
    for variantDict in pwDict["Variants"]:
        if(variantDict[f"{genome1Name}-SubSequence"].isalnum()):
            if(variantDict[f"{genome2Name}-SubSequence"].isalnum()):
                variantDict["Type"] = "Polymorphism"
            elif(all(c == '-' for c in variantDict[f"{genome2Name}-SubSequence"])):
                variantDict["Type"] = "Deletion"
            else:
                raise ValueError(f"ERROR: Mixed variant detected!\n{genome1Name}: {variantDict[f'{genome1Name}-SubSequence']}\n{genome2Name}: {variantDict[f'{genome2Name}-SubSequence']}")
        elif(all(c == '-' for c in variantDict[f"{genome1Name}-SubSequence"])):
            if(variantDict[f"{genome2Name}-SubSequence"].isalnum()):
                variantDict["Type"] = "Insertion"
            else:
                raise ValueError(f"ERROR: Invalid insertion type detected! \n{genome1Name}: {variantDict[f'{genome1Name}-SubSequence']}\n{genome2Name}: {variantDict[f'{genome2Name}-SubSequence']}")




