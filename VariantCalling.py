"""
Author: Alex Somheil
Student ID: Z1888439
Date: 11/18/2023

This file contains various functions for identifying and classifying mutations between two
genomes that have been pairwise aligned and processed by MafProcessing.py.
"""




# This method accepts a full processed pairwise alignment dictionary, and attempts to find positions
# at which mutations occur.
def indexPotentialMutations(pwDict : dict):

    genome1 = pwDict["Genomes"][0]
    genome2 = pwDict["Genomes"][1]

    foundMutations = []
    for block in pwDict["Blocks"]:
        blockSize = len(block[f"{genome1}-AlignmentSequence"])

        genome1_startIndex = block[f"{genome1}-StartPosition"]
        genome2_startIndex = block[f"{genome2}-StartPosition"]

        lastStartIndex = 0
        readingMutation = False
        runningGenome1SubSequence = ""
        runningGenome1GapCount = 0
        runningGenome2SubSequence = ""
        runningGenome2GapCount = 0
        isGapMutation = False
        for index in range(blockSize):
            genome1_Value = block[f"{genome1}-AlignmentSequence"][index]
            genome2_Value = block[f"{genome2}-AlignmentSequence"][index]

            # We make sure to count gaps to get the correct position of the mutations.
            if(genome1_Value == "-"):
                runningGenome1GapCount += 1
            if(genome2_Value == "-"):
                runningGenome2GapCount += 1

            if(genome1_Value == genome2_Value):
                if(readingMutation):
                    readingMutation = False
                    foundMutations.append({"Length" : len(runningGenome1SubSequence), "Genome1" : genome1, "Genome2" : genome2,
                                          f"{genome1}-Location" : (lastStartIndex + genome1_startIndex - runningGenome1GapCount,index + genome1_startIndex - runningGenome1GapCount), f"{genome2}-Location" : (lastStartIndex + genome2_startIndex - runningGenome2GapCount, index + genome2_startIndex - runningGenome2GapCount),
                                          f"{genome1}-SubSequence" : runningGenome1SubSequence,f"{genome2}-SubSequence" : runningGenome2SubSequence})
                    runningGenome1SubSequence = ""
                    runningGenome2SubSequence = ""
            else:
                # We have to test to make sure this is still the same type of mutation and that it shouldn't
                # be treated as two separate mutations.
                if(genome1_Value.isalnum() and genome2_Value.isalnum()):
                    # This special case means a polymorphism is right next to a gap, but we still must treat them
                    # separately.
                    if(readingMutation and isGapMutation):
                        foundMutations.append({"Length": len(runningGenome1SubSequence), "Genome1" : genome1, "Genome2" : genome2,
                                              f"{genome1}-Location": (lastStartIndex + genome1_startIndex - runningGenome1GapCount,index + genome1_startIndex - runningGenome1GapCount),f"{genome2}-Location": (lastStartIndex + genome2_startIndex - runningGenome2GapCount,index + genome2_startIndex - runningGenome2GapCount),
                                              f"{genome1}-SubSequence": runningGenome1SubSequence,f"{genome2}-SubSequence": runningGenome2SubSequence})
                        runningGenome1SubSequence = ""
                        runningGenome2SubSequence = ""
                        lastStartIndex = index
                    isGapMutation = False
                else:
                    # This special case means a polymorphism is right next to a gap, but we still must treat them
                    # separately.
                    if(readingMutation and not isGapMutation):
                        foundMutations.append({"Length": len(runningGenome1SubSequence), "Genome1" : genome1, "Genome2" : genome2,
                                              f"{genome1}-Location": (lastStartIndex + genome1_startIndex - runningGenome1GapCount,index + genome1_startIndex - runningGenome1GapCount),f"{genome2}-Location": (lastStartIndex + genome2_startIndex - runningGenome2GapCount,index + genome2_startIndex - runningGenome2GapCount),
                                              f"{genome1}-SubSequence": runningGenome1SubSequence,f"{genome2}-SubSequence": runningGenome2SubSequence})
                        runningGenome1SubSequence = ""
                        runningGenome2SubSequence = ""
                        lastStartIndex = index
                    isGapMutation = True

                if(not readingMutation):
                    readingMutation = True
                    lastStartIndex = index
                runningGenome1SubSequence += genome1_Value
                runningGenome2SubSequence += genome2_Value
        if (readingMutation):
            foundMutations.append(
                {"Length": len(runningGenome1SubSequence), "Genome1" : genome1, "Genome2" : genome2,
                 f"{genome1}-Location": (lastStartIndex + genome1_startIndex - runningGenome1GapCount, (blockSize - 1) + genome1_startIndex - runningGenome1GapCount),
                 f"{genome2}-Location": (lastStartIndex + genome2_startIndex - runningGenome2GapCount, (blockSize - 1) + genome2_startIndex - runningGenome2GapCount)})

    pwDict["Mutations"] = foundMutations

# This method accepts a single processed pairwise alignment dictionary, and attempts to classify what TYPE of mutation
# it is.
def classifyMutationTypes(pwDict : dict):
    genome1Name = pwDict["Genomes"][0]
    genome2Name = pwDict["Genomes"][1]

    # Identify mutation type.
    for mutationDict in pwDict["Mutations"]:
        genome1SubSeq = mutationDict[f"{genome1Name}-SubSequence"]
        genome2SubSeq = mutationDict[f"{genome2Name}-SubSequence"]
        if(genome1SubSeq.isalnum()):
            if(genome2SubSeq.isalnum()):
                if(len(genome1SubSeq) == 1):
                    if(genome1SubSeq + genome2SubSeq in ("AG", "GA", "CT", "TC")):
                        mutationDict["Type"] = "Transition"
                    else:
                        mutationDict["Type"] = "Transversion"
                else:
                    mutationDict["Type"] = "MNP"
            elif(all(c == '-' for c in genome2SubSeq)):
                mutationDict["Type"] = "Deletion"
            else:
                raise ValueError(f"ERROR: Mixed mutation detected!\n{genome1Name}: {genome1SubSeq}\n{genome2Name}: {genome2SubSeq}")
        elif(all(c == '-' for c in genome1SubSeq)):
            if(genome2SubSeq.isalnum()):
                mutationDict["Type"] = "Insertion"
            else:
                raise ValueError(f"ERROR: Invalid insertion type detected! \n{genome1Name}: {genome1SubSeq}\n{genome2Name}: {genome2SubSeq}")


# This method accepts a single processed pairwise alignment dictionary, and attempts to classify the location
# it occurs on the SARS-COV2 genome, given a regionDict.
def classifyMutationLocations(pwDict : dict,regionDict : dict):
    genome1Name = pwDict["Genomes"][0]

    for mutationDict in pwDict["Mutations"]:
        genome1StartLocation = mutationDict[f"{genome1Name}-Location"][0]

        # First we find the main gene (if applicable)
        foundLocation = None
        for geneName,geneContents in regionDict.items():
            if(geneContents["Location"][0] < genome1StartLocation < geneContents["Location"][1]):
                foundLocation = geneName
                break
        mutationDict["Unit"] = foundLocation if foundLocation is not None else "Other"

        # Now we check if it exists on a subunit of its main gene. Since subunit locations are in reference only
        # to the genome, we have to account for this.
        if(mutationDict["Unit"] != "Other"):
            if(regionDict[mutationDict["Unit"]]["SubUnits"] is not None):
                foundLocation = None
                geneStartOffset = regionDict[mutationDict["Unit"]]["Location"][0]
                for subUnitName,subUnitContents in regionDict[mutationDict["Unit"]]["SubUnits"].items():
                    if((subUnitContents["Location"][0] + geneStartOffset) < genome1StartLocation < (subUnitContents["Location"][1] + geneStartOffset)):
                        foundLocation = subUnitName
                        break
                mutationDict["SubUnit"] = foundLocation if foundLocation is not None else "Other"
            else:
                mutationDict["SubUnit"] = "Other"
        else:
            mutationDict["SubUnit"] = "Other"

