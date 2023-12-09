"""
Author: Alex Somheil
Student ID: Z1888439
Date: 11/20/2023

This file contains functions used to help visualize some of the higher level functional
consequences resulted from comparing SARS-COV2 to its variants.
"""
import Setup as s
maxMutationSize = 10**100

data = s.data
transmissibility = s.transmissibility

# Given a list of mutation dicts, this method returns a similar list ordered
# by size. Optional parameter to omit mutations smaller than a certain size.
def filterMutationsBySize(mutationList : list,minSize : int = 0,maxSize : int = maxMutationSize):
    sortedMutations = [item for item in sorted(mutationList, key=lambda x: x['Length'], reverse=True) if minSize <= item['Length'] <= maxSize]
    return sortedMutations

# Given a list of mutation dicts, this method returns a similar list filtered only for
# a specific unit (and option subunit). Optional booleans to set unit/subunit as blacklist instead of
# whitelist.
def filterMutationsByUnit(mutationList : list, unit : str | list = None, subunit : str | list = None, unitBlackList : bool = False, subunitBlackList : bool = False):
    # Convert unit and subunit to lists if they are not already.
    if(type(unit) is not list):
        unit = [unit] if unit is not None else []
    if(type(subunit) is not list):
        subunit = [subunit] if subunit is not None else []

    filteredMutations = mutationList
    # Filtering logic for units.
    if(unit):
        newFilteredList = []
        for mutation in filteredMutations:
            for thisUnitFilter in unit:
                if(thisUnitFilter == mutation["Unit"]):
                    if(unitBlackList):
                        continue
                    else:
                        newFilteredList.append(mutation)
        filteredMutations = newFilteredList

    # Filtering logic for subunits.
    if(subunit):
        newFilteredList = []
        for mutation in filteredMutations:
            for thisSubunitFilter in subunit:
                if(thisSubunitFilter == mutation["SubUnit"]):
                    if(subunitBlackList):
                        break
                    else:
                        newFilteredList.append(mutation)
                        break
        filteredMutations = newFilteredList

    return sorted(filteredMutations, key=lambda x: x['Unit'], reverse=True)

# Given a list of mutation dicts, this method returns a similar list filtered for mutations that occur between a
# specific minIndex and maxIndex
def filterMutationsByLocation(mutationList : list, minIndex : int = 0, maxIndex : int = maxMutationSize):
    sortedMutations = [item for item in sorted(mutationList, key=lambda x: x[f"{x['Genome1']}-Location"][0], reverse=True) if minIndex <= item[f"{item['Genome1']}-Location"][0] <= maxIndex]
    return sortedMutations

# Given a list of mutation dicts, this method returns a similar list filtered for mutations that occur only for a
# specific genome2 (compared genome)
def filterMutationsByComparedGenome(mutationList : list, compareGenomeName : str):
    sortedMutations = [item for item in mutationList if item['Genome2'].lower() == compareGenomeName.lower()]
    return sortedMutations

# Given a list of mutation dicts, this method returns a similar list filtered for mutations of a certain type of indel
# (insertion, deletion, transversion, etc.)
def filterMutationsByType(mutationList : list, mutationType : str | list):
    if(type(mutationType) is str):
        mutationType = [mutationType]
    sortedMutations = [item for item in mutationList if item['Type'] in mutationType]
    return sortedMutations


# Concatenating all found mutations together into one list for ease of use.
mutations = []
for _pwDictName,_pwDict in s.data["pw"].items():
    mutations += _pwDict["Mutations"]



# Filtering our mutation list to only view mutations that occur on the Spike gene, and within a specific subunit.
spikeMutations = filterMutationsByUnit(mutationList = mutations, unit = "S",)

# Creating lists of mutations that occur on each specific subUnit of the Spike gene.
sub_spikeMutations_NTD = filterMutationsByUnit(mutationList=spikeMutations,subunit="NTD")
sub_spikeMutations_RBD = filterMutationsByUnit(mutationList=spikeMutations,subunit="RBD")
sub_spikeMutations_Cleavage = filterMutationsByUnit(mutationList=spikeMutations,subunit="Cleavage")
sub_spikeMutations_FP = filterMutationsByUnit(mutationList=spikeMutations,subunit="FP")
sub_spikeMutations_IFP = filterMutationsByUnit(mutationList=spikeMutations,subunit="IFP")
sub_spikeMutations_HR1 = filterMutationsByUnit(mutationList=spikeMutations,subunit="HR1")
sub_spikeMutations_HR2 = filterMutationsByUnit(mutationList=spikeMutations,subunit="HR2")
sub_spikeMutations_TM = filterMutationsByUnit(mutationList=spikeMutations,subunit="TM")
sub_spikeMutations_CT = filterMutationsByUnit(mutationList=spikeMutations,subunit="CT")
sub_spikeMutations_Other = filterMutationsByUnit(mutationList=spikeMutations,subunit="Other")

all_mutations_alpha = filterMutationsByComparedGenome(mutationList=mutations,compareGenomeName="alpha")
all_mutations_beta = filterMutationsByComparedGenome(mutationList=mutations,compareGenomeName="beta")
all_mutations_delta = filterMutationsByComparedGenome(mutationList=mutations,compareGenomeName="delta")
all_mutations_gamma = filterMutationsByComparedGenome(mutationList=mutations,compareGenomeName="gamma")
all_mutations_omicronBA1 = filterMutationsByComparedGenome(mutationList=mutations,compareGenomeName="omicronBA1")

all_spikeMutations_alpha = filterMutationsByComparedGenome(mutationList=spikeMutations,compareGenomeName="alpha")
all_spikeMutations_beta = filterMutationsByComparedGenome(mutationList=spikeMutations,compareGenomeName="beta")
all_spikeMutations_delta = filterMutationsByComparedGenome(mutationList=spikeMutations,compareGenomeName="delta")
all_spikeMutations_gamma = filterMutationsByComparedGenome(mutationList=spikeMutations,compareGenomeName="gamma")
all_spikeMutations_omicronBA1 = filterMutationsByComparedGenome(mutationList=spikeMutations,compareGenomeName="omicronBA1")

'''
all_spikeInsertions_alpha = filterMutationsByType(mutationList=all_spikeMutations_alpha,mutationType="Insertion")
all_spikeDeletions_alpha = filterMutationsByType(mutationList=all_spikeMutations_alpha,mutationType="Deletion")
all_spikeSNPs_alpha = filterMutationsByType(mutationList=all_spikeMutations_alpha,mutationType=["Transversion","Transition"])
all_spikeMNPs_alpha = filterMutationsByType(mutationList=all_spikeMutations_alpha,mutationType="MNP")

all_spikeInsertions_beta = filterMutationsByType(mutationList=all_spikeMutations_beta, mutationType="Insertion")
all_spikeDeletions_beta = filterMutationsByType(mutationList=all_spikeMutations_beta, mutationType="Deletion")
all_spikeSNPs_beta = filterMutationsByType(mutationList=all_spikeMutations_beta, mutationType=["Transversion", "Transition"])
all_spikeMNPs_beta = filterMutationsByType(mutationList=all_spikeMutations_beta, mutationType="MNP")

all_spikeInsertions_delta = filterMutationsByType(mutationList=all_spikeMutations_delta, mutationType="Insertion")
all_spikeDeletions_delta = filterMutationsByType(mutationList=all_spikeMutations_delta, mutationType="Deletion")
all_spikeSNPs_delta = filterMutationsByType(mutationList=all_spikeMutations_delta, mutationType=["Transversion", "Transition"])
all_spikeMNPs_delta = filterMutationsByType(mutationList=all_spikeMutations_delta, mutationType="MNP")

all_spikeInsertions_gamma = filterMutationsByType(mutationList=all_spikeMutations_gamma, mutationType="Insertion")
all_spikeDeletions_gamma = filterMutationsByType(mutationList=all_spikeMutations_gamma, mutationType="Deletion")
all_spikeSNPs_gamma = filterMutationsByType(mutationList=all_spikeMutations_gamma, mutationType=["Transversion", "Transition"])
all_spikeMNPs_gamma = filterMutationsByType(mutationList=all_spikeMutations_gamma, mutationType="MNP")

all_spikeInsertions_omicronBA1 = filterMutationsByType(mutationList=all_spikeMutations_omicronBA1, mutationType="Insertion")
all_spikeDeletions_omicronBA1 = filterMutationsByType(mutationList=all_spikeMutations_omicronBA1, mutationType="Deletion")
all_spikeSNPs_omicronBA1 = filterMutationsByType(mutationList=all_spikeMutations_omicronBA1, mutationType=["Transversion", "Transition"])
all_spikeMNPs_omicronBA1 = filterMutationsByType(mutationList=all_spikeMutations_omicronBA1, mutationType="MNP")
'''
