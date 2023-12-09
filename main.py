"""
Author: Alex Somheil
Student ID: Z1888439
Date: 11/18/2023

This file contains various matplotlib graphs and charts for displaying valuable
information found by previous data processing.
"""

import FunctionalConsequences as fc
import matplotlib.pyplot as plt

# Simply display gap rate and sub rate for all maf processed files.
'''
for key,value in b.data["pw"].items():
    print(f"=== {key} ===")
    print(f"Gap Rate: {value['GapRate']}")
    print(f"Sub Rate: {value['SubstitutionRate']}\n")
'''

# Display a bar graph of the frequency of spike protein variations between a certain size range
# across all pairwise alignments
def spikeMutationsFrequencyBySubunit():
    minimumSize = 0
    maximumSize = 1

    #region === Processing ===
    adjusted_spikeMutations_NTD = fc.filterMutationsBySize(variantList = fc.sub_spikeMutations_NTD,minVariantSize=minimumSize,maxVariantSize=maximumSize)
    adjusted_spikeMutations_RBD = fc.filterMutationsBySize(variantList = fc.sub_spikeMutations_RBD,minVariantSize=minimumSize,maxVariantSize=maximumSize)
    adjusted_spikeMutations_Cleavage = fc.filterMutationsBySize(variantList = fc.sub_spikeMutations_Cleavage,minVariantSize=minimumSize,maxVariantSize=maximumSize)
    adjusted_spikeMutations_FP = fc.filterMutationsBySize(variantList = fc.sub_spikeMutations_FP,minVariantSize=minimumSize,maxVariantSize=maximumSize)
    adjusted_spikeMutations_IFP = fc.filterMutationsBySize(variantList = fc.sub_spikeMutations_IFP,minVariantSize=minimumSize,maxVariantSize=maximumSize)
    adjusted_spikeMutations_HR1 = fc.filterMutationsBySize(variantList = fc.sub_spikeMutations_HR1,minVariantSize=minimumSize,maxVariantSize=maximumSize)
    adjusted_spikeMutations_HR2 = fc.filterMutationsBySize(variantList = fc.sub_spikeMutations_HR2,minVariantSize=minimumSize,maxVariantSize=maximumSize)
    adjusted_spikeMutations_TM = fc.filterMutationsBySize(variantList = fc.sub_spikeMutations_TM,minVariantSize=minimumSize,maxVariantSize=maximumSize)
    adjusted_spikeMutations_CT = fc.filterMutationsBySize(variantList = fc.sub_spikeMutations_CT,minVariantSize=minimumSize,maxVariantSize=maximumSize)
    adjusted_spikeMutations_Other = fc.filterMutationsBySize(variantList = fc.sub_spikeMutations_Other,minVariantSize=minimumSize,maxVariantSize=maximumSize)
    spikeGeneSubunitCounts = {"CT" : len(adjusted_spikeMutations_CT),
                              "Cleavage" : len(adjusted_spikeMutations_Cleavage),
                              "FP" : len(adjusted_spikeMutations_FP),
                              "HR1" : len(adjusted_spikeMutations_HR1),
                              "HR2" : len(adjusted_spikeMutations_HR2),
                              "IFP" : len(adjusted_spikeMutations_IFP),
                              "NTD" : len(adjusted_spikeMutations_NTD),
                              "RBD" : len(adjusted_spikeMutations_RBD),
                              "TM" : len(adjusted_spikeMutations_TM)}
    #endregion === Processing ===

    plt.figure(figsize=(10,6))
    plt.bar(list(spikeGeneSubunitCounts.keys()),list(spikeGeneSubunitCounts.values()),color="skyblue")

    plt.xlabel('Subunit of Spike Gene')
    plt.ylabel('Number of Mutations')
    plt.title('Distribution of SNP Mutations Across Spike Protein Subunits')
    plt.xticks(rotation=45)  # Rotate labels for better readability

    plt.show()

# Display a bar graph of the frequency of all mutations by each VOC.
def totalMutationsFrequencyByVariant():
    #region === Processing ===
    spikeGeneVariantCounts = {"Alpha" : len(fc.all_mutations_alpha),
                              "Beta" : len(fc.all_mutations_beta),
                              "Delta" : len(fc.all_mutations_delta),
                              "Gamma" : len(fc.all_mutations_gamma),
                              "OmicronBA1" : len(fc.all_mutations_omicronBA1)}

    #endregion === Processing ===

    plt.figure(figsize=(8,10))
    bars = plt.bar(list(spikeGeneVariantCounts.keys()),list(spikeGeneVariantCounts.values()),color="skyblue")

    plt.xlabel('Variant of Interest')
    plt.ylabel('Number of Total Mutations')
    plt.title('Frequency of Mutations by Variant of Interest')
    plt.xticks(rotation=45)  # Rotate labels for better readability

    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2, yval, int(yval), ha='center', va='bottom')

    plt.show()

# Display a bar graph of the frequency of spike protein mutations by each VOC.
def spikeMutationsFrequencyByVariant():
    #region === Processing ===
    spikeGeneVariantCounts = {"Alpha" : len(fc.all_spikeMutations_alpha),
                              "Beta" : len(fc.all_spikeMutations_beta),
                              "Delta" : len(fc.all_spikeMutations_delta),
                              "Gamma" : len(fc.all_spikeMutations_gamma),
                              "OmicronBA1" : len(fc.all_spikeMutations_omicronBA1)}

    #endregion === Processing ===

    plt.figure(figsize=(10,6))
    bars = plt.bar(list(spikeGeneVariantCounts.keys()),list(spikeGeneVariantCounts.values()),color="skyblue")

    plt.xlabel('Variant of Interest')
    plt.ylabel('Number of Spike Gene Mutations')
    plt.title('Frequency of Spike Gene Mutations by Variant of Interest')
    plt.xticks(rotation=45)  # Rotate labels for better readability

    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2, yval, int(yval), ha='center', va='bottom')

    plt.show()

# Display a line plot showing the transmissibility (R) rate of each VOC.
def transmissibilityByVariant():

    # region === Processing ===
    spikeGeneVariantTransmissibilities = {"Alpha" : fc.transmissibility["alpha"]["AvgR"],
                              "Beta" : fc.transmissibility["beta"]["AvgR"],
                              "Delta" : fc.transmissibility["delta"]["AvgR"],
                              "Gamma" : fc.transmissibility["gamma"]["AvgR"],
                              "OmicronBA1" : fc.transmissibility["omicron"]["AvgR"]}
    # endregion === Processing ===

    plt.figure(figsize=(10, 6))
    linePlot = plt.plot(list(spikeGeneVariantTransmissibilities.keys()), list(spikeGeneVariantTransmissibilities.values()), color="skyblue")

    plt.xlabel('Variant of Interest')
    plt.ylabel('Transmissibility (R) Factor')
    plt.title('Transmissibility Rates of Studied Variants of Interest')
    plt.xticks(rotation=45)  # Rotate labels for better readability

    # Adding the text on each point
    x_data = list(spikeGeneVariantTransmissibilities.keys())
    y_data = list(spikeGeneVariantTransmissibilities.values())
    for i, txt in enumerate(y_data):
        plt.text(x_data[i], y_data[i], f'{txt:.2f}', ha='center', va='bottom')

    plt.show()

# Display a line plot showing the transmissibility (R) rate of each VOC.
def indelTypesByVariant(variantName):

    # region === Processing ===
    numMutationsByRegion = {
        "NTD": len(fc.filterMutationsByComparedGenome(mutationList=fc.sub_spikeMutations_NTD,compareGenomeName=variantName)),
        "RBD": len(fc.filterMutationsByComparedGenome(mutationList=fc.sub_spikeMutations_RBD,compareGenomeName=variantName)),
        "Cleavage": len(fc.filterMutationsByComparedGenome(mutationList=fc.sub_spikeMutations_Cleavage,compareGenomeName=variantName)),
        "FP": len(fc.filterMutationsByComparedGenome(mutationList=fc.sub_spikeMutations_FP,compareGenomeName=variantName)),
        "IFP": len(fc.filterMutationsByComparedGenome(mutationList=fc.sub_spikeMutations_IFP,compareGenomeName=variantName)),
        "HR1": len(fc.filterMutationsByComparedGenome(mutationList=fc.sub_spikeMutations_HR1,compareGenomeName=variantName)),
        "HR2": len(fc.filterMutationsByComparedGenome(mutationList=fc.sub_spikeMutations_HR2,compareGenomeName=variantName)),
        "TM": len(fc.filterMutationsByComparedGenome(mutationList=fc.sub_spikeMutations_TM,compareGenomeName=variantName)),
        "CT": len(fc.filterMutationsByComparedGenome(mutationList=fc.sub_spikeMutations_CT,compareGenomeName=variantName)),
        "Other": len(fc.filterMutationsByComparedGenome(mutationList=fc.sub_spikeMutations_Other,compareGenomeName=variantName)),
    }
    # endregion === Processing ===

    plt.figure(figsize=(10, 6))
    bars = plt.bar(list(numMutationsByRegion.keys()), list(numMutationsByRegion.values()), color="skyblue")

    plt.xlabel('Subregion of Genome')
    plt.ylabel('Frequency of Mutations')
    plt.title(f'Frequency of Mutations by Subregion for the {variantName.capitalize()} VOC')
    plt.xticks(rotation=45)  # Rotate labels for better readability

    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2, yval, int(yval), ha='center', va='bottom')

    plt.show()

indelTypesByVariant("delta")
indelTypesByVariant("omicronBA1")