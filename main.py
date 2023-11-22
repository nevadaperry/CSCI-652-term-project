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
if(True):
    minimumSize = 0
    maximumSize = 1

    #region === Processing ===
    adjusted_spikeVariants_NTD = fc.filterVariantsBySize(variantList = fc.all_spikeVariants_NTD,minVariantSize=minimumSize,maxVariantSize=maximumSize)
    adjusted_spikeVariants_RBD = fc.filterVariantsBySize(variantList = fc.all_spikeVariants_RBD,minVariantSize=minimumSize,maxVariantSize=maximumSize)
    adjusted_spikeVariants_Cleavage = fc.filterVariantsBySize(variantList = fc.all_spikeVariants_Cleavage,minVariantSize=minimumSize,maxVariantSize=maximumSize)
    adjusted_spikeVariants_FP = fc.filterVariantsBySize(variantList = fc.all_spikeVariants_FP,minVariantSize=minimumSize,maxVariantSize=maximumSize)
    adjusted_spikeVariants_IFP = fc.filterVariantsBySize(variantList = fc.all_spikeVariants_IFP,minVariantSize=minimumSize,maxVariantSize=maximumSize)
    adjusted_spikeVariants_HR1 = fc.filterVariantsBySize(variantList = fc.all_spikeVariants_HR1,minVariantSize=minimumSize,maxVariantSize=maximumSize)
    adjusted_spikeVariants_HR2 = fc.filterVariantsBySize(variantList = fc.all_spikeVariants_HR2,minVariantSize=minimumSize,maxVariantSize=maximumSize)
    adjusted_spikeVariants_TM = fc.filterVariantsBySize(variantList = fc.all_spikeVariants_TM,minVariantSize=minimumSize,maxVariantSize=maximumSize)
    adjusted_spikeVariants_CT = fc.filterVariantsBySize(variantList = fc.all_spikeVariants_CT,minVariantSize=minimumSize,maxVariantSize=maximumSize)
    adjusted_spikeVariants_Other = fc.filterVariantsBySize(variantList = fc.all_spikeVariants_Other,minVariantSize=minimumSize,maxVariantSize=maximumSize)
    spikeGeneSubunitCounts = {"CT" : len(adjusted_spikeVariants_CT),
                              "Cleavage" : len(adjusted_spikeVariants_Cleavage),
                              "FP" : len(adjusted_spikeVariants_FP),
                              "HR1" : len(adjusted_spikeVariants_HR1),
                              "HR2" : len(adjusted_spikeVariants_HR2),
                              "IFP" : len(adjusted_spikeVariants_IFP),
                              "NTD" : len(adjusted_spikeVariants_NTD),
                              "RBD" : len(adjusted_spikeVariants_RBD),
                              "TM" : len(adjusted_spikeVariants_TM)}
    #endregion === Processing ===

    plt.figure(figsize=(10,6))
    plt.bar(list(spikeGeneSubunitCounts.keys()),list(spikeGeneSubunitCounts.values()),color="skyblue")

    plt.xlabel('Subunit of Spike Gene')
    plt.ylabel('Number of Variants')
    plt.title('Distribution of SNP Variants Across Spike Protein Subunits')
    plt.xticks(rotation=45)  # Rotate labels for better readability

    plt.show()