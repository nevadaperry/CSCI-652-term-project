"""
Author: Alex Somheil
Student ID: Z1888439
Date: 11/20/2023

This file contains functions used to help visualize some of the higher level functional
consequences resulted from comparing SARS-COV2 to its variants.
"""
import Setup as s
maxVariantSize = 10**100

# Given a list of variant dicts, this method returns a similar list ordered
# by size. Optional parameter to omit variants smaller than a certain size.
def filterVariantsBySize(variantList : list,minVariantSize : int = 0,maxVariantSize : int = maxVariantSize):
    sortedVariants = [item for item in sorted(variantList, key=lambda x: x['Length'], reverse=True) if minVariantSize <= item['Length'] <= maxVariantSize]
    return sortedVariants

# Given a list of variant dicts, this method returns a similar list filtered only for
# a specific unit (and option subunit). Optional booleans to set unit/subunit as blacklist instead of
# whitelist.
def filterVariantsByUnit(variantList : list, unit : str | list = None, subunit : str | list = None, unitBlackList : bool = False, subunitBlackList : bool = False):
    # Convert unit and subunit to lists if they are not already.
    if(type(unit) is not list):
        unit = [unit] if unit is not None else []
    if(type(subunit) is not list):
        subunit = [subunit] if subunit is not None else []

    filteredVariants = variantList
    # Filtering logic for units.
    if(unit):
        newFilteredList = []
        for variant in filteredVariants:
            for thisUnitFilter in unit:
                if(thisUnitFilter == variant["Unit"]):
                    if(unitBlackList):
                        continue
                    else:
                        newFilteredList.append(variant)
        filteredVariants = newFilteredList

    # Filtering logic for subunits.
    if(subunit):
        newFilteredList = []
        for variant in filteredVariants:
            for thisSubunitFilter in subunit:
                if(thisSubunitFilter == variant["SubUnit"]):
                    if(subunitBlackList):
                        break
                    else:
                        newFilteredList.append(variant)
                        break
        filteredVariants = newFilteredList

    return sorted(filteredVariants, key=lambda x: x['Unit'], reverse=True)

# Given a list of variant dicts, this method returns a similar list filtered for variants that occur between a
# specific minIndex and maxIndex
def filterVariantsByLocation(variantList : list, minIndex : int = 0, maxIndex : int = maxVariantSize):
    sortedVariants = [item for item in sorted(variantList, key=lambda x: x[f"{x['Genome1']}-Location"][0], reverse=True) if minIndex <= item[f"{item['Genome1']}-Location"][0] <= maxIndex]
    return sortedVariants

# Given a list of variants dicts, this method returns a similar list filtered for variants that occur only for a
# specific genome2 (compared genome)
def filterVariantsByComparedGenome(variantList : list, compareGenomeName : str):
    sortedVariants = [item for item in variantList if item['Genome2'].lower() == compareGenomeName.lower()]
    return sortedVariants


# Concatenating all found variants together into one list for ease of use.
variants = []
for _pwDictName,_pwDict in s.data["pw"].items():
    variants += _pwDict["Variants"]

# Filtering our variant list to only view variants that occur on the Spike gene, and within a specific subunit.
spikeVariants = filterVariantsByUnit(variantList = variants, unit = "S",)

# Creating lists of variants that occur on each specific subUnit of the Spike gene.
sub_spikeVariants_NTD = filterVariantsByUnit(variantList=spikeVariants,subunit="NTD")
sub_spikeVariants_RBD = filterVariantsByUnit(variantList=spikeVariants,subunit="RBD")
sub_spikeVariants_Cleavage = filterVariantsByUnit(variantList=spikeVariants,subunit="Cleavage")
sub_spikeVariants_FP = filterVariantsByUnit(variantList=spikeVariants,subunit="FP")
sub_spikeVariants_IFP = filterVariantsByUnit(variantList=spikeVariants,subunit="IFP")
sub_spikeVariants_HR1 = filterVariantsByUnit(variantList=spikeVariants,subunit="HR1")
sub_spikeVariants_HR2 = filterVariantsByUnit(variantList=spikeVariants,subunit="HR2")
sub_spikeVariants_TM = filterVariantsByUnit(variantList=spikeVariants,subunit="TM")
sub_spikeVariants_CT = filterVariantsByUnit(variantList=spikeVariants,subunit="CT")
sub_spikeVariants_Other = filterVariantsByUnit(variantList=spikeVariants,subunit="Other")

all_spikeVariants_alpha = filterVariantsByComparedGenome(variantList=spikeVariants,compareGenomeName="alpha")
all_spikeVariants_beta = filterVariantsByComparedGenome(variantList=spikeVariants,compareGenomeName="beta")
all_spikeVariants_delta = filterVariantsByComparedGenome(variantList=spikeVariants,compareGenomeName="delta")
all_spikeVariants_epsilon = filterVariantsByComparedGenome(variantList=spikeVariants,compareGenomeName="epsilon")
all_spikeVariants_eta = filterVariantsByComparedGenome(variantList=spikeVariants,compareGenomeName="eta")
all_spikeVariants_gamma = filterVariantsByComparedGenome(variantList=spikeVariants,compareGenomeName="gamma")
all_spikeVariants_iota = filterVariantsByComparedGenome(variantList=spikeVariants,compareGenomeName="iota")
all_spikeVariants_kappa = filterVariantsByComparedGenome(variantList=spikeVariants,compareGenomeName="kappa")
all_spikeVariants_lambda = filterVariantsByComparedGenome(variantList=spikeVariants,compareGenomeName="lambda")
all_spikeVariants_mu = filterVariantsByComparedGenome(variantList=spikeVariants,compareGenomeName="mu")
all_spikeVariants_omicron = filterVariantsByComparedGenome(variantList=spikeVariants,compareGenomeName="omicron")
all_spikeVariants_omicronBA1 = filterVariantsByComparedGenome(variantList=spikeVariants,compareGenomeName="omicronBA1")
all_spikeVariants_omicronBA2 = filterVariantsByComparedGenome(variantList=spikeVariants,compareGenomeName="omicronBA2")
all_spikeVariants_omicronBA3 = filterVariantsByComparedGenome(variantList=spikeVariants,compareGenomeName="omicronBA3")
all_spikeVariants_omicronBA4 = filterVariantsByComparedGenome(variantList=spikeVariants,compareGenomeName="omicronBA4")
all_spikeVariants_omicronBA5 = filterVariantsByComparedGenome(variantList=spikeVariants,compareGenomeName="omicronBA5")
all_spikeVariants_omicronBQ1 = filterVariantsByComparedGenome(variantList=spikeVariants,compareGenomeName="omicronBQ1")
all_spikeVariants_omicronXAK = filterVariantsByComparedGenome(variantList=spikeVariants,compareGenomeName="omicronXAK")
all_spikeVariants_omicronXBB1 = filterVariantsByComparedGenome(variantList=spikeVariants,compareGenomeName="omicronXBB1")
all_spikeVariants_theta = filterVariantsByComparedGenome(variantList=spikeVariants,compareGenomeName="theta")
all_spikeVariants_zeta = filterVariantsByComparedGenome(variantList=spikeVariants,compareGenomeName="zeta")

