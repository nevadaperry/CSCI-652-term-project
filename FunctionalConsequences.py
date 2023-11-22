"""
Author: Alex Somheil
Student ID: Z1888439
Date: 11/20/2023

This file contains functions used to help visualize some of the higher level functional
consequences resulted from comparing SARS-COV2 to its variants.
"""
import Setup as s


# Given a list of variant dicts, this method returns a similar list ordered
# by size. Optional parameter to omit variants smaller than a certain size.
def sortVariantsBySize(variantList : list,omitVariantsSmallerThan : int = 0):
    sortedVariants = [item for item in sorted(variantList, key=lambda x: x['Length'], reverse=True) if item['Length'] >= omitVariantsSmallerThan]
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

# Concatenating all found variants together into one list for ease of use.
variants = []
for _pwDictName,_pwDict in s.data["pw"].items():
    variants += _pwDict["Variants"]

# Sorting variant list by size
variants = sortVariantsBySize(variantList = variants,omitVariantsSmallerThan=0)

# Filtering our variant list to only view variants that occur on the Spike gene, and within a specific subunit.
spikeVariants = filterVariantsByUnit(variantList = variants, unit = "S",)


spikeVariants_NTD = filterVariantsByUnit(variantList=spikeVariants,subunit="NTD")
spikeVariants_RBD = filterVariantsByUnit(variantList=spikeVariants,subunit="RBD")
spikeVariants_Cleavage = filterVariantsByUnit(variantList=spikeVariants,subunit="Cleavage")
spikeVariants_FP = filterVariantsByUnit(variantList=spikeVariants,subunit="FP")
spikeVariants_IFP = filterVariantsByUnit(variantList=spikeVariants,subunit="IFP")
spikeVariants_HR1 = filterVariantsByUnit(variantList=spikeVariants,subunit="HR1")
spikeVariants_HR2 = filterVariantsByUnit(variantList=spikeVariants,subunit="HR2")
spikeVariants_TM = filterVariantsByUnit(variantList=spikeVariants,subunit="TM")
spikeVariants_CT = filterVariantsByUnit(variantList=spikeVariants,subunit="CT")
spikeVariants_Other = filterVariantsByUnit(variantList=spikeVariants,subunit="Other")