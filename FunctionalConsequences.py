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

def filterVariantsByUnit(variantList: list, unit, subunit=None, unitBlacklist=False, subunitBlackList=False):
    # Convert unit and subunit to lists if they are not already.
    if not isinstance(unit, list):
        unit = [unit] if unit is not None else []
    if not isinstance(subunit, list):
        subunit = [subunit] if subunit is not None else []

    # Filtering logic for units.
    if unitBlacklist:
        filteredVariants = [item for item in variantList if item['Gene'] not in unit]
    else:
        filteredVariants = [item for item in variantList if item['Gene'] in unit]

    # Filtering logic for subunits.
    if subunit:
        if subunitBlackList:
            sortedVariants = [item for item in filteredVariants if item['SubUnit'] not in subunit]
        else:
            sortedVariants = [item for item in filteredVariants if item['SubUnit'] in subunit]
    else:
        sortedVariants = filteredVariants

    return sorted(sortedVariants, key=lambda x: x['Gene'], reverse=True)


everyVariant = []
for _pwDictName,_pwDict in s.data["pw"].items():
    everyVariant += _pwDict["Variants"]

everyVariant = sortVariantsBySize(variantList = everyVariant,omitVariantsSmallerThan=0)
everyVariant = filterVariantsByUnit(variantList = everyVariant, unit = "S",subunit="Other",subunitBlackList=True)


