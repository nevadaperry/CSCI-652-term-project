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
def filterVariantsByUnit(variantList : list, unit : str | list, subunit : str | list = None, unitBlacklist : bool = False, subunitBlackList : bool = False):
    # Convert unit and subunit to lists if they are not already.
    if(type(unit) is list):
        unit = [unit] if unit is not None else []
    if(type(subunit) is list):
        subunit = [subunit] if subunit is not None else []

    # Filtering logic for units.
    if unitBlacklist:
        filteredVariants = [item for item in variantList if item['Gene'] not in unit]
    else:
        filteredVariants = [item for item in variantList if item['Gene'] in unit]

    # Filtering logic for subunits.
    if subunit:
        if subunitBlackList:
            filteredVariants = [item for item in filteredVariants if item['SubUnit'] not in subunit]
        else:
            filteredVariants = [item for item in filteredVariants if item['SubUnit'] in subunit]

    return sorted(filteredVariants, key=lambda x: x['Gene'], reverse=True)

# Concatenating all found variants together into one list for ease of use.
everyVariant = []
for _pwDictName,_pwDict in s.data["pw"].items():
    everyVariant += _pwDict["Variants"]

# Filtering out variants of size smaller than X
everyVariant = sortVariantsBySize(variantList = everyVariant,omitVariantsSmallerThan=0)

# Filtering our variant list to only view variants that occur on the Spike gene, and within a specific subunit.
everyVariant = filterVariantsByUnit(variantList = everyVariant, unit = "S",subunit="Other",subunitBlackList=True)


