"""
Author: Alex Somheil
Student ID: Z1888439
Date: 11/20/2023

This file contains functions used to help visualize some of the higher level functional
consequences resulted from comparing SARS-COV2 to its variants.
"""
import Setup as s


# Simple method to return a list of all variants found in a pairwise alignment dict, ordered
# by size. Optional parameter to omit variants smaller than a certain size.
def sortVariantsBySize(pwDict : dict,omitVariantsSmallerThan : int = 0):
    sortedVariants = [item for item in sorted(pwDict["Variants"], key=lambda x: x['Length'], reverse=True) if item['Length'] >= omitVariantsSmallerThan]
    return sortedVariants




everyVariant = []
for _pwDictName,_pwDict in s.data["pw"].items():
    everyVariant += sortVariantsBySize(pwDict = _pwDict,omitVariantsSmallerThan=10)
everyVariant = [item for item in sorted(everyVariant, key=lambda x: x['Length'], reverse=True)]

