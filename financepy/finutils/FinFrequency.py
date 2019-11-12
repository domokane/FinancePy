# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 10:52:29 2018

@author: Dominic O'Kane
"""

from enum import Enum


class FinFrequencyTypes(Enum):
    ANNUAL = 1
    SEMI_ANNUAL = 2
    QUARTERLY = 3
    MONTHLY = 4

###############################################################################


def FinFrequency(frequencyType):

    if frequencyType in FinFrequencyTypes:
        if frequencyType == FinFrequencyTypes.ANNUAL:
            return 1
        elif frequencyType == FinFrequencyTypes.SEMI_ANNUAL:
            return 2
        elif frequencyType == FinFrequencyTypes.QUARTERLY:
            return 4
        elif frequencyType == FinFrequencyTypes.MONTHLY:
            return 12
    elif type(frequencyType) is int:
        return frequencyType
    else:
        raise ValueError("Unknown frequency type")

###############################################################################
