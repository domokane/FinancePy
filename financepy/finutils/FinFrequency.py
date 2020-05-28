##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################



from enum import Enum


class FinFrequencyTypes(Enum):
    ANNUAL = 1
    SEMI_ANNUAL = 2
    QUARTERLY = 3
    MONTHLY = 4

###############################################################################


def FinFrequency(frequencyType):
    ''' This is a function that takes in a Frequency Type and returns an
    integer for the number of times a year a payment occurs.'''
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
