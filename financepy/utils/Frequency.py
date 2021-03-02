##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ..utils.FinError import FinError

###############################################################################

from enum import Enum


class FinFrequencyTypes(Enum):
    SIMPLE = 0
    ANNUAL = 1
    SEMI_ANNUAL = 2
    TRI_ANNUAL = 3
    QUARTERLY = 4
    MONTHLY = 12
    CONTINUOUS = 99

###############################################################################


def Frequency(freqType):
    ''' This is a function that takes in a Frequency Type and returns an
    integer for the number of times a year a payment occurs.'''
    if isinstance(freqType, FinFrequencyTypes) is False:
        print("FinFrequency:", freqType)
        raise FinError("Unknown frequency type")

    if freqType == FinFrequencyTypes.CONTINUOUS:
        return -1
    elif freqType == FinFrequencyTypes.SIMPLE:
        return 0
    elif freqType == FinFrequencyTypes.ANNUAL:
        return 1
    elif freqType == FinFrequencyTypes.SEMI_ANNUAL:
        return 2
    elif freqType == FinFrequencyTypes.TRI_ANNUAL:
        return 3
    elif freqType == FinFrequencyTypes.QUARTERLY:
        return 4
    elif freqType == FinFrequencyTypes.MONTHLY:
        return 12

###############################################################################
