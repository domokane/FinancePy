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


def Frequency(freq_type):
    """ This is a function that takes in a Frequency Type and returns an
    integer for the number of times a year a payment occurs."""
    if isinstance(freq_type, FinFrequencyTypes) is False:
        print("FinFrequency:", freq_type)
        raise FinError("Unknown frequency type")

    if freq_type == FinFrequencyTypes.CONTINUOUS:
        return -1
    elif freq_type == FinFrequencyTypes.SIMPLE:
        return 0
    elif freq_type == FinFrequencyTypes.ANNUAL:
        return 1
    elif freq_type == FinFrequencyTypes.SEMI_ANNUAL:
        return 2
    elif freq_type == FinFrequencyTypes.TRI_ANNUAL:
        return 3
    elif freq_type == FinFrequencyTypes.QUARTERLY:
        return 4
    elif freq_type == FinFrequencyTypes.MONTHLY:
        return 12

###############################################################################
