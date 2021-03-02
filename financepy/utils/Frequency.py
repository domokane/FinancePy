##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ..utils.FinError import FinError

###############################################################################

from enum import Enum


class FrequencyTypes(Enum):
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
    if isinstance(freq_type, FrequencyTypes) is False:
        print("FinFrequency:", freq_type)
        raise FinError("Unknown frequency type")

    if freq_type == FrequencyTypes.CONTINUOUS:
        return -1
    elif freq_type == FrequencyTypes.SIMPLE:
        return 0
    elif freq_type == FrequencyTypes.ANNUAL:
        return 1
    elif freq_type == FrequencyTypes.SEMI_ANNUAL:
        return 2
    elif freq_type == FrequencyTypes.TRI_ANNUAL:
        return 3
    elif freq_type == FrequencyTypes.QUARTERLY:
        return 4
    elif freq_type == FrequencyTypes.MONTHLY:
        return 12

###############################################################################
