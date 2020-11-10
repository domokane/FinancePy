##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from enum import Enum
from ..finutils.FinError import FinError
from ..finutils.FinDate import FinDate
from ..finutils.FinGlobalVariables import gSmall

import numpy as np

###############################################################################


class FinFrequencyTypes(Enum):
    SIMPLE = 0
    ANNUAL = 1
    SEMI_ANNUAL = 2
    TRI_ANNUAL = 3
    QUARTERLY = 4
    MONTHLY = 12
    CONTINUOUS = 99

###############################################################################


def FinFrequency(freqType):
    ''' This is a function that takes in a Frequency Type and returns an
    integer for the number of times a year a payment occurs.'''
    if freqType in FinFrequencyTypes:
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
    elif isinstance(freqType, int):
        return freqType
    else:
        print("Frequency type", freqType)
        raise FinError("Unknown frequency type")

###############################################################################
