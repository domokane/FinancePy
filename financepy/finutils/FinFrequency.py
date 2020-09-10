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
    CONTINUOUS = 1
    SIMPLE = 2
    ANNUAL = 3
    SEMI_ANNUAL = 4
    QUARTERLY = 5
    MONTHLY = 6

###############################################################################


def FinFrequency(frequencyType):
    ''' This is a function that takes in a Frequency Type and returns an
    integer for the number of times a year a payment occurs.'''
    if frequencyType in FinFrequencyTypes:
        if frequencyType == FinFrequencyTypes.CONTINUOUS:
            return -1
        elif frequencyType == FinFrequencyTypes.SIMPLE:
            return 0
        elif frequencyType == FinFrequencyTypes.ANNUAL:
            return 1
        elif frequencyType == FinFrequencyTypes.SEMI_ANNUAL:
            return 2
        elif frequencyType == FinFrequencyTypes.QUARTERLY:
            return 4
        elif frequencyType == FinFrequencyTypes.MONTHLY:
            return 12
    elif isinstance(frequencyType, int):
        return frequencyType
    else:
        print("Frequency type", frequencyType)
        raise FinError("Unknown frequency type")

###############################################################################
