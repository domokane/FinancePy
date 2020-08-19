##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from enum import Enum
from ..finutils.FinError import FinError

import numpy as np

###############################################################################


class FinFrequencyTypes(Enum):
    CONTINUOUS = -1
    SIMPLE = 0
    ANNUAL = 1
    SEMI_ANNUAL = 2
    QUARTERLY = 3
    MONTHLY = 4

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
        raise FinError("Unknown frequency type")

###############################################################################


def zeroToDf(r: float,
             t: (float, np.ndarray),
             frequencyType: FinFrequencyTypes):
    ''' Convert a zero with a specified compounding frequency to a discount
    factor. '''

    if isinstance(t, np.ndarray):
        t = np.maximum(t, 1e-6)
    else:
        t = max(t, 1e-6)

    f = FinFrequency(frequencyType)

    if frequencyType == FinFrequencyTypes.CONTINUOUS:
        df = np.exp(-r*t)
    elif frequencyType == FinFrequencyTypes.SIMPLE:
        df = 1.0 / (1.0 + r * t)
    else:
        df = 1.0 / np.power(1.0 + r/f, f * t)

    return df

###############################################################################


def dfToZero(df: float,
             t: float,
             frequencyType: FinFrequencyTypes):
    ''' Convert a discount factor to a zero rate with a specific compounding
    frequency. '''

    if isinstance(t, np.ndarray):
        t = np.maximum(t, 1e-6)
    else:
        t = max(t, 1e-6)

    f = FinFrequency(frequencyType)

    if frequencyType == FinFrequencyTypes.CONTINUOUS:
        r = -np.log(df)/t
    elif frequencyType == FinFrequencyTypes.SIMPLE:
        r = (1.0/df - 1.0)/t
    else:
        r = (df**(-1.0/t)-1.0)*f

    return r

###############################################################################
