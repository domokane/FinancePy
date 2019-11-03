# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 10:52:29 2018

@author: Dominic O'Kane
"""

from . import FinError

from enum import Enum


class FinFrequencyTypes(Enum):
    ANNUAL = 1
    SEMI_ANNUAL = 2
    QUARTERLY = 3
    MONTHLY = 4

###############################################################################


def FinFrequency(frequencyType):

    if frequencyType == FinFrequencyTypes.ANNUAL:
        frequency = 1
    elif frequencyType == FinFrequencyTypes.SEMI_ANNUAL:
        frequency = 2
    elif frequencyType == FinFrequencyTypes.QUARTERLY:
        frequency = 4
    elif frequencyType == FinFrequencyTypes.MONTHLY:
        frequency = 12
    else:
        raise FinError("Frequency Type not found:" + str(frequencyType))

    return frequency

###############################################################################

    def __str__(self):
        return str(self._name)

###############################################################################
