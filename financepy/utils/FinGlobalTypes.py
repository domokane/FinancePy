##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# This is an exhaustive list of all option types

from enum import Enum

###############################################################################

class FinLongShort(Enum):
    LONG = 1
    SHORT = 2

###############################################################################

class FinOptionTypes(Enum):
    EUROPEAN_CALL = 1
    EUROPEAN_PUT = 2
    AMERICAN_CALL = 3
    AMERICAN_PUT = 4
    DIGITAL_CALL = 5
    DIGITAL_PUT = 6
    ASIAN_CALL = 7
    ASIAN_PUT = 8
    COMPOUND_CALL = 9
    COMPOUND_PUT = 10

###############################################################################

class FinCapFloorTypes(Enum):
    CAP = 1
    FLOOR = 2

###############################################################################

class FinSwapTypes(Enum):
    PAY = 1
    RECEIVE = 2

###############################################################################

class FinExerciseTypes(Enum):
    EUROPEAN = 1
    BERMUDAN = 2
    AMERICAN = 3

###############################################################################

class FinSolverTypes(Enum):
    CONJUGATE_GRADIENT = 0
    NELDER_MEAD = 1
    NELDER_MEAD_NUMBA = 2
