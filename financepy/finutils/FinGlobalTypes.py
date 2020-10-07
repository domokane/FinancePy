##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# This is an exhaustive list of all option types

from enum import Enum

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
    PAYER = 1
    RECEIVER = 2

###############################################################################

class FinExerciseTypes(Enum):
    EUROPEAN = 1
    BERMUDAN = 2
    AMERICAN = 3

###############################################################################
