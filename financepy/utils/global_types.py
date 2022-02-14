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


class OptionTypes(Enum):
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


class SwapTypes(Enum):
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


###############################################################################

class TouchOptionTypes(Enum):
    DOWN_AND_IN_CASH_AT_HIT = 1,         # S0>H pays $1 at hit time from above
    UP_AND_IN_CASH_AT_HIT = 2,           # S0<H pays $1 at hit time from below
    DOWN_AND_IN_CASH_AT_EXPIRY = 3,      # S0>H pays $1 at T if hit from below
    UP_AND_IN_CASH_AT_EXPIRY = 4,        # S0<H pays $1 at T if hit from below
    DOWN_AND_OUT_CASH_OR_NOTHING = 5,    # S0>H pays $1 at T if S>H for all t<T
    UP_AND_OUT_CASH_OR_NOTHING = 6,      # S0<H pays $1 at T if S<H for all t<T
    DOWN_AND_IN_ASSET_AT_HIT = 7,        # S0>H pays H at hit time from above
    UP_AND_IN_ASSET_AT_HIT = 8,          # S0>H pays H at hit time from below
    DOWN_AND_IN_ASSET_AT_EXPIRY = 9,     # S0>H pays S(T) at T if S<H for t < T
    UP_AND_IN_ASSET_AT_EXPIRY = 10,      # S0<H pays S(T) at T if S>H for t < T
    DOWN_AND_OUT_ASSET_OR_NOTHING = 11,  # S0>H pays S(T) at T if S>H for t < T
    UP_AND_OUT_ASSET_OR_NOTHING = 12     # S0<H pays S(T) at T if S<H for t < T
