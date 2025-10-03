##############################################################################
# Copyright (C) 2018 2019 2020 Dominic O'Kane
##############################################################################

# This is an exhaustive list of all option types

from enum import Enum

########################################################################################


class FinLongShort(Enum):
    """Enumeration of long/short positions."""

    LONG = 1
    SHORT = 2


########################################################################################


class DoubleBarrierTypes(Enum):
    """Enumeration of long/short positions."""

    KNOCK_IN = 1
    KNOCK_OUT = 2


########################################################################################


class OptionTypes(Enum):
    """Enumeration of option types."""

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


########################################################################################


class EquityBarrierTypes(Enum):
    """Enumeration of equity barrier types."""

    DOWN_AND_OUT_CALL = 1
    DOWN_AND_IN_CALL = 2
    UP_AND_OUT_CALL = 3
    UP_AND_IN_CALL = 4
    UP_AND_OUT_PUT = 5
    UP_AND_IN_PUT = 6
    DOWN_AND_OUT_PUT = 7
    DOWN_AND_IN_PUT = 8


########################################################################################


class FinCapFloorTypes(Enum):
    """Enumeration of cap/floor types."""

    CAP = 1
    FLOOR = 2


########################################################################################


class SwapTypes(Enum):
    """Enumeration of swap types."""

    PAY = 1
    RECEIVE = 2


########################################################################################


class ReturnTypes(Enum):
    """Enumeration of return types."""

    TOTAL_RETURN = 1
    PRICE_RETURN = 2


########################################################################################


class FinExerciseTypes(Enum):
    """Enumeration of exercise types."""

    EUROPEAN = 1
    BERMUDAN = 2
    AMERICAN = 3


########################################################################################


class FinSolverTypes(Enum):
    """Enumeration of solver types."""

    CONJUGATE_GRADIENT = 0
    NELDER_MEAD = 1
    NELDER_MEAD_NUMBA = 2


########################################################################################


class TouchOptionTypes(Enum):
    """Enumeration of touch option types."""

    DOWN_AND_IN_CASH_AT_HIT = 1  # s0>H pays $1 at hit time from above
    UP_AND_IN_CASH_AT_HIT = 2  # s0<H pays $1 at hit time from below
    DOWN_AND_IN_CASH_AT_EXPIRY = 3  # s0>H pays $1 at T if hit from below
    UP_AND_IN_CASH_AT_EXPIRY = 4  # s0<H pays $1 at T if hit from below
    DOWN_AND_OUT_CASH_OR_NOTHING = 5  # s0>H pays $1 at T if S>H for all t<T
    UP_AND_OUT_CASH_OR_NOTHING = 6  # s0<H pays $1 at T if S<H for all t<T
    DOWN_AND_IN_ASSET_AT_HIT = 7  # s0>H pays H at hit time from above
    UP_AND_IN_ASSET_AT_HIT = 8  # s0>H pays H at hit time from below
    DOWN_AND_IN_ASSET_AT_EXPIRY = 9  # s0>H pays S(T) at T if S<H for t < T
    UP_AND_IN_ASSET_AT_EXPIRY = 10  # s0<H pays S(T) at T if S>H for t < T
    DOWN_AND_OUT_ASSET_OR_NOTHING = 11  # s0>H pays S(T) at T if S>H for t < T
    UP_AND_OUT_ASSET_OR_NOTHING = 12  # s0<H pays S(T) at T if S<H for t < T
