##############################################################################
# Copyright (C) 2018 2019 2020 Dominic O'Kane
##############################################################################

# This is a centralised exhaustive list of all FinancePy types

from enum import Enum

########################################################################################


class LMMModelTypes(Enum):

    ONE_FACTOR = 1
    HW_M_FACTOR = 2
    FULL_N_FACTOR = 3


class BlackTypes(Enum):
    ANALYTICAL = 1
    CRR_TREE = 2


class HestonNumericalSchemeTypes(Enum):
    EULER = 1
    EULERLOG = 2
    QUADEXP = 3

class CIRNumericalSchemeTypes(Enum):
    EULER = 1
    LOGNORMAL = 2
    MILSTEIN = 3
    KAHLJACKEL = 4
    EXACT = 5


class GBMNumericalSchemeTypes(Enum):
    NORMAL = 1
    ANTITHETIC = 2

class VolFuncTypes(Enum):

    CLARK = 0
    SABR = 1
    SABR_BETA_ONE = 2
    SABR_BETA_HALF = 3
    BBG = 4
    CLARK5 = 5
    SVI = 6
    SSVI = 7


class YTMCalcType(Enum):
    ZERO = 0
    UK_DMO = 1
    US_STREET = 2
    US_TREASURY = 3
    CFETS = 4  # China Foreign Exchange Trade System
    CALCULUS = 5  # Using Calculus for duration


class InterpTypes(Enum):
    FLAT_FWD_RATES = 1
    LINEAR_DISCOUNT = 2
    LINEAR_ZERO_RATES = 3
    LINEAR_ONFWD_RATES = 4
    FINCUBIC_ZERO_RATES = 5
    NATCUBIC_LOG_DISCOUNT = 6
    NATCUBIC_ZERO_RATES = 7
    PCHIP_ZERO_RATES = 8
    PCHIP_LOG_DISCOUNT = 9
    TENSION_ZERO_RATES = 10


class BlackScholesTypes(Enum):

    DEFAULT = 0
    ANALYTICAL = 1
    CRR_TREE = 2
    BARONE_ADESI = 3
    LSMC = 4
    BJERKSUND_STENSLAND = 5
    FINITE_DIFFERENCE = 6
    PSOR = 7


class ProcessTypes(Enum):

    GBM_PROCESS = 1
    CIR_PROCESS = 2
    HESTON_PROCESS = 3
    VASICEK_PROCESS = 4
    CEV_PROCESS = 5
    JUMP_DIFFUSION_PROCESS = 6
    
class VasicekNumericalSchemeTypes(Enum):
    NORMAL = 1
    ANTITHETIC = 2

class FXBarrierTypes(Enum):
    DOWN_AND_OUT_CALL = 1
    DOWN_AND_IN_CALL = 2
    UP_AND_OUT_CALL = 3
    UP_AND_IN_CALL = 4
    UP_AND_OUT_PUT = 5
    UP_AND_IN_PUT = 6
    DOWN_AND_OUT_PUT = 7
    DOWN_AND_IN_PUT = 8


class OISCompoundingTypes(Enum):
    COMPOUNDED = 1
    OVERNIGHT_COMPOUNDED_ANNUAL_RATE = 2
    AVERAGED = 3
    AVERAGED_DAILY = 4


class HWEuropeanCalcTypes(Enum):

    JAMSHIDIAN = 1
    EXPIRY_ONLY = 2
    EXPIRY_TREE = 3


class FXATMMethodTypes(Enum):
    SPOT = 1  # Spot FX Rate
    FWD = 2  # At the money forward
    FWD_DELTA_NEUTRAL = 3  # K = F*exp(0.5*sigma*sigma*T)
    FWD_DELTA_NEUTRAL_PREM_ADJ = 4  # K = F*exp(-0.5*sigma*sigma*T)


class FXDeltaMethodTypes(Enum):
    SPOT_DELTA = 1
    FORWARD_DELTA = 2
    SPOT_DELTA_PREM_ADJ = 3
    FORWARD_DELTA_PREM_ADJ = 4


class LongShortTypes(Enum):
    """Enumeration of long/short positions."""

    LONG = 1
    SHORT = 2


########################################################################################


class DoubleBarrierTypes(Enum):
    """Enumeration of long/short positions."""

    KNOCK_IN = 1
    KNOCK_OUT = 2


########################################################################################

class DigitalOptionTypes(Enum):
    CASH_OR_NOTHING = 1
    ASSET_OR_NOTHING = 2


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


class BarrierTypes(Enum):
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


class CapFloorTypes(Enum):
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


class ExerciseTypes(Enum):
    """Enumeration of exercise types."""

    EUROPEAN = 1
    BERMUDAN = 2
    AMERICAN = 3


########################################################################################


class SolverTypes(Enum):
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
