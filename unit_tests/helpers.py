###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.global_types import SwapTypes
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.credit.cds_curve import CDSCurve
from financepy.products.rates.ibor_fra import IborFRA
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.products.rates.ibor_swap import IborSwap
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.credit.cds import CDS

from os.path import dirname, join


def build_Ibor_Curve(tradeDate):

    value_dt = tradeDate.add_days(1)
    dc_type = DayCountTypes.ACT_360

    depos = []
    fras = []
    swaps = []

    dc_type = DayCountTypes.THIRTY_E_360_ISDA
    fixed_freq = FrequencyTypes.SEMI_ANNUAL
    settle_dt = value_dt

    maturity_dt = settle_dt.add_months(12)
    swap1 = IborSwap(
        settle_dt, maturity_dt, SwapTypes.PAY, 0.0502, fixed_freq, dc_type
    )
    swaps.append(swap1)

    maturity_dt = settle_dt.add_months(24)
    swap2 = IborSwap(
        settle_dt, maturity_dt, SwapTypes.PAY, 0.0502, fixed_freq, dc_type
    )
    swaps.append(swap2)

    maturity_dt = settle_dt.add_months(36)
    swap3 = IborSwap(
        settle_dt, maturity_dt, SwapTypes.PAY, 0.0501, fixed_freq, dc_type
    )
    swaps.append(swap3)

    maturity_dt = settle_dt.add_months(48)
    swap4 = IborSwap(
        settle_dt, maturity_dt, SwapTypes.PAY, 0.0502, fixed_freq, dc_type
    )
    swaps.append(swap4)

    maturity_dt = settle_dt.add_months(60)
    swap5 = IborSwap(
        settle_dt, maturity_dt, SwapTypes.PAY, 0.0501, fixed_freq, dc_type
    )
    swaps.append(swap5)

    libor_curve = IborSingleCurve(value_dt, depos, fras, swaps)

    return libor_curve


def buildIssuerCurve(tradeDate, libor_curve):

    value_dt = tradeDate.add_days(1)

    cdsMarketContracts = []

    cds_cpn = 0.0048375
    maturity_dt = Date(29, 6, 2010)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    recovery_rate = 0.40

    issuer_curve = CDSCurve(
        value_dt, cdsMarketContracts, libor_curve, recovery_rate
    )

    return issuer_curve


def buildFlatIssuerCurve(tradeDate, libor_curve, spread, recovery_rate):

    value_dt = tradeDate.add_days(1)

    cdsMarketContracts = []

    maturity_dt = Date(29, 6, 2010)
    cds = CDS(value_dt, maturity_dt, spread)
    cdsMarketContracts.append(cds)

    issuer_curve = CDSCurve(
        value_dt, cdsMarketContracts, libor_curve, recovery_rate
    )

    return issuer_curve


def buildFullIssuerCurve(value_dt):

    dc_type = DayCountTypes.ACT_360
    depos = []
    irBump = 0.0

    m = 1.0  # 0.00000000000

    spot_days = 0
    settle_dt = value_dt.add_days(spot_days)

    maturity_dt = settle_dt.add_months(1)
    depo1 = IborDeposit(settle_dt, maturity_dt, m * 0.0016, dc_type)

    maturity_dt = settle_dt.add_months(2)
    depo2 = IborDeposit(settle_dt, maturity_dt, m * 0.0020, dc_type)

    maturity_dt = settle_dt.add_months(3)
    depo3 = IborDeposit(settle_dt, maturity_dt, m * 0.0024, dc_type)

    maturity_dt = settle_dt.add_months(6)
    depo4 = IborDeposit(settle_dt, maturity_dt, m * 0.0033, dc_type)

    maturity_dt = settle_dt.add_months(12)
    depo5 = IborDeposit(settle_dt, maturity_dt, m * 0.0056, dc_type)

    depos.append(depo1)
    depos.append(depo2)
    depos.append(depo3)
    depos.append(depo4)
    depos.append(depo5)

    fras = []

    spot_days = 2
    settle_dt = value_dt.add_days(spot_days)

    swaps = []
    dc_type = DayCountTypes.THIRTY_E_360_ISDA
    fixed_freq = FrequencyTypes.SEMI_ANNUAL

    maturity_dt = settle_dt.add_months(24)
    swap1 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0044 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap1)

    maturity_dt = settle_dt.add_months(36)
    swap2 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0078 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap2)

    maturity_dt = settle_dt.add_months(48)
    swap3 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0119 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap3)

    maturity_dt = settle_dt.add_months(60)
    swap4 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0158 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap4)

    maturity_dt = settle_dt.add_months(72)
    swap5 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0192 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap5)

    maturity_dt = settle_dt.add_months(84)
    swap6 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0219 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap6)

    maturity_dt = settle_dt.add_months(96)
    swap7 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0242 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap7)

    maturity_dt = settle_dt.add_months(108)
    swap8 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0261 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap8)

    maturity_dt = settle_dt.add_months(120)
    swap9 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0276 + irBump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap9)

    libor_curve = IborSingleCurve(value_dt, depos, fras, swaps)

    cdsMarketContracts = []
    cds_cpn = 0.005743
    maturity_dt = value_dt.next_cds_date(6)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    cds_cpn = 0.007497
    maturity_dt = value_dt.next_cds_date(12)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    cds_cpn = 0.011132
    maturity_dt = value_dt.next_cds_date(24)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    cds_cpn = 0.013932
    maturity_dt = value_dt.next_cds_date(36)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    cds_cpn = 0.015764
    maturity_dt = value_dt.next_cds_date(48)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    cds_cpn = 0.017366
    maturity_dt = value_dt.next_cds_date(60)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    cds_cpn = 0.020928
    maturity_dt = value_dt.next_cds_date(84)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    cds_cpn = 0.022835
    maturity_dt = value_dt.next_cds_date(120)
    cds = CDS(value_dt, maturity_dt, cds_cpn)
    cdsMarketContracts.append(cds)

    recovery_rate = 0.40

    issuer_curve = CDSCurve(
        value_dt, cdsMarketContracts, libor_curve, recovery_rate
    )

    return libor_curve, issuer_curve


def loadHomogeneousSpreadCurves(
    value_dt,
    libor_curve,
    cdsSpread3Y,
    cdsSpread5Y,
    cdsSpread7Y,
    cdsSpread10Y,
    num_credits,
):

    maturity3Y = value_dt.next_cds_date(36)
    maturity5Y = value_dt.next_cds_date(60)
    maturity7Y = value_dt.next_cds_date(84)
    maturity10Y = value_dt.next_cds_date(120)

    recovery_rate = 0.40

    cds3Y = CDS(value_dt, maturity3Y, cdsSpread3Y)
    cds5Y = CDS(value_dt, maturity5Y, cdsSpread5Y)
    cds7Y = CDS(value_dt, maturity7Y, cdsSpread7Y)
    cds10Y = CDS(value_dt, maturity10Y, cdsSpread10Y)

    contracts = [cds3Y, cds5Y, cds7Y, cds10Y]

    issuer_curve = CDSCurve(value_dt, contracts, libor_curve, recovery_rate)

    issuer_curves = []
    for _ in range(0, num_credits):
        issuer_curves.append(issuer_curve)

    return issuer_curves


def loadHeterogeneousSpreadCurves(value_dt, libor_curve):

    maturity3Y = value_dt.next_cds_date(36)
    maturity5Y = value_dt.next_cds_date(60)
    maturity7Y = value_dt.next_cds_date(84)
    maturity10Y = value_dt.next_cds_date(120)

    path = dirname(__file__)
    filename = "CDX_NA_IG_S7_SPREADS.csv"
    full_filename_path = join(path, "data", filename)
    f = open(full_filename_path, "r")

    data = f.readlines()
    issuer_curves = []

    for row in data[1:]:

        splitRow = row.split(",")
        spd3Y = float(splitRow[1]) / 10000.0
        spd5Y = float(splitRow[2]) / 10000.0
        spd7Y = float(splitRow[3]) / 10000.0
        spd10Y = float(splitRow[4]) / 10000.0
        recovery_rate = float(splitRow[5])

        cds3Y = CDS(value_dt, maturity3Y, spd3Y)
        cds5Y = CDS(value_dt, maturity5Y, spd5Y)
        cds7Y = CDS(value_dt, maturity7Y, spd7Y)
        cds10Y = CDS(value_dt, maturity10Y, spd10Y)
        cds_contracts = [cds3Y, cds5Y, cds7Y, cds10Y]

        issuer_curve = CDSCurve(
            value_dt, cds_contracts, libor_curve, recovery_rate
        )

        issuer_curves.append(issuer_curve)

    return issuer_curves


def loadHomogeneousCDSCurves(
    value_dt,
    libor_curve,
    cdsSpread3Y,
    cdsSpread5Y,
    cdsSpread7Y,
    cdsSpread10Y,
    num_credits,
):

    maturity3Y = value_dt.next_cds_date(36)
    maturity5Y = value_dt.next_cds_date(60)
    maturity7Y = value_dt.next_cds_date(84)
    maturity10Y = value_dt.next_cds_date(120)

    recovery_rate = 0.40

    cds3Y = CDS(value_dt, maturity3Y, cdsSpread3Y)
    cds5Y = CDS(value_dt, maturity5Y, cdsSpread5Y)
    cds7Y = CDS(value_dt, maturity7Y, cdsSpread7Y)
    cds10Y = CDS(value_dt, maturity10Y, cdsSpread10Y)

    contracts = [cds3Y, cds5Y, cds7Y, cds10Y]

    issuer_curve = CDSCurve(value_dt, contracts, libor_curve, recovery_rate)

    issuer_curves = []
    for _ in range(0, num_credits):
        issuer_curves.append(issuer_curve)

    return issuer_curves


def buildIborSingleCurve(value_dt, last_tenor="51Y"):

    settle_dt = value_dt.add_days(2)
    dc_type = DayCountTypes.ACT_360

    depos = []
    fras = []
    swaps = []

    maturity_dt = settle_dt.add_months(1)
    depo1 = IborDeposit(value_dt, maturity_dt, -0.00251, dc_type)
    depos.append(depo1)

    # Series of 1M futures
    start_dt = settle_dt.next_imm_date()
    end_dt = start_dt.add_months(1)
    fra = IborFRA(start_dt, end_dt, -0.0023, dc_type)
    fras.append(fra)

    start_dt = start_dt.add_months(1)
    end_dt = start_dt.add_months(1)
    fra = IborFRA(start_dt, end_dt, -0.00234, dc_type)
    fras.append(fra)

    start_dt = start_dt.add_months(1)
    end_dt = start_dt.add_months(1)
    fra = IborFRA(start_dt, end_dt, -0.00225, dc_type)
    fras.append(fra)

    start_dt = start_dt.add_months(1)
    end_dt = start_dt.add_months(1)
    fra = IborFRA(start_dt, end_dt, -0.00226, dc_type)
    fras.append(fra)

    start_dt = start_dt.add_months(1)
    end_dt = start_dt.add_months(1)
    fra = IborFRA(start_dt, end_dt, -0.00219, dc_type)
    fras.append(fra)

    start_dt = start_dt.add_months(1)
    end_dt = start_dt.add_months(1)
    fra = IborFRA(start_dt, end_dt, -0.00213, dc_type)
    fras.append(fra)

    start_dt = start_dt.add_months(1)
    end_dt = start_dt.add_months(1)
    fra = IborFRA(start_dt, end_dt, -0.00186, dc_type)
    fras.append(fra)

    start_dt = start_dt.add_months(1)
    end_dt = start_dt.add_months(1)
    fra = IborFRA(start_dt, end_dt, -0.00189, dc_type)
    fras.append(fra)

    start_dt = start_dt.add_months(1)
    end_dt = start_dt.add_months(1)
    fra = IborFRA(start_dt, end_dt, -0.00175, dc_type)
    fras.append(fra)

    start_dt = start_dt.add_months(1)
    end_dt = start_dt.add_months(1)
    fra = IborFRA(start_dt, end_dt, -0.00143, dc_type)
    fras.append(fra)

    start_dt = start_dt.add_months(1)
    end_dt = start_dt.add_months(1)
    fra = IborFRA(start_dt, end_dt, -0.00126, dc_type)
    fras.append(fra)

    start_dt = start_dt.add_months(1)
    end_dt = start_dt.add_months(1)
    fra = IborFRA(start_dt, end_dt, -0.00126, dc_type)
    fras.append(fra)

    ###########################################################################
    ###########################################################################
    ###########################################################################
    ###########################################################################

    fixed_freq = FrequencyTypes.ANNUAL
    dc_type = DayCountTypes.THIRTY_E_360
    fixed_leg_type = SwapTypes.PAY

    #######################################
    maturity_dt = settle_dt.add_months(24)
    swap_rate = -0.001506
    swap1 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap1)

    #######################################
    maturity_dt = settle_dt.add_months(36)
    swap_rate = -0.000185
    swap2 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap2)

    #######################################
    maturity_dt = settle_dt.add_months(48)
    swap_rate = 0.001358
    swap3 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap3)

    #######################################
    maturity_dt = settle_dt.add_months(60)
    swap_rate = 0.0027652
    swap4 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap4)

    #######################################
    maturity_dt = settle_dt.add_months(72)
    swap_rate = 0.0041539
    swap5 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap5)

    #######################################
    maturity_dt = settle_dt.add_months(84)
    swap_rate = 0.0054604
    swap6 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap6)

    #######################################
    maturity_dt = settle_dt.add_months(96)
    swap_rate = 0.006674
    swap7 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap7)

    #######################################
    maturity_dt = settle_dt.add_months(108)
    swap_rate = 0.007826
    swap8 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap8)

    #######################################
    maturity_dt = settle_dt.add_months(120)
    swap_rate = 0.008821
    swap9 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap9)

    #######################################
    maturity_dt = settle_dt.add_months(132)
    swap_rate = 0.0097379
    swap10 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap10)

    #######################################
    maturity_dt = settle_dt.add_months(144)
    swap_rate = 0.0105406
    swap11 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap11)

    #######################################
    maturity_dt = settle_dt.add_months(180)
    swap_rate = 0.0123927
    swap12 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap12)

    #######################################
    maturity_dt = settle_dt.add_months(240)
    swap_rate = 0.0139882
    swap13 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap13)

    #######################################
    maturity_dt = settle_dt.add_months(300)
    swap_rate = 0.0144972
    swap14 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap14)

    #######################################
    maturity_dt = settle_dt.add_months(360)
    swap_rate = 0.0146081
    swap15 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap15)

    #######################################
    maturity_dt = settle_dt.add_months(420)
    swap_rate = 0.01461897
    swap16 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap16)

    #######################################
    maturity_dt = settle_dt.add_months(480)
    swap_rate = 0.014567455
    swap17 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap17)

    #######################################
    maturity_dt = settle_dt.add_months(540)
    swap_rate = 0.0140826
    swap18 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap18)

    #######################################
    maturity_dt = settle_dt.add_months(600)
    swap_rate = 0.01436822
    swap19 = IborSwap(
        settle_dt, maturity_dt, fixed_leg_type, swap_rate, fixed_freq, dc_type
    )
    swaps.append(swap19)

    ########################################

    last_maturity_date = settle_dt.add_tenor(last_tenor)
    swaps_to_use = [
        s for s in swaps if s.fixed_leg.maturity_dt <= last_maturity_date
    ]

    libor_curve = IborSingleCurve(value_dt, depos, fras, swaps_to_use)

    return libor_curve
