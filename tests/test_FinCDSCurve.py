###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.global_types import SwapTypes
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.products.credit.cds_curve import CDSCurve
from financepy.products.rates.ibor_swap import IborSwap
from financepy.products.credit.cds import CDS


def test_FinCDSCurve():

    curve_date = Date(20, 12, 2018)

    swaps = []
    depos = []
    fras = []

    fixedDCC = DayCountTypes.ACT_365F
    fixedFreq = FrequencyTypes.SEMI_ANNUAL
    fixed_coupon = 0.05

    for i in range(1, 11):

        maturity_date = curve_date.add_months(12 * i)
        swap = IborSwap(curve_date,
                        maturity_date,
                        SwapTypes.PAY,
                        fixed_coupon,
                        fixedFreq,
                        fixedDCC)
        swaps.append(swap)

    libor_curve = IborSingleCurve(curve_date, depos, fras, swaps)

    cds_contracts = []

    for i in range(1, 11):
        maturity_date = curve_date.add_months(12 * i)
        cds = CDS(curve_date, maturity_date, 0.005 + 0.001 * (i - 1))
        cds_contracts.append(cds)

    recovery_rate = 0.40

    issuer_curve = CDSCurve(curve_date,
                            cds_contracts,
                            libor_curve,
                            recovery_rate,
                            use_cache=False)

    assert round(issuer_curve._times[0], 4) == 0.0
    assert round(issuer_curve._times[5], 4) == 5.0027
    assert round(issuer_curve._times[9], 4) == 9.0055
    assert round(issuer_curve._values[0], 4) == 1.0
    assert round(issuer_curve._values[5], 4) == 0.9249
    assert round(issuer_curve._values[9], 4) == 0.8071

    i = 1
    maturity_date = curve_date.add_months(12 * i)
    cds = CDS(curve_date, maturity_date, 0.005 + 0.001 * (i - 1))
    v = cds.value(curve_date, issuer_curve, recovery_rate)
    assert round(v['dirty_pv'] * 1000, 4) == 5.6028
    assert round(v['clean_pv'] * 1000, 4) == 5.6028

    i = 5
    maturity_date = curve_date.add_months(12 * i)
    cds = CDS(curve_date, maturity_date, 0.005 + 0.001 * (i - 1))
    v = cds.value(curve_date, issuer_curve, recovery_rate)
    assert round(v['dirty_pv'] * 1000, 4) == 8.3480
    assert round(v['clean_pv'] * 1000, 4) == 8.3480

    i = 10
    maturity_date = curve_date.add_months(12 * i)
    cds = CDS(curve_date, maturity_date, 0.005 + 0.001 * (i - 1))
    v = cds.value(curve_date, issuer_curve, recovery_rate)
    assert round(v['dirty_pv'] * 1000, 4) == -1.3178
    assert round(v['clean_pv'] * 1000, 4) == -1.3178

