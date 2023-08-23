###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from helpers import build_Ibor_Curve, buildFlatIssuerCurve
from financepy.utils.global_types import SwapTypes
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.credit.cds_curve import CDSCurve
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.products.rates.ibor_swap import IborSwap
from financepy.products.credit.cds import CDS
from financepy.products.credit.cds_index_option import CDSIndexOption
from financepy.products.credit.cds_index_portfolio import CDSIndexPortfolio
import os
import numpy as np


def test_dirty_priceCDSIndexOption():

    tradeDate = Date(1, 8, 2007)
    step_in_date = tradeDate.add_days(1)
    valuation_date = step_in_date

    libor_curve = build_Ibor_Curve(tradeDate)

    maturity3Y = tradeDate.next_cds_date(36)
    maturity5Y = tradeDate.next_cds_date(60)
    maturity7Y = tradeDate.next_cds_date(84)
    maturity10Y = tradeDate.next_cds_date(120)

    path = os.path.join(os.path.dirname(__file__),
                        './/data//CDX_NA_IG_S7_SPREADS.csv')
    f = open(path, 'r')
    data = f.readlines()
    f.close()
    issuer_curves = []

    for row in data[1:]:

        splitRow = row.split(",")
        creditName = splitRow[0]
        spd3Y = float(splitRow[1]) / 10000.0
        spd5Y = float(splitRow[2]) / 10000.0
        spd7Y = float(splitRow[3]) / 10000.0
        spd10Y = float(splitRow[4]) / 10000.0
        recovery_rate = float(splitRow[5])

        cds3Y = CDS(step_in_date, maturity3Y, spd3Y)
        cds5Y = CDS(step_in_date, maturity5Y, spd5Y)
        cds7Y = CDS(step_in_date, maturity7Y, spd7Y)
        cds10Y = CDS(step_in_date, maturity10Y, spd10Y)
        cds_contracts = [cds3Y, cds5Y, cds7Y, cds10Y]

        issuer_curve = CDSCurve(valuation_date,
                                cds_contracts,
                                libor_curve,
                                recovery_rate)

        issuer_curves.append(issuer_curve)

    ##########################################################################
    ##########################################################################

    indexUpfronts = [0.0, 0.0, 0.0, 0.0]
    indexMaturityDates = [Date(20, 12, 2009),
                          Date(20, 12, 2011),
                          Date(20, 12, 2013),
                          Date(20, 12, 2016)]
    indexRecovery = 0.40

    index_coupon = 0.004
    volatility = 0.50
    expiry_date = Date(1, 2, 2008)
    maturity_date = Date(20, 12, 2011)
    notional = 10000.0
    tolerance = 1e-6

    index_strike_results = [
        (20.0, 20.0, [16.1, 6.2, -70.7, 22.9, -60.6, 16.1, 6.1]),
        (25.0, 30.0, [11.9, 16.8, -35.3, 28.6, -40.3, 11.9, 16.7]),
        (50.0, 40.0, [63.6, 4.6, 0.0, 57.4, 60.5, 63.4, 4.6]),
    ]

    for index, strike, results in index_strike_results:

        #######################################################################

        cds_contracts = []
        for dt in indexMaturityDates:
            cds = CDS(valuation_date, dt, index / 10000.0)
            cds_contracts.append(cds)

        index_curve = CDSCurve(valuation_date, cds_contracts,
                               libor_curve, indexRecovery)

        indexSpreads = [index / 10000.0] * 4

        indexPortfolio = CDSIndexPortfolio()
        adjustedIssuerCurves = indexPortfolio.hazard_rate_adjust_intrinsic(
            valuation_date,
            issuer_curves,
            indexSpreads,
            indexUpfronts,
            indexMaturityDates,
            indexRecovery,
            tolerance)

        #######################################################################

        option = CDSIndexOption(expiry_date,
                                maturity_date,
                                index_coupon,
                                strike / 10000.0,
                                notional)

        v_pay_1, v_rec_1, strikeValue, mu, expH = option.value_anderson(
            valuation_date, adjustedIssuerCurves, indexRecovery, volatility)

        v_pay_2, v_rec_2 = option.value_adjusted_black(valuation_date,
                                                       index_curve,
                                                       indexRecovery,
                                                       libor_curve,
                                                       volatility)

        assert round(v_pay_1, 1) == results[0]
        assert round(v_rec_1, 1) == results[1]
        assert round(strikeValue, 1) == results[2]
        assert round(mu, 1) == results[3]
        assert round(expH, 1) == results[4]
        assert round(v_pay_2, 1) == results[5]
        assert round(v_rec_2, 1) == results[6]
