###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from helpers import build_Ibor_Curve
from financepy.utils.date import Date
from financepy.products.credit.cds_curve import CDSCurve
from financepy.products.credit.cds import CDS
from financepy.products.credit.cds_index_portfolio import CDSIndexPortfolio
import os


def test_CDSIndexPortfolio():
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

    # Now determine the average spread of the index
    cdsIndex = CDSIndexPortfolio()

    averageSpd3Y = cdsIndex.average_spread(valuation_date,
                                           step_in_date,
                                           maturity3Y,
                                           issuer_curves) * 10000.0

    averageSpd5Y = cdsIndex.average_spread(valuation_date,
                                           step_in_date,
                                           maturity5Y,
                                           issuer_curves) * 10000.0

    averageSpd7Y = cdsIndex.average_spread(valuation_date,
                                           step_in_date,
                                           maturity7Y,
                                           issuer_curves) * 10000.0

    averageSpd10Y = cdsIndex.average_spread(valuation_date,
                                            step_in_date,
                                            maturity10Y,
                                            issuer_curves) * 10000.0

    assert round(averageSpd3Y, 4) == 19.8221
    assert round(averageSpd5Y, 4) == 36.0357
    assert round(averageSpd7Y, 4) == 50.1336
    assert round(averageSpd10Y, 4) == 63.6622

    # Now determine the intrinsic spread of the index to the same maturity
    # dates. As the single name CDS contracts
    cdsIndex = CDSIndexPortfolio()

    intrinsicSpd3Y = cdsIndex.intrinsic_spread(valuation_date,
                                               step_in_date,
                                               maturity3Y,
                                               issuer_curves) * 10000.0

    intrinsicSpd5Y = cdsIndex.intrinsic_spread(valuation_date,
                                               step_in_date,
                                               maturity5Y,
                                               issuer_curves) * 10000.0

    intrinsicSpd7Y = cdsIndex.intrinsic_spread(valuation_date,
                                               step_in_date,
                                               maturity7Y,
                                               issuer_curves) * 10000.0

    intrinsicSpd10Y = cdsIndex.intrinsic_spread(valuation_date,
                                                step_in_date,
                                                maturity10Y,
                                                issuer_curves) * 10000.0

    assert round(intrinsicSpd3Y, 4) == 19.6789
    assert round(intrinsicSpd5Y, 4) == 35.5393
    assert round(intrinsicSpd7Y, 4) == 49.0120
    assert round(intrinsicSpd10Y, 4) == 61.4139
