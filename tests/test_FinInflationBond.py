##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.market.curves.discount_curve_zeros import DiscountCurveZeros
from financepy.products.inflation.FinInflationIndexCurve import FinInflationIndexCurve
from financepy.products.bonds import YTMCalcType
from financepy.products.inflation.FinInflationBond import FinInflationBond
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes


def test_FinInflationBondBBG():
    # https://data.bloomberglp.com/bat/sites/3/2017/07/SF-2017_Paul-Fjeldsted.pdf
    # Look for CPI Bond example

    settlement_date = Date(21, 7, 2017)
    issue_date = Date(15, 7, 2010)
    maturity_date = Date(15, 7, 2020)
    coupon = 0.0125
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    face = 100.0
    baseCPIValue = 218.08532

    bond = FinInflationBond(issue_date,
                            maturity_date,
                            coupon,
                            freq_type,
                            accrual_type,
                            baseCPIValue)

    clean_price = 104.03502

    yld = bond.current_yield(clean_price)
    assert round(yld * 100, 4) == 1.2015

    # Inherited functions that just calculate real yield without CPI adjustments
    ytm = bond.yield_to_maturity(settlement_date,
                                 clean_price,
                                 YTMCalcType.UK_DMO)

    assert round(ytm, 4) == -0.0010

    ytm = bond.yield_to_maturity(settlement_date,
                                 clean_price,
                                 YTMCalcType.US_STREET)

    assert round(ytm, 4) == -0.0010

    ytm = bond.yield_to_maturity(settlement_date,
                                 clean_price,
                                 YTMCalcType.US_TREASURY)

    assert round(ytm, 4) == -0.0010

    dirty_price = bond.dirty_price_from_ytm(settlement_date, ytm)
    assert round(dirty_price, 4) == 104.0554

    clean_price = bond.clean_price_from_ytm(settlement_date, ytm)
    assert round(clean_price, 4) == 104.0350

    accddays = bond._accrued_days
    assert accddays == 6.0

    accrued_interest = bond._accrued_interest
    assert round(accrued_interest, 4) == 0.0204

    # Inflation functions that calculate nominal yield with CPI adjustment
    refCPIValue = 244.65884

    clean_price = bond.clean_price_from_ytm(settlement_date, ytm)
    assert round(clean_price, 4) == 104.0350

    inflationAccd = bond.calc_inflation_accrued_interest(settlement_date,
                                                         refCPIValue)

    assert round(inflationAccd * 100, 4) == 2.2864

    lastCpnCPIValue = 244.61839

    clean_price = bond.flat_price_from_yield_to_maturity(settlement_date, ytm,
                                                         lastCpnCPIValue,
                                                         YTMCalcType.US_TREASURY)

    assert round(clean_price, 4) == 116.6923

    principal = bond.inflation_principal(settlement_date,
                                         ytm,
                                         refCPIValue,
                                         YTMCalcType.US_TREASURY)

    assert round(principal, 4) == 116.7116

    duration = bond.dollar_duration(settlement_date, ytm)
    assert round(duration, 4) == 305.9297

    modified_duration = bond.modified_duration(settlement_date, ytm)
    assert round(modified_duration, 4) == 2.9401

    macauley_duration = bond.macauley_duration(settlement_date, ytm)
    assert round(macauley_duration, 4) == 2.9386

    conv = bond.convexity_from_ytm(settlement_date, ytm)
    assert round(conv, 4) == 0.1020


def test_FinInflationBondStack():
    # https://stackoverflow.com/questions/57676724/failing-to-obtain-correct-accrued-interest-with-quantlib-inflation-bond-pricer-i

    issue_date = Date(25, 9, 2013)
    maturity_date = Date(22, 3, 2068)
    coupon = 0.00125
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    face = 100.0
    baseCPIValue = 249.70

    bond = FinInflationBond(issue_date,
                            maturity_date,
                            coupon,
                            freq_type,
                            accrual_type,
                            face,
                            baseCPIValue)

    clean_price = 104.03502

    yld = bond.current_yield(clean_price)
    assert round(yld * 1000, 4) == 1.2015
