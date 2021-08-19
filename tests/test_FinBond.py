##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date
from financepy.utils.math import ONE_MILLION
from financepy.products.bonds.bond import Bond
from financepy.products.bonds.bond import YTMCalcType


def test_bondtutor_example():
    #  EXAMPLE FROM http://bondtutor.com/btchp4/topic6/topic6.htm

    accrualConvention = DayCountTypes.ACT_ACT_ICMA
    y = 0.062267
    settlement_date = Date(19, 4, 1994)
    issue_date = Date(15, 7, 1990)
    maturity_date = Date(15, 7, 1997)
    coupon = 0.085
    face = ONE_MILLION
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bond = Bond(issue_date, maturity_date,
                coupon, freq_type, accrualConvention, face)

    full_price = bond.full_price_from_ytm(settlement_date, y)
    assert round(full_price, 4) == 108.7696
    clean_price = bond.clean_price_from_ytm(settlement_date, y)
    assert round(clean_price, 4) == 106.5625
    accrued_interest = bond._accrued_interest
    assert round(accrued_interest, 4) == 22071.8232
    ytm = bond.yield_to_maturity(settlement_date, clean_price)
    assert round(ytm, 4) == 0.0622

    bump = 1e-4
    priceBumpedUp = bond.full_price_from_ytm(settlement_date, y + bump)
    assert round(priceBumpedUp, 4) == 108.7395

    priceBumpedDn = bond.full_price_from_ytm(settlement_date, y - bump)
    assert round(priceBumpedDn, 4) == 108.7998

    durationByBump = -(priceBumpedUp - full_price) / bump
    assert round(durationByBump, 4) == 301.1932

    duration = bond.dollar_duration(settlement_date, y)
    assert round(duration, 4) == 301.2458
    assert round(duration - durationByBump, 4) == 0.0526

    modified_duration = bond.modified_duration(settlement_date, y)
    assert round(modified_duration, 4) == 2.7696

    macauley_duration = bond.macauley_duration(settlement_date, y)
    assert round(macauley_duration, 4) == 2.8558

    conv = bond.convexity_from_ytm(settlement_date, y)
    assert round(conv, 4) == 0.0967


def test_bloomberg_us_treasury_example():
    settlement_date = Date(21, 7, 2017)
    issue_date = Date(15, 5, 2010)
    maturity_date = Date(15, 5, 2027)
    coupon = 0.02375
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA
    face = 100.0

    bond = Bond(issue_date,
                maturity_date,
                coupon,
                freq_type,
                accrual_type,
                face)

    clean_price = 99.7808417

    yld = bond.current_yield(clean_price)
    assert round(yld, 4) == 0.0238

    ytm = bond.yield_to_maturity(settlement_date, clean_price,
                                 YTMCalcType.UK_DMO)
    assert round(ytm, 4) == 0.0240

    ytm = bond.yield_to_maturity(settlement_date, clean_price,
                                 YTMCalcType.US_STREET)
    assert round(ytm, 4) == 0.0240

    ytm = bond.yield_to_maturity(settlement_date, clean_price,
                                 YTMCalcType.US_TREASURY)
    assert round(ytm, 4) == 0.0240

    full_price = bond.full_price_from_ytm(settlement_date, ytm)
    assert round(full_price, 4) == 100.2149

    clean_price = bond.clean_price_from_ytm(settlement_date, ytm)
    assert round(clean_price, 4) == 99.7825

    accrued_interest = bond._accrued_interest
    assert round(accrued_interest, 4) == 0.4324

    accddays = bond._accrued_days
    assert round(accddays, 4) == 67.0

    duration = bond.dollar_duration(settlement_date, ytm)
    assert round(duration, 4) == 869.0934

    modified_duration = bond.modified_duration(settlement_date, ytm)
    assert round(modified_duration, 4) == 8.6723

    macauley_duration = bond.macauley_duration(settlement_date, ytm)
    assert round(macauley_duration, 4) == 8.7764

    conv = bond.convexity_from_ytm(settlement_date, ytm)
    assert round(conv, 4) == 0.8517


def test_bloomberg_apple_corp_example():
    settlement_date = Date(21, 7, 2017)
    issue_date = Date(13, 5, 2012)
    maturity_date = Date(13, 5, 2022)
    coupon = 0.027
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.THIRTY_E_360_ISDA
    face = 100.0

    bond = Bond(issue_date, maturity_date,
                coupon, freq_type, accrual_type, face)

    clean_price = 101.581564

    yld = bond.current_yield(clean_price)
    assert round(yld, 4) == 0.0266

    ytm = bond.yield_to_maturity(settlement_date, clean_price,
                                 YTMCalcType.UK_DMO)
    assert round(ytm, 4) == 0.0235

    ytm = bond.yield_to_maturity(settlement_date, clean_price,
                                 YTMCalcType.US_STREET)
    assert round(ytm, 4) == 0.0235

    ytm = bond.yield_to_maturity(settlement_date, clean_price,
                                 YTMCalcType.US_TREASURY)
    assert round(ytm, 4) == 0.0235

    full_price = bond.full_price_from_ytm(settlement_date, ytm)
    assert round(full_price, 4) == 102.0932

    clean_price = bond.clean_price_from_ytm(settlement_date, ytm)
    assert round(clean_price, 4) == 101.5832

    accddays = bond._accrued_days
    assert accddays == 68

    accrued_interest = bond._accrued_interest
    assert round(accrued_interest, 4) == 0.51

    duration = bond.dollar_duration(settlement_date, ytm)
    assert round(duration, 4) == 456.5778

    modified_duration = bond.modified_duration(settlement_date, ytm)
    assert round(modified_duration, 4) == 4.4722

    macauley_duration = bond.macauley_duration(settlement_date, ytm)
    assert round(macauley_duration, 4) == 4.5247

    conv = bond.convexity_from_ytm(settlement_date, ytm)
    assert round(conv, 4) == 0.2302
