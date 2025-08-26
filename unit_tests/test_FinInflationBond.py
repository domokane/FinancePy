# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.market.curves.discount_curve_zeros import DiscountCurveZeros
from financepy.products.inflation.FinInflationIndexCurve import FinInflationIndexCurve
from financepy.products.bonds import YTMCalcType
from financepy.products.inflation.FinInflationBond import FinInflationBond
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes

########################################################################################


def test_fin_inflation_bond_bbg():

    # https://data.bloomberglp.com/bat/sites/3/2017/07/SF-2017_Paul-Fjeldsted.pdf
    # Look for CPI Bond example

    settle_dt = Date(21, 7, 2017)
    issue_dt = Date(15, 7, 2010)
    maturity_dt = Date(15, 7, 2020)
    coupon = 0.0125
    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    base_cpi_value = 218.08532
    ex_dividend_days = 0

    bond = FinInflationBond(
        issue_dt,
        maturity_dt,
        coupon,
        freq_type,
        dc_type,
        ex_dividend_days,
        base_cpi_value,
    )

    clean_price = 104.03502

    yld = bond.current_yield(clean_price)
    assert round(yld * 100, 4) == 1.2015

    # Inherited functions that just calculate real yield without CPI adjustments
    ytm = bond.yield_to_maturity(settle_dt, clean_price, YTMCalcType.UK_DMO)

    assert round(ytm, 4) == -0.0010

    ytm = bond.yield_to_maturity(settle_dt, clean_price, YTMCalcType.US_STREET)

    assert round(ytm, 4) == -0.0010

    ytm = bond.yield_to_maturity(settle_dt, clean_price, YTMCalcType.US_TREASURY)

    assert round(ytm, 4) == -0.0010

    dirty_price = bond.dirty_price_from_ytm(settle_dt, ytm)
    assert round(dirty_price, 4) == 104.0554

    clean_price = bond.clean_price_from_ytm(settle_dt, ytm)
    assert round(clean_price, 4) == 104.0350

    accddays = bond.accrued_days
    assert accddays == 6.0

    accrued_interest = bond.accrued_int
    assert round(accrued_interest, 4) == 0.0204

    # Inflation functions that calculate nominal yield with CPI adjustment
    ref_cpi_value = 244.65884

    clean_price = bond.clean_price_from_ytm(settle_dt, ytm)
    assert round(clean_price, 4) == 104.0350

    face = 100.0
    inflation_accd = bond.inflation_accrued_interest(settle_dt, face, ref_cpi_value)

    assert round(inflation_accd * 100, 4) == 2.2864

    last_cpn_cpi_value = 244.61839

    clean_price = bond.flat_price_from_yield_to_maturity(
        settle_dt, ytm, last_cpn_cpi_value, YTMCalcType.US_TREASURY
    )

    assert round(clean_price, 4) == 116.6923

    principal = bond.inflation_principal(
        settle_dt, face, ytm, ref_cpi_value, YTMCalcType.US_TREASURY
    )

    assert round(principal, 4) == 116.7342

    duration = bond.dollar_duration(settle_dt, ytm)
    assert round(duration, 4) == 305.9297

    modified_duration = bond.modified_duration(settle_dt, ytm)
    assert round(modified_duration, 4) == 2.9401

    macauley_duration = bond.macauley_duration(settle_dt, ytm)
    assert round(macauley_duration, 4) == 2.9386

    conv = bond.convexity_from_ytm(settle_dt, ytm)
    assert round(conv, 4) == 0.1020


########################################################################################


def test_fin_inflation_bond_stack():

    # https://stackoverflow.com/questions/57676724/failing-to-obtain-correct-accrued-interest-with-quantlib-inflation-bond-pricer-i

    issue_dt = Date(25, 9, 2013)
    maturity_dt = Date(22, 3, 2068)
    coupon = 0.00125
    freq_type = FrequencyTypes.SEMI_ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ICMA
    base_cpi_value = 249.70
    ex_dividend_days = 0

    bond = FinInflationBond(
        issue_dt,
        maturity_dt,
        coupon,
        freq_type,
        dc_type,
        ex_dividend_days,
        base_cpi_value,
    )

    clean_price = 104.03502

    yld = bond.current_yield(clean_price)
    assert round(yld * 1000, 4) == 1.2015
