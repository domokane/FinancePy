import numpy as np
import matplotlib.pyplot as plt

from financepy.market.curves.discount_curve_pwf_onf import DiscountCurvePWFONF
from financepy.utils.global_vars import G_BASIS_POINT
from financepy.utils.date import Date

# when set to True this file can be run standalone and will produce some useful output.
# Set to False to use as part of a testing framework
diagnostics_mode = False

########################################################################################


def test_fin_discount_curve_pcfonf_01():

    start_dt = Date(1, 1, 2015)
    knot_dts = [
        Date(1, 1, 2015),
        Date(1, 6, 2016),
        Date(1, 12, 2017),
        Date(1, 4, 2018),
        Date(1, 8, 2019),
    ]
    ondfwd_rates = [0, 0.02, 0.04, 0.06, 0.08]

    curve = DiscountCurvePWFONF(
        start_dt,
        knot_dts,
        ondfwd_rates,
    )

    test_dates = [
        Date(1, 6, 2015),
        Date(1, 2, 2017),
        Date(1, 2, 2018),
        Date(1, 8, 2018),
        Date(1, 12, 2019),
    ]
    expected_onfwd = [0.02, 0.04, 0.06, 0.08, 0.08]
    actual_onfwd = curve.fwd(test_dates)

    one_bp = 1e-4
    for d, e, a in zip(test_dates, expected_onfwd, actual_onfwd):
        assert (
            abs(e - a) < one_bp / 100
        ), f"Mismatch for date {d}, expected = {e}, actual = {a}"


########################################################################################


def test_fin_discount_curve_pcfonf_02():

    start_dt = Date(1, 1, 2015)
    knot_dts = [Date(1, 6, 2017), Date(1, 6, 2018), Date(2, 6, 2018)]
    ondfwd_rates = [0, 0.01, 0.0]

    curve = DiscountCurvePWFONF(
        start_dt,
        knot_dts,
        ondfwd_rates,
    )

    test_dates = [
        Date(1, 6, 2015),
        Date(1, 12, 2017),
        Date(1, 8, 2018),
    ]
    expected_onfwd = [0.0, 0.01, 0.0]
    actual_onfwd = curve.fwd(test_dates)

    for d, e, a in zip(test_dates, expected_onfwd, actual_onfwd):
        assert (
            abs(e - a) < G_BASIS_POINT / 100
        ), f"Mismatch for date {d}, expected = {e}, actual = {a}"


########################################################################################


def test_fin_discount_curve_pcfonf_03():

    valuation_date = Date(1, 1, 2015)
    start_dt = Date(1, 6, 2017)
    end_dt = Date(1, 6, 2018)
    level = 0.01

    curve = DiscountCurvePWFONF.brick_wall_curve(
        valuation_date,
        start_dt,
        end_dt,
        level,
    )

    test_dates = [
        Date(1, 6, 2015),
        Date(1, 12, 2017),
        Date(1, 8, 2018),
    ]
    expected_onfwd = [0.0, 0.01, 0.0]
    actual_onfwd = curve.fwd(test_dates)

    for d, e, a in zip(test_dates, expected_onfwd, actual_onfwd):
        assert (
            abs(e - a) < G_BASIS_POINT / 100
        ), f"Mismatch for date {d}, expected = {e}, actual = {a}"


########################################################################################


def test_fin_discount_curve_pcfonf_04():

    valuation_date = Date(1, 1, 2015)
    start_dt = valuation_date
    end_dt = Date(1, 6, 2018)
    level = 0.01

    curve = DiscountCurvePWFONF.brick_wall_curve(
        valuation_date,
        start_dt,
        end_dt,
        level,
    )

    test_dates = [
        Date(1, 6, 2015),
        Date(1, 12, 2017),
        Date(1, 8, 2018),
    ]
    expected_onfwd = [0.01, 0.01, 0.0]
    actual_onfwd = curve.fwd(test_dates)

    for d, e, a in zip(test_dates, expected_onfwd, actual_onfwd):
        assert (
            abs(e - a) < G_BASIS_POINT / 100
        ), f"Mismatch for date {d}, expected = {e}, actual = {a}"


########################################################################################


def test_fin_discount_curve_pcfonf_05():

    valuation_date = Date(1, 1, 2015)
    start_dt = valuation_date.add_months(-6)
    end_dt = Date(1, 6, 2018)
    level = 0.01

    curve = DiscountCurvePWFONF.brick_wall_curve(
        valuation_date,
        start_dt,
        end_dt,
        level,
    )

    test_dates = [
        Date(1, 6, 2015),
        Date(1, 12, 2017),
        Date(1, 8, 2018),
    ]
    expected_onfwd = [0.01, 0.01, 0.0]
    actual_onfwd = curve.fwd(test_dates)

    for d, e, a in zip(test_dates, expected_onfwd, actual_onfwd):
        assert (
            abs(e - a) < G_BASIS_POINT / 100
        ), f"Mismatch for date {d}, expected = {e}, actual = {a}"


########################################################################################


def test_fin_discount_curve_pcfonf_flat():

    valuation_date = Date(1, 1, 2015)
    level = 0.01
    curve = DiscountCurvePWFONF.flat_curve(
        valuation_date,
        level,
    )

    test_dates = [
        Date(1, 6, 2015),
        Date(1, 12, 2017),
        Date(1, 8, 2018),
    ]
    expected_onfwd = [0.01, 0.01, 0.01]
    actual_onfwd = curve.fwd(test_dates)

    if diagnostics_mode:
        years = np.linspace(0, 10, 100)
        dates = valuation_date.add_years(years)
        onrates = curve.fwd(dates)
        plt.plot(years, onrates)
        plt.show()

    for d, e, a in zip(test_dates, expected_onfwd, actual_onfwd):
        assert (
            abs(e - a) < G_BASIS_POINT / 100
        ), f"Mismatch for date {d}, expected = {e}, actual = {a}"


########################################################################################

if diagnostics_mode and __name__ == "__main__":
    test_fin_discount_curve_pcfonf_flat()
