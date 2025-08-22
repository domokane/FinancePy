# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from FinTestCases import FinTestCases, global_test_case_mode
from financepy.utils.global_vars import G_DAYS_IN_YEARS
from financepy.market.curves.discount_curve_poly import DiscountCurvePoly
from financepy.market.curves.discount_curve_zeros import DiscountCurveZeros
from financepy.market.curves.discount_curve_pwl import DiscountCurvePWL
from financepy.market.curves.discount_curve_pwf import DiscountCurvePWF
from financepy.market.curves.discount_curve_pwf_onf import DiscountCurvePWFONF
from financepy.market.curves.discount_curve_nss import DiscountCurveNSS
from financepy.market.curves.discount_curve_ns import DiscountCurveNS
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.market.curves.interpolator import InterpTypes
from financepy.utils.date import Date, set_date_format, DateFormatTypes
import matplotlib.pyplot as plt
import numpy as np

import sys

sys.path.append("..")


test_cases = FinTestCases(__file__, global_test_case_mode)

set_date_format(DateFormatTypes.UK_LONG)

plot_graphs = False

# TODO: Add other discount discount


########################################################################################


def test_FinDiscountCurves():

    # Create a curve from times and discount factors
    value_dt = Date(1, 1, 2018)
    years = [1.0, 2.0, 3.0, 4.0, 5.0]
    dates = value_dt.add_years(years)
    years2 = []

    for dt in dates:
        y = (dt - value_dt) / G_DAYS_IN_YEARS
        years2.append(y)

    rates = np.array([0.05, 0.06, 0.065, 0.07, 0.075])
    discount_factors = np.exp(-np.array(rates) * np.array(years2))
    curves_list = []

    fin_discount_curve = DiscountCurve(
        value_dt, dates, discount_factors, InterpTypes.FLAT_FWD_RATES
    )
    curves_list.append(fin_discount_curve)

    fin_discount_curve_flat = DiscountCurveFlat(value_dt, 0.05)
    curves_list.append(fin_discount_curve_flat)

    fin_discount_curve_ns = DiscountCurveNS(value_dt, 0.0305, -0.01, 0.08, 10.0)
    curves_list.append(fin_discount_curve_ns)

    fin_discount_curve_nss = DiscountCurveNSS(
        value_dt, 0.035, -0.02, 0.09, 0.1, 1.0, 2.0
    )
    curves_list.append(fin_discount_curve_nss)

    fin_discount_curve_poly = DiscountCurvePoly(value_dt, [0.05, 0.002, -0.00005])
    curves_list.append(fin_discount_curve_poly)

    fin_discount_curve_PWF = DiscountCurvePWF(value_dt, dates, rates)
    curves_list.append(fin_discount_curve_PWF)

    fin_discount_curve_PWL = DiscountCurvePWL(value_dt, dates, rates)
    curves_list.append(fin_discount_curve_PWL)

    fin_discount_curve_zeros = DiscountCurveZeros(value_dt, dates, rates)
    curves_list.append(fin_discount_curve_zeros)

    fin_discount_curve_PWFONF = DiscountCurvePWFONF(value_dt, dates, rates)
    curves_list.append(fin_discount_curve_PWFONF)

    curve_names = []
    for curve in curves_list:
        curve_names.append(type(curve).__name__)

    test_cases.banner("SINGLE CALLS NO VECTORS")
    test_cases.header("CURVE", "DATE", "ZERO", "DF", "CCFWD", "MMFWD", "SWAP")

    years = np.linspace(1, 10, 10)
    fwd_maturity_dts = value_dt.add_years(years)

    test_cases.banner("######################################################")
    test_cases.banner("SINGLE CALLS")
    test_cases.banner("######################################################")

    for name, curve in zip(curve_names, curves_list):
        for fwd_maturity_dt in fwd_maturity_dts:
            tenor = "3M"
            zero_rate = curve.zero_rate(fwd_maturity_dt)
            fwd = curve.fwd(fwd_maturity_dt)
            fwd_rate = curve.fwd_rate(fwd_maturity_dt, tenor)
            swap_rate = curve.swap_rate(value_dt, fwd_maturity_dt)
            df = curve.df(fwd_maturity_dt)

            test_cases.print(
                "%-20s" % name,
                "%-12s" % fwd_maturity_dt,
                "%7.6f" % (zero_rate),
                "%8.7f" % (df),
                "%7.6f" % (fwd),
                "%7.6f" % (fwd_rate),
                "%7.6f" % (swap_rate),
            )

    # Examine vectorisation
    test_cases.banner("######################################################")
    test_cases.banner("VECTORISATIONS")
    test_cases.banner("######################################################")

    for name, curve in zip(curve_names, curves_list):
        tenor = "3M"
        zero_rate = curve.zero_rate(fwd_maturity_dts)
        fwd = curve.fwd(fwd_maturity_dts)
        fwd_rate = curve.fwd_rate(fwd_maturity_dts, tenor)
        swap_rate = curve.swap_rate(value_dt, fwd_maturity_dts)
        df = curve.df(fwd_maturity_dts)

        for i in range(0, len(fwd_maturity_dts)):
            test_cases.print(
                "%-20s" % name,
                "%-12s" % fwd_maturity_dts[i],
                "%7.6f" % (zero_rate[i]),
                "%8.7f" % (df[i]),
                "%7.6f" % (fwd[i]),
                "%7.6f" % (fwd_rate[i]),
                "%7.6f" % (swap_rate[i]),
            )

    if plot_graphs:

        years = np.linspace(0, 10, 121)
        years2 = years + 1.0
        fwd_dts = value_dt.add_years(years)
        fwd_dts_2 = value_dt.add_years(years2)

        plt.figure()
        for name, curve in zip(curve_names, curves_list):
            zero_rates = curve.zero_rate(fwd_dts)
            plt.plot(years, zero_rates, label=name)
        plt.legend()
        plt.title("Zero Rates")

        plt.figure()
        for name, curve in zip(curve_names, curves_list):
            fwd_rates = curve.fwd(fwd_dts)
            plt.plot(years, fwd_rates, label=name)
        plt.legend()
        plt.title("CC Fwd Rates")

        plt.figure()
        for name, curve in zip(curve_names, curves_list):
            fwd_rates = curve.fwd_rate(fwd_dts, fwd_dts_2)
            plt.plot(years, fwd_rates, label=name)
        plt.legend()
        plt.title("CC Fwd Rates")

        plt.figure()
        for name, curve in zip(curve_names, curves_list):
            dfs = curve.df(fwd_dts)
            plt.plot(years, dfs, label=name)
        plt.legend()
        plt.title("Discount Factors")


########################################################################################


test_FinDiscountCurves()
test_cases.compare_test_cases()

########################################################################################

