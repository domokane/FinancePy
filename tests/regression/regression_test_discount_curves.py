# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import matplotlib.pyplot as plt
import numpy as np

import add_fp_to_path

from financepy.utils.global_vars import G_DAYS_IN_YEAR
from financepy.market.curves.poly_discount_curve import PolyDiscountCurve
from financepy.market.curves.zero_rates_discount_curve import ZeroRatesDiscountCurve
from financepy.market.curves.pwl_discount_curve import PWLDiscountCurve
from financepy.market.curves.pwf_discount_curve import PWFDiscountCurve
from financepy.market.curves.pwf_onf_discount_curve import PWFONFDiscountCurve
from financepy.market.curves.nss_discount_curve import NSSDiscountCurve
from financepy.market.curves.ns_discount_curve import NSDiscountCurve
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.market.curves.interpolator import InterpTypes
from financepy.utils.date_format import set_date_format, DateFormatTypes
from financepy.utils.date import Date
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.compounding import CompoundingTypes
from financepy.utils.day_count import DayCountTypes

from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

PLOT_GRAPHS = False

# TODO: Add other discount discount

########################################################################################


def test_fin_discount_curves():

    # Create a curve from times and discount factors
    value_dt = Date(1, 1, 2018)
    years = [1.0, 2.0, 3.0, 4.0, 5.0]
    dates = value_dt.add_years(years)
    years2 = []

    for dt in dates:
        y = (dt - value_dt) / G_DAYS_IN_YEAR
        years2.append(y)

    rates = np.array([0.05, 0.06, 0.065, 0.07, 0.075])
    discount_factors = np.exp(-np.array(rates) * np.array(years2))
    curves_list = []

    fin_discount_curve = DiscountCurve(
        value_dt, dates, discount_factors, InterpTypes.FLAT_FWD_RATES
    )
    curves_list.append(fin_discount_curve)

    fin_flat_discount_curve = FlatDiscountCurve(value_dt, 0.05)
    curves_list.append(fin_flat_discount_curve)

    fin_discount_curve_ns = NSDiscountCurve(
        value_dt, 0.0305, -0.01, 0.08, 10.0
    )
    curves_list.append(fin_discount_curve_ns)

    fin_discount_curve_nss = NSSDiscountCurve(
        value_dt, 0.035, -0.02, 0.09, 0.1, 1.0, 2.0
    )
    curves_list.append(fin_discount_curve_nss)

    fin_discount_curve_poly = PolyDiscountCurve(
        value_dt, [0.05, 0.002, -0.00005]
    )
    curves_list.append(fin_discount_curve_poly)

    fin_discount_curve_pwf = PWFDiscountCurve(value_dt, dates, rates)
    curves_list.append(fin_discount_curve_pwf)

    fin_discount_curve_pwl = PWLDiscountCurve(value_dt, dates, rates)
    curves_list.append(fin_discount_curve_pwl)

    fin_discount_curve_zeros = ZeroRatesDiscountCurve(value_dt, dates, rates)
    curves_list.append(fin_discount_curve_zeros)

    fin_discount_curve_pwfonf = PWFONFDiscountCurve(value_dt, dates, rates)
    curves_list.append(fin_discount_curve_pwfonf)

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

    accrual_dc_type = DayCountTypes.ACT_360

    for name, curve in zip(curve_names, curves_list):

        #        print(name)

        for fwd_maturity_dt in fwd_maturity_dts:

            tenor = "3M"
            zero_rate = curve.zero_rate(fwd_maturity_dt)
            cc_fwd = curve.fwd_rate(
                fwd_maturity_dt,
                tenor,
                accrual_dc_type,
                CompoundingTypes.CONTINUOUS,
            )
            mm_fwd = curve.fwd_rate(fwd_maturity_dt, tenor)
            swap_rate = curve.swap_rate(value_dt, fwd_maturity_dt)
            df = curve.df(fwd_maturity_dt)

            test_cases.print(
                "%-20s" % name,
                "%-12s" % fwd_maturity_dt,
                "%7.6f" % (zero_rate),
                "%8.7f" % (df),
                "%7.6f" % (cc_fwd),
                "%7.6f" % (mm_fwd),
                "%7.6f" % (swap_rate),
            )

        # print(curve)
        # bumped_curve = curve.bump_parallel(0.0001)
        # print("===============>>>>BUMPED")
        # print(bumped_curve)
        # print("====================================================")
        # print("====================================================")

    # Examine vectorisation
    test_cases.banner("######################################################")
    test_cases.banner("VECTORISATIONS")
    test_cases.banner("######################################################")

    for name, curve in zip(curve_names, curves_list):
        tenor = "3M"
        zero_rate = curve.zero_rate(fwd_maturity_dts)
        fwd = curve.fwd_rate_inst(fwd_maturity_dts)
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

    if PLOT_GRAPHS:

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
            fwd_rates = curve.fwd_rate_inst(fwd_dts)
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


def test_bump():

    # Create a curve from times and discount factors
    value_dt = Date(1, 1, 2018)
    years = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    dates = value_dt.add_years(years)
    years2 = []

    for dt in dates:
        y = (dt - value_dt) / G_DAYS_IN_YEAR
        years2.append(y)

    zeros = np.array(
        [0.05, 0.06, 0.065, 0.07, 0.075, 0.08, 0.081, 0.082, 0.083, 0.084]
    )
    dfs = np.exp(-np.array(zeros) * np.array(years2))

    curve = DiscountCurve(value_dt, dates, dfs, InterpTypes.FLAT_FWD_RATES)

    #    print(curve)
    #    print(zeros)

    bp = 0.0001
    bumped_curve = curve.bump_parallel(bp)
    bumped_zeros = bumped_curve.zero_rate(dates)


#    print(bumped_curve)
#    print(bumped_zeros)


########################################################################################


test_fin_discount_curves()
test_bump()
test_cases.compare_test_cases()
