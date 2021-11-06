###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.math import scale
from financepy.market.curves.discount_curve_ns import DiscountCurveNS
from financepy.utils.date import Date
import numpy as np


tau = 2.0
times = np.linspace(0.0, 10.0, 5)
curve_date = Date(6, 6, 2019)
dates = curve_date.add_years(times)


def test_factor_loading_zero_rates():
    curve = DiscountCurveNS(curve_date, 1, 0, 0, tau)
    factor1loading = curve.zero_rate(dates)
    curve = DiscountCurveNS(curve_date, 0, 1, 0, tau)
    factor2loading = curve.zero_rate(dates)
    curve = DiscountCurveNS(curve_date, 0, 0, 1, tau)
    factor3loading = curve.zero_rate(dates)

    assert [round(x, 4) for x in factor1loading] == [
        1.0, 0.9852, 0.9855, 0.9856, 0.9855]
    assert [round(x, 4) for x in factor2loading] == [
        1.0001, 0.5622, 0.3618,  0.2566, 0.1958]
    assert [round(x, 4) for x in factor3loading] == [
        0.0001, 0.2801, 0.2809, 0.2334, 0.1891]


def test_beta1():
    beta2 = -0.02
    beta3 = 0.02

    beta1 = 0.03
    curve = DiscountCurveNS(curve_date, beta1, beta2, beta3, tau)
    zero_rates = curve.zero_rate(dates)
    assert [round(x, 4) for x in zero_rates] == [
        0.01, 0.0239, 0.0279, 0.0291, 0.0294]

    beta1 = 0.04
    curve = DiscountCurveNS(curve_date, beta1, beta2, beta3, tau)
    zero_rates = curve.zero_rate(dates)
    assert [round(x, 4) for x in zero_rates] == [
        0.02, 0.0338, 0.0378, 0.039, 0.0393]

    beta1 = 0.05
    curve = DiscountCurveNS(curve_date, beta1, beta2, beta3, tau)
    zero_rates = curve.zero_rate(dates)
    assert [round(x, 4) for x in zero_rates] == [
        0.03, 0.0436, 0.0477, 0.0488, 0.0491]

    beta1 = 0.06
    curve = DiscountCurveNS(curve_date, beta1, beta2, beta3, tau)
    zero_rates = curve.zero_rate(dates)
    assert [round(x, 4) for x in zero_rates] == [
        0.04, 0.0535, 0.0575, 0.0587, 0.059]

    beta1 = 0.07
    curve = DiscountCurveNS(curve_date, beta1, beta2, beta3, tau)
    zero_rates = curve.zero_rate(dates)
    assert [round(x, 4) for x in zero_rates] == [
        0.05, 0.0633, 0.0674, 0.0685, 0.0689]


def test_beta2():
    beta1 = 0.06
    beta3 = 0.02

    beta2 = -0.04
    curve = DiscountCurveNS(curve_date, beta1, beta2, beta3, tau)
    zero_rates = curve.zero_rate(dates)
    assert [round(x, 4) for x in zero_rates] == [
        0.02, 0.0422, 0.0503, 0.0535, 0.0551]

    beta2 = -0.02
    curve = DiscountCurveNS(curve_date, beta1, beta2, beta3, tau)
    zero_rates = curve.zero_rate(dates)
    assert [round(x, 4) for x in zero_rates] == [
        0.04, 0.0535, 0.0575, 0.0587, 0.059]

    beta2 = 0.00
    curve = DiscountCurveNS(curve_date, beta1, beta2, beta3, tau)
    zero_rates = curve.zero_rate(dates)
    assert [round(x, 4) for x in zero_rates] == [
        0.06, 0.0647, 0.0648, 0.0638, 0.0629]

    beta2 = 0.02
    curve = DiscountCurveNS(curve_date, beta1, beta2, beta3, tau)
    zero_rates = curve.zero_rate(dates)
    assert [round(x, 4) for x in zero_rates] == [
        0.08, 0.076, 0.072, 0.0689, 0.0668]

    beta2 = 0.04
    curve = DiscountCurveNS(curve_date, beta1, beta2, beta3, tau)
    zero_rates = curve.zero_rate(dates)
    assert [round(x, 4) for x in zero_rates] == [
        0.1, 0.0872, 0.0792, 0.0741, 0.0707]


def test_beta3():
    beta1 = 0.06
    beta2 = -0.02

    beta3 = -0.02
    curve = DiscountCurveNS(curve_date, beta1, beta2, beta3, tau)
    zero_rates = curve.zero_rate(dates)
    assert [round(x, 4) for x in zero_rates] == [
        0.04, 0.0423, 0.0463, 0.0493, 0.0514]

    beta3 = 0.00
    curve = DiscountCurveNS(curve_date, beta1, beta2, beta3, tau)
    zero_rates = curve.zero_rate(dates)
    assert [round(x, 4) for x in zero_rates] == [
        0.04, 0.0479, 0.0519, 0.054, 0.0552]

    beta3 = 0.02
    curve = DiscountCurveNS(curve_date, beta1, beta2, beta3, tau)
    zero_rates = curve.zero_rate(dates)
    assert [round(x, 4) for x in zero_rates] == [
        0.04, 0.0535, 0.0575, 0.0587, 0.059]

    beta3 = 0.04
    curve = DiscountCurveNS(curve_date, beta1, beta2, beta3, tau)
    zero_rates = curve.zero_rate(dates)
    assert [round(x, 4) for x in zero_rates] == [
        0.04, 0.0591, 0.0631, 0.0633, 0.0628]

    beta3 = 0.06
    curve = DiscountCurveNS(curve_date, beta1, beta2, beta3, tau)
    zero_rates = curve.zero_rate(dates)
    assert [round(x, 4) for x in zero_rates] == [
        0.04, 0.0647, 0.0688, 0.068, 0.0666]
