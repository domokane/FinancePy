###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.math import scale
from financepy.market.curves.discount_curve_nss import DiscountCurveNSS
from financepy.utils.date import Date
import numpy as np


tau1 = 2.0
tau2 = 0.5
times = np.linspace(0.0, 10.0, 5)
start_date = Date(1, 1, 2020)
dates = start_date.add_years(times)


def test_factor_loading_zero_rates():
    curve1 = DiscountCurveNSS(start_date, 1., 0., 0., 0., tau1, tau2)
    factor1loading = curve1.zero_rate(dates)
    curve2 = DiscountCurveNSS(start_date, 0., 1., 0., 0., tau1, tau2)
    factor2loading = curve2.zero_rate(dates)
    curve3 = DiscountCurveNSS(start_date, 0., 0., 1., 0., tau1, tau2)
    factor3loading = curve3.zero_rate(dates)
    curve4 = DiscountCurveNSS(start_date, 0., 0., 0., 1., tau1, tau2)
    factor4loading = curve4.zero_rate(dates)

    assert [round(x, 4) for x in factor1loading] == [
        1.0, 0.9852, 0.9852, 0.9856, 0.9855]
    assert [round(x, 4) for x in factor2loading] == [
        1.0001, 0.5628, 0.3617, 0.2568, 0.1958]
    assert [round(x, 4) for x in factor3loading] == [.0001,
                                                     0.2800, 0.2809, 0.2335, 0.1891]
    assert [round(x, 4) for x in factor4loading] == [-0.0,
                                                     0.1893, 0.0985, 0.0657, 0.0493]


def test_beta1():
    beta2 = -0.02
    beta3 = -0.02
    beta4 = 0.08

    beta1 = 0.03
    curve1 = DiscountCurveNSS(start_date,
                              beta1, beta2, beta3, beta4, tau1, tau2)
    zero_rates = curve1.zero_rate(dates)
    assert [round(x, 4) for x in zero_rates] == [
        0.01, 0.0278, 0.0246, 0.025, 0.0258]

    beta1 = 0.04
    curve2 = DiscountCurveNSS(start_date,
                              beta1, beta2, beta3, beta4, tau1, tau2)
    zero_rates = curve2.zero_rate(dates)
    assert [round(x, 4) for x in zero_rates] == [
        0.02, 0.0377, 0.0344, 0.0349, 0.0357]

    beta1 = 0.05
    curve3 = DiscountCurveNSS(start_date,
                              beta1, beta2, beta3, beta4, tau1, tau2)
    zero_rates = curve3.zero_rate(dates)
    assert [round(x, 4) for x in zero_rates] == [
        0.03, 0.0476, 0.0443, 0.0447, 0.0455]

    beta1 = 0.06
    curve4 = DiscountCurveNSS(start_date,
                              beta1, beta2, beta3, beta4, tau1, tau2)
    zero_rates = curve4.zero_rate(dates)
    assert [round(x, 4) for x in zero_rates] == [
        0.04, 0.0574, 0.0541, 0.0546, 0.0554]

    beta1 = 0.07
    curve5 = DiscountCurveNSS(start_date,
                              beta1, beta2, beta3, beta4, tau1, tau2)
    zero_rates = curve5.zero_rate(dates)
    assert [round(x, 4) for x in zero_rates] == [
        0.05, 0.0673, 0.064, 0.0644, 0.0652]
