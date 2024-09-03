###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.models.hw_tree import HWTree
from financepy.models.bk_tree import BKTree
from financepy.models.bdt_tree import BDTTree
from financepy.utils.global_types import OptionTypes
from financepy.products.bonds.bond_option import BondOption
from financepy.utils.global_vars import g_days_in_year
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.bonds.bond import Bond
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.utils.date import Date
import numpy as np


settle_dt = Date(1, 12, 2019)
issue_dt = Date(1, 12, 2018)
maturity_dt = settle_dt.add_tenor("10Y")
coupon = 0.05
freq_type = FrequencyTypes.SEMI_ANNUAL
dc_type = DayCountTypes.ACT_ACT_ICMA
bond = Bond(issue_dt, maturity_dt, coupon, freq_type, dc_type)

t_mat = (maturity_dt - settle_dt) / g_days_in_year
times = np.linspace(0, t_mat, 20)
dates = settle_dt.add_years(times)
dfs = np.exp(-0.05*times)
discount_curve = DiscountCurve(settle_dt, dates, dfs)

expiry_dt = settle_dt.add_tenor("18m")
face = 100.0

num_time_steps = 100


def test_european_call_bk():
    option_type = OptionTypes.EUROPEAN_CALL
    strike_price = 100

    bond_option = BondOption(
        bond, expiry_dt, strike_price, option_type)

    sigma = 0.20
    a = 0.1
    num_time_steps = 20
    model = BKTree(sigma, a, num_time_steps)

    v = bond_option.value(settle_dt, discount_curve, model)

    assert round(v, 4) == 1.7055


def test_american_call_bk():
    option_type = OptionTypes.AMERICAN_CALL
    strike_price = 100

    bond_option = BondOption(
        bond, expiry_dt, strike_price, option_type)

    sigma = 0.01
    a = 0.1
    model = BKTree(sigma, a)

    v = bond_option.value(settle_dt, discount_curve, model)

    assert round(v, 4) == 0.0069


def test_european_put_bk():
    option_type = OptionTypes.EUROPEAN_PUT
    strike_price = 100

    bond_option = BondOption(
        bond, expiry_dt, strike_price, option_type)

    sigma = 0.01
    a = 0.1
    model = BKTree(sigma, a)

    v = bond_option.value(settle_dt, discount_curve, model)

    assert round(v, 4) == 0.4060


def test_american_put_bk():
    option_type = OptionTypes.AMERICAN_PUT
    strike_price = 100

    bond_option = BondOption(
        bond, expiry_dt, strike_price, option_type)

    sigma = 0.02
    a = 0.1
    model = BKTree(sigma, a)

    v = bond_option.value(settle_dt, discount_curve, model)

    assert round(v, 4) == 0.5331


def test_european_call_bdt():
    option_type = OptionTypes.EUROPEAN_CALL
    strike_price = 100

    bond_option = BondOption(
        bond, expiry_dt, strike_price, option_type)

    sigma = 0.20
    model = BDTTree(sigma, num_time_steps)

    v = bond_option.value(settle_dt, discount_curve, model)

    assert round(v, 4) == 2.9156


def test_american_call_bdt():
    option_type = OptionTypes.AMERICAN_CALL
    strike_price = 100

    bond_option = BondOption(
        bond, expiry_dt, strike_price, option_type)

    sigma = 0.20
    model = BDTTree(sigma, num_time_steps)

    v = bond_option.value(settle_dt, discount_curve, model)

    assert round(v, 4) == 3.0939


def test_european_put_bdt():
    option_type = OptionTypes.EUROPEAN_PUT
    strike_price = 100

    bond_option = BondOption(
        bond, expiry_dt, strike_price, option_type)

    sigma = 0.01
    model = BDTTree(sigma, num_time_steps)

    v = bond_option.value(settle_dt, discount_curve, model)

    assert round(v, 4) == 0.4326


def test_american_put_bdt():
    option_type = OptionTypes.AMERICAN_PUT
    strike_price = 100

    bond_option = BondOption(
        bond, expiry_dt, strike_price, option_type)

    sigma = 0.02
    model = BDTTree(sigma, num_time_steps)

    v = bond_option.value(settle_dt, discount_curve, model)

    assert round(v, 4) == 0.6141


# Results different from TestFinBondOptionHWModel.py
# because t_mat != 10.0
def test_european_call_hw():
    option_type = OptionTypes.EUROPEAN_CALL
    strike_price = 100
    num_time_steps = 100

    bond_option = BondOption(
        bond, expiry_dt, strike_price, option_type)

    sigma = 0.01
    a = 0.1
    model = HWTree(sigma, a, num_time_steps)

    v = bond_option.value(settle_dt, discount_curve, model)

    assert round(v, 4) == 1.8809


def test_american_call_hw():
    option_type = OptionTypes.AMERICAN_CALL
    strike_price = 100

    bond_option = BondOption(
        bond, expiry_dt, strike_price, option_type)

    sigma = 0.01
    a = 0.1
    model = HWTree(sigma, a)

    v = bond_option.value(settle_dt, discount_curve, model)

    assert round(v, 4) == 2.0443


def test_european_put_hw():
    option_type = OptionTypes.EUROPEAN_PUT
    strike_price = 100

    bond_option = BondOption(
        bond, expiry_dt, strike_price, option_type)

    sigma = 0.01
    a = 0.1
    model = HWTree(sigma, a)

    v = bond_option.value(settle_dt, discount_curve, model)

    assert round(v, 4) == 2.2767


def test_american_put_hw():
    option_type = OptionTypes.AMERICAN_PUT
    strike_price = 100

    bond_option = BondOption(
        bond, expiry_dt, strike_price, option_type)

    sigma = 0.02
    a = 0.1
    model = HWTree(sigma, a)

    v = bond_option.value(settle_dt, discount_curve, model)

    assert round(v, 4) == 4.7948
