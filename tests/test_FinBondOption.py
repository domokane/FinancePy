###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.models.hw_tree import HWTree
from financepy.models.bk_tree import BKTree
from financepy.models.bdt_tree import BDTTree
from financepy.utils.global_types import OptionTypes
from financepy.products.bonds.bond_option import BondOption
from financepy.utils.global_vars import gDaysInYear
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.bonds.bond import Bond
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.utils.date import Date
import numpy as np


settlement_date = Date(1, 12, 2019)
issue_date = Date(1, 12, 2018)
maturity_date = settlement_date.add_tenor("10Y")
coupon = 0.05
freq_type = FrequencyTypes.SEMI_ANNUAL
accrual_type = DayCountTypes.ACT_ACT_ICMA
bond = Bond(issue_date, maturity_date, coupon, freq_type, accrual_type)

tmat = (maturity_date - settlement_date) / gDaysInYear
times = np.linspace(0, tmat, 20)
dates = settlement_date.add_years(times)
dfs = np.exp(-0.05*times)
discount_curve = DiscountCurve(settlement_date, dates, dfs)

expiry_date = settlement_date.add_tenor("18m")
face = 100.0

num_time_steps = 100


def test_european_call_bk():
    option_type = OptionTypes.EUROPEAN_CALL
    strike_price = 100

    bond_option = BondOption(
        bond, expiry_date, strike_price, option_type)

    sigma = 0.20
    a = 0.1
    num_time_steps = 20
    model = BKTree(sigma, a, num_time_steps)

    v = bond_option.value(settlement_date, discount_curve, model)

    assert round(v, 4) == 1.7055


def test_american_call_bk():
    option_type = OptionTypes.AMERICAN_CALL
    strike_price = 100

    bond_option = BondOption(
        bond, expiry_date, strike_price, option_type)

    sigma = 0.01
    a = 0.1
    model = BKTree(sigma, a)

    v = bond_option.value(settlement_date, discount_curve, model)

    assert round(v, 4) == 0.0069


def test_european_put_bk():
    option_type = OptionTypes.EUROPEAN_PUT
    strike_price = 100

    bond_option = BondOption(
        bond, expiry_date, strike_price, option_type)

    sigma = 0.01
    a = 0.1
    model = BKTree(sigma, a)

    v = bond_option.value(settlement_date, discount_curve, model)

    assert round(v, 4) == 0.4060


def test_american_put_bk():
    option_type = OptionTypes.AMERICAN_PUT
    strike_price = 100

    bond_option = BondOption(
        bond, expiry_date, strike_price, option_type)

    sigma = 0.02
    a = 0.1
    model = BKTree(sigma, a)

    v = bond_option.value(settlement_date, discount_curve, model)

    assert round(v, 4) == 0.5331


def test_european_call_bdt():
    option_type = OptionTypes.EUROPEAN_CALL
    strike_price = 100

    bond_option = BondOption(
        bond, expiry_date, strike_price, option_type)

    sigma = 0.20
    model = BDTTree(sigma, num_time_steps)

    v = bond_option.value(settlement_date, discount_curve, model)

    assert round(v, 4) == 2.9156


def test_american_call_bdt():
    option_type = OptionTypes.AMERICAN_CALL
    strike_price = 100

    bond_option = BondOption(
        bond, expiry_date, strike_price, option_type)

    sigma = 0.20
    model = BDTTree(sigma, num_time_steps)

    v = bond_option.value(settlement_date, discount_curve, model)

    assert round(v, 4) == 3.0939


def test_european_put_bdt():
    option_type = OptionTypes.EUROPEAN_PUT
    strike_price = 100

    bond_option = BondOption(
        bond, expiry_date, strike_price, option_type)

    sigma = 0.01
    model = BDTTree(sigma, num_time_steps)

    v = bond_option.value(settlement_date, discount_curve, model)

    assert round(v, 4) == 0.4326


def test_american_put_bdt():
    option_type = OptionTypes.AMERICAN_PUT
    strike_price = 100

    bond_option = BondOption(
        bond, expiry_date, strike_price, option_type)

    sigma = 0.02
    model = BDTTree(sigma, num_time_steps)

    v = bond_option.value(settlement_date, discount_curve, model)

    assert round(v, 4) == 0.6141


# Results different from TestFinBondOptionHWModel.py
# because tmat != 10.0
def test_european_call_hw():
    option_type = OptionTypes.EUROPEAN_CALL
    strike_price = 100
    num_time_steps = 100

    bond_option = BondOption(
        bond, expiry_date, strike_price, option_type)

    sigma = 0.01
    a = 0.1
    model = HWTree(sigma, a, num_time_steps)

    v = bond_option.value(settlement_date, discount_curve, model)

    assert round(v, 4) == 1.8809


def test_american_call_hw():
    option_type = OptionTypes.AMERICAN_CALL
    strike_price = 100

    bond_option = BondOption(
        bond, expiry_date, strike_price, option_type)

    sigma = 0.01
    a = 0.1
    model = HWTree(sigma, a)

    v = bond_option.value(settlement_date, discount_curve, model)

    assert round(v, 4) == 2.0443


def test_european_put_hw():
    option_type = OptionTypes.EUROPEAN_PUT
    strike_price = 100

    bond_option = BondOption(
        bond, expiry_date, strike_price, option_type)

    sigma = 0.01
    a = 0.1
    model = HWTree(sigma, a)

    v = bond_option.value(settlement_date, discount_curve, model)

    assert round(v, 4) == 2.2767


def test_american_put_hw():
    option_type = OptionTypes.AMERICAN_PUT
    strike_price = 100

    bond_option = BondOption(
        bond, expiry_date, strike_price, option_type)

    sigma = 0.02
    a = 0.1
    model = HWTree(sigma, a)

    v = bond_option.value(settlement_date, discount_curve, model)

    assert round(v, 4) == 4.7948
