###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.bonds.bond import Bond
from financepy.products.bonds.bond_future import BondFuture

freq = FrequencyTypes.SEMI_ANNUAL
basis = DayCountTypes.ACT_ACT_ICMA
issue_date = Date(15, 2, 2004)


def test_bond_future_1():
    bond = Bond(issue_date, Date(15, 8, 2011), 0.0500, freq, basis)

    assert bond._maturity_date == Date(15, 8, 2011)
    assert bond._coupon * 100 == 5.0

    first_delivery_date = Date(1, 3, 2002)
    last_delivery_date = Date(28, 3, 2002)
    contract_size = 100000
    contractCoupon = 0.06
    bondFutureContract = BondFuture("TYH2",
                                    first_delivery_date,
                                    last_delivery_date,
                                    contract_size,
                                    contractCoupon)

    cf = bondFutureContract.conversion_factor(bond)

    assert round(cf, 4) == 92.9688


def test_bond_future_2():
    bond = Bond(issue_date, Date(15, 8, 2027), 0.0225, freq, basis)
    assert bond._maturity_date == Date(15, 8, 2027)

    settlement_date = Date(10, 10, 2017)
    price = 99 + 1 / 32

    yld = bond.yield_to_maturity(settlement_date, price)

    assert round(yld, 4) == 0.0236

    first_delivery_date = Date(1, 12, 2017)
    last_delivery_date = Date(28, 12, 2017)

    contract_size = 100000
    contractCoupon = 0.06
    bondFutureContract = BondFuture("TYZ7",
                                    first_delivery_date,
                                    last_delivery_date,
                                    contract_size,
                                    contractCoupon)

    cf = bondFutureContract.conversion_factor(bond)

    assert round(cf, 4) == 73.1429

    futures_price = 125.265625

    pip = bondFutureContract.principal_invoice_price(bond, futures_price)

    assert round(pip, 4) == 9162291.0800

    tia = bondFutureContract.total_invoice_amount(
        settlement_date, bond, futures_price)

    assert round(tia, 4) == 9162294.5

###############################################################################

def test_future_bond_ctd():

    first_delivery_date = Date(1, 12, 2017)
    last_delivery_date = Date(28, 12, 2017)

    contract_size = 100000
    contractCoupon = 0.06

    bondFutureContract = BondFuture("TYZ7",
                                    first_delivery_date,
                                    last_delivery_date,
                                    contract_size,
                                    contractCoupon)

    bonds = []
    prices = []
    bond = Bond(issue_date, Date(15, 8, 2027), 0.0225, freq, basis)
    bonds.append(bond)
    prices.append(99 + 1 / 32)
    bond = Bond(issue_date, Date(15, 5, 2027), 0.02375, freq, basis)
    bonds.append(bond)
    prices.append(100 + 5 / 32 + 1 / 64)
    bond = Bond(issue_date, Date(15, 2, 2027), 0.0225, freq, basis)
    bonds.append(bond)
    prices.append(99 + 5 / 32 + 1 / 64)
    bond = Bond(issue_date, Date(15, 11, 2026), 0.02, freq, basis)
    bonds.append(bond)
    prices.append(97 + 7 / 32 + 1 / 64)
    bond = Bond(issue_date, Date(15, 8, 2026), 0.015, freq, basis)
    bonds.append(bond)
    prices.append(93 + 14 / 32)
    bond = Bond(issue_date, Date(15, 5, 2026), 0.01625, freq, basis)
    bonds.append(bond)
    prices.append(94 + 21 / 32 + 1 / 64)
    bond = Bond(issue_date, Date(15, 2, 2026), 0.01625, freq, basis)
    bonds.append(bond)
    prices.append(94 + 29 / 32)
    bond = Bond(issue_date, Date(15, 11, 2025), 0.0225, freq, basis)
    bonds.append(bond)
    prices.append(99 + 25 / 32)
    bond = Bond(issue_date, Date(15, 8, 2025), 0.02, freq, basis)
    bonds.append(bond)
    prices.append(98 + 3 / 32)
    bond = Bond(issue_date, Date(15, 5, 2025), 0.02125, freq, basis)
    bonds.append(bond)
    prices.append(99 + 5 / 32 + 1 / 64)
    bond = Bond(issue_date, Date(15, 2, 2025), 0.02, freq, basis)
    bonds.append(bond)
    prices.append(98 + 14 / 32 + 1 / 64)
    bond = Bond(issue_date, Date(15, 11, 2024), 0.0225, freq, basis)
    bonds.append(bond)
    prices.append(100 + 9 / 32 + 1 / 64)
    bond = Bond(issue_date, Date(15, 8, 2024), 0.02375, freq, basis)
    bonds.append(bond)
    prices.append(101 + 7 / 32 + 1 / 64)
    bond = Bond(issue_date, Date(15, 8, 2024), 0.01875, freq, basis)
    bonds.append(bond)
    # There may be an error in the document says 98-01+
    prices.append(98 + 1 / 32)

    futures_price = 125.265625

    ctd = bondFutureContract.cheapest_to_deliver(bonds,
                                                 prices,
                                                 futures_price)

    assert round(ctd._coupon, 4) == 0.0238
