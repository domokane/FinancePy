# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import pytest

from financepy.products.bonds.bond import Bond
from financepy.products.bonds.bond_future import BondFuture
from financepy.utils.calendar import BusDayAdjustTypes, CalendarTypes
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes

freq = FrequencyTypes.SEMI_ANNUAL
basis = DayCountTypes.ACT_ACT_ICMA
issue_dt = Date(15, 2, 2004)

########################################################################################


def test_bond_future_1():

    bond = Bond(issue_dt, Date(15, 8, 2011), 0.0500, freq, basis)

    assert bond.maturity_dt == Date(15, 8, 2011)
    assert bond.cpn * 100 == 5.0

    first_delivery_date = Date(1, 3, 2002)
    last_delivery_date = Date(28, 3, 2002)
    contract_size = 100000
    contract_cpn = 0.06
    bond_fut_contract = BondFuture(
        "TYH2",
        first_delivery_date,
        last_delivery_date,
        contract_size,
        contract_cpn,
    )

    cf = bond_fut_contract.conversion_factor(bond)

    assert round(cf, 4) == 0.9297


########################################################################################


def test_bond_future_2():

    bond = Bond(issue_dt, Date(15, 8, 2027), 0.0225, freq, basis)
    assert bond.maturity_dt == Date(15, 8, 2027)

    settle_dt = Date(10, 10, 2017)
    price = 99 + 1 / 32

    yld = bond.yield_to_maturity(settle_dt, price)

    assert round(yld, 4) == 0.0236

    first_delivery_date = Date(1, 12, 2017)
    last_delivery_date = Date(28, 12, 2017)

    contract_size = 100000
    contract_cpn = 0.06
    bond_fut_contract = BondFuture(
        "TYZ7",
        first_delivery_date,
        last_delivery_date,
        contract_size,
        contract_cpn,
    )

    cf = bond_fut_contract.conversion_factor(bond)

    assert round(cf, 4) == 0.7314

    futures_price = 125.265625

    pip = bond_fut_contract.principal_invoice(bond, futures_price)

    assert round(pip, 4) == 91619.2781

    tia = bond_fut_contract.total_invoice_amount(settle_dt, bond, futures_price)

    assert round(tia, 4) == 91961.6694


########################################################################################


def test_future_bond_ctd():

    first_delivery_date = Date(1, 12, 2017)
    last_delivery_date = Date(28, 12, 2017)

    contract_size = 100000
    contract_cpn = 0.06

    bond_fut_contract = BondFuture(
        "TYZ7",
        first_delivery_date,
        last_delivery_date,
        contract_size,
        contract_cpn,
    )

    bonds = []
    prices = []
    bond = Bond(issue_dt, Date(15, 8, 2027), 0.0225, freq, basis)
    bonds.append(bond)
    prices.append(99 + 1 / 32)
    bond = Bond(issue_dt, Date(15, 5, 2027), 0.02375, freq, basis)
    bonds.append(bond)
    prices.append(100 + 5 / 32 + 1 / 64)
    bond = Bond(issue_dt, Date(15, 2, 2027), 0.0225, freq, basis)
    bonds.append(bond)
    prices.append(99 + 5 / 32 + 1 / 64)
    bond = Bond(issue_dt, Date(15, 11, 2026), 0.02, freq, basis)
    bonds.append(bond)
    prices.append(97 + 7 / 32 + 1 / 64)
    bond = Bond(issue_dt, Date(15, 8, 2026), 0.015, freq, basis)
    bonds.append(bond)
    prices.append(93 + 14 / 32)
    bond = Bond(issue_dt, Date(15, 5, 2026), 0.01625, freq, basis)
    bonds.append(bond)
    prices.append(94 + 21 / 32 + 1 / 64)
    bond = Bond(issue_dt, Date(15, 2, 2026), 0.01625, freq, basis)
    bonds.append(bond)
    prices.append(94 + 29 / 32)
    bond = Bond(issue_dt, Date(15, 11, 2025), 0.0225, freq, basis)
    bonds.append(bond)
    prices.append(99 + 25 / 32)
    bond = Bond(issue_dt, Date(15, 8, 2025), 0.02, freq, basis)
    bonds.append(bond)
    prices.append(98 + 3 / 32)
    bond = Bond(issue_dt, Date(15, 5, 2025), 0.02125, freq, basis)
    bonds.append(bond)
    prices.append(99 + 5 / 32 + 1 / 64)
    bond = Bond(issue_dt, Date(15, 2, 2025), 0.02, freq, basis)
    bonds.append(bond)
    prices.append(98 + 14 / 32 + 1 / 64)
    bond = Bond(issue_dt, Date(15, 11, 2024), 0.0225, freq, basis)
    bonds.append(bond)
    prices.append(100 + 9 / 32 + 1 / 64)
    bond = Bond(issue_dt, Date(15, 8, 2024), 0.02375, freq, basis)
    bonds.append(bond)
    prices.append(101 + 7 / 32 + 1 / 64)
    bond = Bond(issue_dt, Date(15, 8, 2024), 0.01875, freq, basis)
    bonds.append(bond)
    # There may be an error in the document says 98-01+
    prices.append(98 + 1 / 32)

    futures_price = 125.265625

    rows = []

    for bond, clean_price in zip(bonds, prices):
        cf = bond_fut_contract.conversion_factor(bond)
        gross_basis = bond_fut_contract.gross_basis(
            bond,
            clean_price,
            futures_price,
        )

        rows.append(
            (
                bond.maturity_dt,
                bond.cpn,
                clean_price,
                cf,
                gross_basis,
            )
        )

    rows = sorted(rows, key=lambda x: x[-1])

    for row in rows:
        print(row)

    ctd = bond_fut_contract.ctd(bonds, prices, futures_price)

    assert (
        ctd
        is bonds[
            min(
                range(len(bonds)),
                key=lambda i: bond_fut_contract.gross_basis(
                    bonds[i],
                    prices[i],
                    futures_price,
                ),
            )
        ]
    )


@pytest.fixture
def repo_bond():
    return Bond(
        Date(15, 2, 2020),
        Date(15, 2, 2032),
        0.06,
        FrequencyTypes.SEMI_ANNUAL,
        DayCountTypes.ACT_ACT_ICMA,
        cal_type=CalendarTypes.NONE,
        bd_type=BusDayAdjustTypes.NONE,
    )


@pytest.mark.parametrize("repo_rate", [0.01, 0.05])
@pytest.mark.parametrize(
    "settle_dt,delivery_dt",
    [
        (Date(1, 5, 2026), Date(1, 8, 2026)),  # no intervening coupon
        (Date(1, 8, 2026), Date(15, 8, 2026)),  # coupon on delivery
        (Date(15, 8, 2026), Date(1, 12, 2026)),  # coupon on settlement
        (Date(1, 8, 2026), Date(1, 12, 2026)),  # early coupon
        (Date(1, 3, 2026), Date(1, 9, 2026)),  # late coupon
        (Date(1, 1, 2026), Date(1, 3, 2027)),  # multiple coupons
    ],
)
def test_implied_repo_inverts_forward_price(
    repo_bond, settle_dt, delivery_dt, repo_rate
):
    future = BondFuture("TEST", delivery_dt, delivery_dt, 100_000, 0.06)
    clean_price = 99.0
    forward_price = repo_bond.forward_price(
        settle_dt, delivery_dt, clean_price, repo_rate
    )
    futures_price = forward_price / future.conversion_factor(repo_bond)

    implied = future.implied_repo_rate(repo_bond, settle_dt, clean_price, futures_price)

    assert implied == pytest.approx(repo_rate, abs=1e-12, rel=0.0)


@pytest.mark.parametrize("repo_rate", [-0.005, 0.0, 0.03])
def test_implied_repo_balances_coupon_reinvestment_cashflows(repo_bond, repo_rate):
    settle_dt = Date(1, 8, 2026)
    coupon_dt = Date(15, 8, 2026)
    delivery_dt = Date(1, 12, 2026)
    future = BondFuture("TEST", delivery_dt, delivery_dt, 100_000, 0.06)
    clean_price = 99.0
    full_price = clean_price + repo_bond.accrued_interest(settle_dt)
    holding_years = (delivery_dt - settle_dt) / 360.0
    reinvestment_years = (delivery_dt - coupon_dt) / 360.0
    coupon = 3.0  # six-percent semiannual coupon per 100 face

    # Financing the purchase equals the delivery invoice plus the reinvested coupon.
    invoice = full_price * (1.0 + repo_rate * holding_years) - coupon * (
        1.0 + repo_rate * reinvestment_years
    )
    futures_price = (
        invoice - repo_bond.accrued_interest(delivery_dt)
    ) / future.conversion_factor(repo_bond)
    implied = future.implied_repo_rate(repo_bond, settle_dt, clean_price, futures_price)

    assert implied == pytest.approx(repo_rate, abs=1e-12, rel=0.0)
