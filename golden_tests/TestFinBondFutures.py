###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.bonds.bond import Bond
from financepy.products.bonds.bond_future import BondFuture


test_cases = FinTestCases(__file__, globalTestCaseMode)


def test_BondFuture():

    # Example taken from Martellini and Priaulet page 360
    freq = FrequencyTypes.SEMI_ANNUAL
    basis = DayCountTypes.ACT_ACT_ICMA
    issue_dt = Date(15, 2, 2004)

    bond1 = Bond(issue_dt, Date(15, 8, 2011), 0.0500, freq, basis)
    bond2 = Bond(issue_dt, Date(15, 2, 2011), 0.0500, freq, basis)
    bond3 = Bond(issue_dt, Date(15, 8, 2010), 0.0575, freq, basis)
    bond4 = Bond(issue_dt, Date(15, 2, 2010), 0.0650, freq, basis)
    bond5 = Bond(issue_dt, Date(15, 8, 2009), 0.0600, freq, basis)
    bond6 = Bond(issue_dt, Date(15, 5, 2009), 0.0550, freq, basis)
    bond7 = Bond(issue_dt, Date(15, 11, 2008), 0.0475, freq, basis)

    bonds = []
    bonds.append(bond1)
    bonds.append(bond2)
    bonds.append(bond3)
    bonds.append(bond4)
    bonds.append(bond5)
    bonds.append(bond6)
    bonds.append(bond7)

    first_delivery_dt = Date(1, 3, 2002)
    last_delivery_dt = Date(28, 3, 2002)
    contract_size = 100000
    contractCoupon = 0.06

    bondFutureContract = BondFuture("TYH2",
                                    first_delivery_dt,
                                    last_delivery_dt,
                                    contract_size,
                                    contractCoupon)

    settle_dt = Date(10, 12, 2001)

    # Get the Conversion Factors
    test_cases.header("Bond Maturity", "Coupon", "Conversion Factor")
    for bond in bonds:
        cf = bondFutureContract.conversion_factor(bond)
        test_cases.print(bond.maturity_dt, bond.cpn * 100, cf)

    # Example from
    # https://www.cmegroup.com/education/files/understanding-treasury-futures.pdf

    test_cases.banner("EXAMPLE FROM CME")
    test_cases.banner("================")
    settle_dt = Date(10, 10, 2017)

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

    bonds.reverse()
    prices.reverse()

    test_cases.header("BOND MATURITY", "COUPON", "PRICE")
    for bond, clean_price in zip(bonds, prices):
        test_cases.print(str(bond.maturity_dt), str(bond.cpn), clean_price)

    test_cases.header("BOND MATURITY", "COUPON", "YIELD")
    for bond, clean_price in zip(bonds, prices):
        yld = bond.yield_to_maturity(settle_dt, clean_price)
        test_cases.print(str(bond.maturity_dt), str(bond.cpn), yld)

    first_delivery_dt = Date(1, 12, 2017)
    last_delivery_dt = Date(28, 12, 2017)

    contract_size = 100000
    contractCoupon = 0.06

    bondFutureContract = BondFuture("TYZ7",
                                    first_delivery_dt,
                                    last_delivery_dt,
                                    contract_size,
                                    contractCoupon)

    test_cases.header("BOND MATURITY", "COUPON", "CF")
    for bond in bonds:
        cf = bondFutureContract.conversion_factor(bond)
        test_cases.print(str(bond.maturity_dt), str(bond.cpn), cf)

    # Get the Invoice Prices
    futures_price = 125.265625

    test_cases.header("BOND MATURITY", "PRINCIPAL INVOICE PRICE")
    for bond in bonds:
        pip = bondFutureContract.principal_invoice_price(bond, futures_price)
        test_cases.print(str(bond.maturity_dt), pip)

    test_cases.header("BOND MATURITY", "TOTAL INVOICE AMOUNT")
    for bond in bonds:
        tia = bondFutureContract.total_invoice_amount(
            settle_dt, bond, futures_price)
        test_cases.print(str(bond.maturity_dt), tia)

    ctd = bondFutureContract.cheapest_to_deliver(bonds,
                                                 prices,
                                                 futures_price)

    test_cases.header("CTD MATURITY", "CTD COUPON")
    test_cases.print(str(ctd.maturity_dt), ctd.cpn)

##########################################################################


test_BondFuture()
test_cases.compareTestCases()
