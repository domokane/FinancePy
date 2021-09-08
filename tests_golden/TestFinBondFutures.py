###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.bonds.bond import Bond
from financepy.products.bonds.bond_future import BondFuture
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)


def test_BondFuture():

    # Example taken from Martellini and Priaulet page 360
    freq = FrequencyTypes.SEMI_ANNUAL
    basis = DayCountTypes.ACT_ACT_ICMA
    issue_date = Date(15, 2, 2004)

    bond1 = Bond(issue_date, Date(15, 8, 2011), 0.0500, freq, basis)
    bond2 = Bond(issue_date, Date(15, 2, 2011), 0.0500, freq, basis)
    bond3 = Bond(issue_date, Date(15, 8, 2010), 0.0575, freq, basis)
    bond4 = Bond(issue_date, Date(15, 2, 2010), 0.0650, freq, basis)
    bond5 = Bond(issue_date, Date(15, 8, 2009), 0.0600, freq, basis)
    bond6 = Bond(issue_date, Date(15, 5, 2009), 0.0550, freq, basis)
    bond7 = Bond(issue_date, Date(15, 11, 2008), 0.0475, freq, basis)

    bonds = []
    bonds.append(bond1)
    bonds.append(bond2)
    bonds.append(bond3)
    bonds.append(bond4)
    bonds.append(bond5)
    bonds.append(bond6)
    bonds.append(bond7)

    first_delivery_date = Date(1, 3, 2002)
    last_delivery_date = Date(28, 3, 2002)
    contract_size = 100000
    contractCoupon = 0.06

    bondFutureContract = BondFuture("TYH2",
                                    first_delivery_date,
                                    last_delivery_date,
                                    contract_size,
                                    contractCoupon)

    settlement_date = Date(10, 12, 2001)

    # Get the Conversion Factors
    testCases.header("Bond Maturity", "Coupon", "Conversion Factor")
    for bond in bonds:
        cf = bondFutureContract.conversion_factor(bond)
        testCases.print(bond._maturity_date, bond._coupon * 100, cf)

    # Example from
    # https://www.cmegroup.com/education/files/understanding-treasury-futures.pdf

    testCases.banner("EXAMPLE FROM CME")
    testCases.banner("================")
    settlement_date = Date(10, 10, 2017)

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

    testCases.header("BOND MATURITY", "YIELD")
    for bond, clean_price in zip(bonds, prices):
        yld = bond.yield_to_maturity(settlement_date, clean_price)
        testCases.print(str(bond._maturity_date), yld)

    first_delivery_date = Date(1, 12, 2017)
    last_delivery_date = Date(28, 12, 2017)

    contract_size = 100000
    contractCoupon = 0.06

    bondFutureContract = BondFuture("TYZ7",
                                    first_delivery_date,
                                    last_delivery_date,
                                    contract_size,
                                    contractCoupon)

    testCases.header("BOND MATURITY", "CF")
    for bond in bonds:
        cf = bondFutureContract.conversion_factor(bond)
        testCases.print(str(bond._maturity_date), cf)

    # Get the Invoice Prices
    futures_price = 125.265625

    testCases.header("BOND MATURITY", "PRINCIPAL INVOICE PRICE")
    for bond in bonds:
        pip = bondFutureContract.principal_invoice_price(bond, futures_price)
        testCases.print(str(bond._maturity_date), pip)

    testCases.header("BOND MATURITY", "TOTAL INVOICE AMOUNT")
    for bond in bonds:
        tia = bondFutureContract.total_invoice_amount(
            settlement_date, bond, futures_price)
        testCases.print(str(bond._maturity_date), tia)

    ctd = bondFutureContract.cheapest_to_deliver(bonds,
                                                 prices,
                                                 futures_price)

    testCases.header("CTD MATURITY", "CTD COUPON")
    testCases.print(str(ctd._maturity_date), ctd._coupon)

##########################################################################


test_BondFuture()
testCases.compareTestCases()
