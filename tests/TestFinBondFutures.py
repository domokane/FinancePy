# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2019

@author: Dominic O'Kane
"""
from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.products.bonds.FinBondFuture import FinBondFuture
from financepy.products.bonds.FinBond import FinBond
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinDate import FinDate
import sys
sys.path.append("..//..")


testCases = FinTestCases(__file__, globalTestCaseMode)


def test_FinBondFuture():

        # Example taken from Martellini and Priaulet page 360
    freq = FinFrequencyTypes.SEMI_ANNUAL
    basis = FinDayCountTypes.ACT_ACT_ICMA

    bond1 = FinBond(FinDate(2011, 8, 15), 0.0500, freq, basis)
    bond2 = FinBond(FinDate(2011, 2, 15), 0.0500, freq, basis)
    bond3 = FinBond(FinDate(2010, 8, 15), 0.0575, freq, basis)
    bond4 = FinBond(FinDate(2010, 2, 15), 0.0650, freq, basis)
    bond5 = FinBond(FinDate(2009, 8, 15), 0.0600, freq, basis)
    bond6 = FinBond(FinDate(2009, 5, 15), 0.0550, freq, basis)
    bond7 = FinBond(FinDate(2008, 11, 15), 0.0475, freq, basis)

    bonds = []
    bonds.append(bond1)
    bonds.append(bond2)
    bonds.append(bond3)
    bonds.append(bond4)
    bonds.append(bond5)
    bonds.append(bond6)
    bonds.append(bond7)

    firstDeliveryDate = FinDate(2002, 3, 1)
    lastDeliveryDate = FinDate(2002, 3, 28)
    contractSize = 100000
    contractCoupon = 0.06

    bondFutureContract = FinBondFuture("TYH2",
                                       firstDeliveryDate,
                                       lastDeliveryDate,
                                       contractSize,
                                       contractCoupon)

    settlementDate = FinDate(2001, 12, 10)

    # Get the Conversion Factors
    testCases.header("Bond Maturity", "Coupon", "Conversion Factor")
    for bond in bonds:
        cf = bondFutureContract.conversionFactor(bond)
        testCases.print(bond._maturityDate, bond._coupon * 100, cf)

    # Example from
    # https://www.cmegroup.com/education/files/understanding-treasury-futures.pdf

    testCases.banner("EXAMPLE FROM CME")
    testCases.banner("================")
    settlementDate = FinDate(2017, 10, 10)

    bonds = []
    prices = []
    bond = FinBond(FinDate(2027, 8, 15), 0.0225, freq, basis)
    bonds.append(bond)
    prices.append(99 + 1 / 32)
    bond = FinBond(FinDate(2027, 5, 15), 0.02375, freq, basis)
    bonds.append(bond)
    prices.append(100 + 5 / 32 + 1 / 64)
    bond = FinBond(FinDate(2027, 2, 15), 0.0225, freq, basis)
    bonds.append(bond)
    prices.append(99 + 5 / 32 + 1 / 64)
    bond = FinBond(FinDate(2026, 11, 15), 0.02, freq, basis)
    bonds.append(bond)
    prices.append(97 + 7 / 32 + 1 / 64)
    bond = FinBond(FinDate(2026, 8, 15), 0.015, freq, basis)
    bonds.append(bond)
    prices.append(93 + 14 / 32)
    bond = FinBond(FinDate(2026, 5, 15), 0.01625, freq, basis)
    bonds.append(bond)
    prices.append(94 + 21 / 32 + 1 / 64)
    bond = FinBond(FinDate(2026, 2, 15), 0.01625, freq, basis)
    bonds.append(bond)
    prices.append(94 + 29 / 32)
    bond = FinBond(FinDate(2025, 11, 15), 0.0225, freq, basis)
    bonds.append(bond)
    prices.append(99 + 25 / 32)
    bond = FinBond(FinDate(2025, 8, 15), 0.02, freq, basis)
    bonds.append(bond)
    prices.append(98 + 3 / 32)
    bond = FinBond(FinDate(2025, 5, 15), 0.02125, freq, basis)
    bonds.append(bond)
    prices.append(99 + 5 / 32 + 1 / 64)
    bond = FinBond(FinDate(2025, 2, 15), 0.02, freq, basis)
    bonds.append(bond)
    prices.append(98 + 14 / 32 + 1 / 64)
    bond = FinBond(FinDate(2024, 11, 15), 0.0225, freq, basis)
    bonds.append(bond)
    prices.append(100 + 9 / 32 + 1 / 64)
    bond = FinBond(FinDate(2024, 8, 15), 0.02375, freq, basis)
    bonds.append(bond)
    prices.append(101 + 7 / 32 + 1 / 64)
    bond = FinBond(FinDate(2024, 8, 15), 0.01875, freq, basis)
    bonds.append(bond)
    # There may be an error in the document says 98-01+
    prices.append(98 + 1 / 32)

    testCases.header("BOND MATURITY", "YIELD")
    for bond, cleanPrice in zip(bonds, prices):
        yld = bond.yieldToMaturity(settlementDate, cleanPrice)
        testCases.print(str(bond._maturityDate), yld)

    firstDeliveryDate = FinDate(2017, 12, 1)
    lastDeliveryDate = FinDate(2017, 12, 28)

    contractSize = 100000
    contractCoupon = 0.06

    bondFutureContract = FinBondFuture("TYZ7",
                                       firstDeliveryDate,
                                       lastDeliveryDate,
                                       contractSize,
                                       contractCoupon)

    testCases.header("BOND MATURITY", "CF")
    for bond in bonds:
        cf = bondFutureContract.conversionFactor(bond)
        testCases.print(str(bond._maturityDate), cf)

    # Get the Invoice Prices
    futuresPrice = 125.265625

    testCases.header("BOND MATURITY", "PRINCIPAL INVOICE PRICE")
    for bond in bonds:
        pip = bondFutureContract.principalInvoicePrice(bond, futuresPrice)
        testCases.print(str(bond._maturityDate), pip)

    testCases.header("BOND MATURITY", "TOTAL INVOICE AMOUNT")
    for bond in bonds:
        tia = bondFutureContract.totalInvoiceAmount(
            settlementDate, bond, futuresPrice)
        testCases.print(str(bond._maturityDate), tia)

    ctd = bondFutureContract.cheapestToDeliver(bonds,
                                               prices,
                                               futuresPrice)

    testCases.header("CTD MATURITY", "CTD COUPON")
    testCases.print(str(ctd._maturityDate), ctd._coupon)

##########################################################################


test_FinBondFuture()
testCases.compareTestCases()
