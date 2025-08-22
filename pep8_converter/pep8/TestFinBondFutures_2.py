# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import sys

from FinTestCases import FinTestCases, global_test_case_mode
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.bonds.bond import Bond, YTMCalcType
from financepy.products.bonds.bond_future import BondFuture

import pandas as pd

sys.path.append("..")

test_cases = FinTestCases(__file__, global_test_case_mode)




########################################################################################


def test_BondFutures():

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
    contract_cpn = 0.06

    bfut = BondFuture(
        "TYH2",
        first_delivery_dt,
        last_delivery_dt,
        contract_size,
        contract_cpn,
    )

    settle_dt = Date(10, 12, 2001)

    # Get the Conversion Factors
    test_cases.header("Bond Maturity", "Coupon", "Conversion Factor")
    for bond in bonds:
        cf = bfut.conversion_factor(bond)
        test_cases.print(bond.maturity_dt, bond.cpn * 100, cf)


########################################################################################


########################################################################################


def test_BondFutures_CME_two_bond_examples():

    # Example from
    # https://www.cmegroup.com/education/files/understanding-treasury-futures.pdf

    freq = FrequencyTypes.SEMI_ANNUAL
    basis = DayCountTypes.ACT_ACT_ICMA
    issue_dt = Date(15, 2, 2004)

    test_cases.banner("EXAMPLE FROM CME")
    test_cases.banner("================")
    settle_dt = Date(10, 10, 2017)

    first_del_dt = Date(1, 12, 2017)
    last_del_dt = Date(29, 12, 2017)

    fut_size = 100000
    fut_coupon = 0.06

    bfut = BondFuture("TYZ7", first_del_dt, last_del_dt, fut_size, fut_coupon)

    # Get the Invoice Prices
    fut_price = 125.265625

    # REPRODUCING EXAMPLES IN REPORT TABLE PAGE 6
    bond1 = Bond(issue_dt, Date(15, 8, 2024), 0.02375, freq, basis)
    bond2 = Bond(issue_dt, Date(15, 8, 2024), 0.01875, freq, basis)

    print("Futures Price       %12.6f %12.6f" % (fut_price, fut_price))

    cf1 = bfut.conversion_factor(bond1)
    cf2 = bfut.conversion_factor(bond2)
    print("x CF                %12.4f %12.4f" % (cf1, cf2))
    print("x 1000              %12.2f %12.2f" % (1000, 1000))

    pip1 = bfut.principal_invoice(bond1, fut_price)
    pip2 = bfut.principal_invoice(bond2, fut_price)
    print("Principal invoice   %12.2f %12.2f " % (pip1, pip2))

    price1 = 101 + 7 / 32 + 1 / 64
    price2 = 98 + 1 / 32 + 0 / 64
    cash1 = price1 * 1000
    cash2 = price2 * 1000
    print("Cash Price          %12.2f %12.2f" % (-cash1, -cash2))

    tia1 = bfut.total_invoice_amount(settle_dt, bond1, fut_price)
    tia2 = bfut.total_invoice_amount(settle_dt, bond2, fut_price)
    print("Total Invoice price %12.2f %12.2f" % (tia1, tia2))

    delgainloss1 = bfut.delivery_gain_loss(bond1, price1, fut_price)
    delgainloss2 = bfut.delivery_gain_loss(bond2, price2, fut_price)
    print("Delivery Gain/Loss  %12.2f %12.2f" % (delgainloss1, delgainloss2))

    basis1 = bfut.basis(bond1, fut_price, price1)
    basis2 = bfut.basis(bond2, fut_price, price2)
    print("Basis (32nds)       %12.2f %12.2f" % (basis1 * 32, basis2 * 32))

    basis1 = bfut.basis_net_of_carry(bond1, settle_dt, fut_price, price1)
    basis2 = bfut.basis_net_of_carry(bond2, settle_dt, fut_price, price2)
    print("Basis (32nds)       %12.2f %12.2f" % (basis1, basis2))

    ###########################################################################


########################################################################################


def test_BondFutures_CME_table():

    freq = FrequencyTypes.SEMI_ANNUAL
    basis = DayCountTypes.ACT_ACT_ICMA
    issue_dt = Date(15, 2, 2004)
    yield_convention = YTMCalcType.US_TREASURY

    test_cases.banner("TABLE 3 EXAMPLE FROM CME")
    test_cases.banner("========================")
    settle_dt = Date(10, 10, 2017)

    bonds = []
    prices = []
    clean_prices = []

    if 1 == 1:
        # bond 1
        bond = Bond(issue_dt, Date(15, 8, 2027), 0.0225, freq, basis)
        bonds.append(bond)
        prices.append(99 + 1 / 32)
        clean_prices.append(99.0391)

        # bond 2
        bond = Bond(issue_dt, Date(15, 5, 2027), 0.02375, freq, basis)
        bonds.append(bond)
        prices.append(100 + 5 / 32 + 1 / 64)
        clean_prices.append(100.168)

        # bond 3
        bond = Bond(issue_dt, Date(15, 2, 2027), 0.0225, freq, basis)
        bonds.append(bond)
        prices.append(99 + 5 / 32 + 1 / 64)
        clean_prices.append(99.1641)

        # bond 4
        bond = Bond(issue_dt, Date(15, 11, 2026), 0.02, freq, basis)
        bonds.append(bond)
        prices.append(97 + 7 / 32 + 1 / 64)
        clean_prices.append(97.2305)

        # bond 5
        bond = Bond(issue_dt, Date(15, 8, 2026), 0.015, freq, basis)
        bonds.append(bond)
        prices.append(93 + 14 / 32)
        clean_prices.append(93.4414)

        # bond 6
        bond = Bond(issue_dt, Date(15, 5, 2026), 0.01625, freq, basis)
        bonds.append(bond)
        prices.append(94 + 21 / 32 + 1 / 64)
        clean_prices.append(94.6641)

        # bond 7
        bond = Bond(issue_dt, Date(15, 2, 2026), 0.01625, freq, basis)
        bonds.append(bond)
        prices.append(94 + 29 / 32)
        clean_prices.append(94.9063)

        # bond 8
        bond = Bond(issue_dt, Date(15, 11, 2025), 0.0225, freq, basis)
        bonds.append(bond)
        prices.append(99 + 25 / 32)
        clean_prices.append(99.7813)

        # bond 9
        bond = Bond(issue_dt, Date(15, 8, 2025), 0.02, freq, basis)
        bonds.append(bond)
        prices.append(98 + 3 / 32)
        clean_prices.append(98.0938)

        # bond 10
        bond = Bond(issue_dt, Date(15, 5, 2025), 0.02125, freq, basis)
        bonds.append(bond)
        prices.append(99 + 5 / 32 + 1 / 64)
        clean_prices.append(99.1719)

        # bond 11
        bond = Bond(issue_dt, Date(15, 2, 2025), 0.02, freq, basis)
        bonds.append(bond)
        prices.append(98 + 14 / 32 + 1 / 64)
        clean_prices.append(98.4531)

        # bond 12
        bond = Bond(issue_dt, Date(15, 11, 2024), 0.0225, freq, basis)
        bonds.append(bond)
        prices.append(100 + 9 / 32 + 1 / 64)
        clean_prices.append(100.3008)

        # bond 13
        bond = Bond(issue_dt, Date(30, 9, 2024), 0.02125, freq, basis)
        bonds.append(bond)
        prices.append(99.6016)
        clean_prices.append(99.6016)

        # bond 14
        bond = Bond(issue_dt, Date(31, 8, 2024), 0.01875, freq, basis)
        bonds.append(bond)
        prices.append(98 + 1 / 32)
        clean_prices.append(98.0508)

        # bond 15
        bond = Bond(issue_dt, Date(15, 8, 2024), 0.02375, freq, basis)
        bonds.append(bond)
        prices.append(101 + 7 / 32 + 1 / 64)
        clean_prices.append(101.2266)

        # bond 16
        bond = Bond(issue_dt, Date(31, 7, 2024), 0.02125, freq, basis)
        bonds.append(bond)
        prices.append(99.6758)
        clean_prices.append(99.6758)

        # bond 17
        bond = Bond(issue_dt, Date(30, 6, 2024), 0.02, freq, basis)
        bonds.append(bond)
        prices.append(98.9336)
        clean_prices.append(98.9336)

    # bonds.reverse()
    # prices.reverse()
    # new_prices.reverse()

    test_cases.header("BOND MATURITY", "COUPON", "PRICE")
    for bond, clean_price in zip(bonds, clean_prices):
        test_cases.print(str(bond.maturity_dt), str(bond.cpn), clean_price)

    test_cases.header("BOND MATURITY", "COUPON", "YIELD")
    for bond, clean_price in zip(bonds, clean_prices):
        yld = bond.yield_to_maturity(settle_dt, clean_price)
        test_cases.print(str(bond.maturity_dt), str(bond.cpn), yld)

    first_delivery_dt = Date(1, 12, 2017)
    last_delivery_dt = Date(29, 12, 2017)

    contract_size = 100000
    contract_cpn = 0.06

    bfut = BondFuture(
        "TYZ7",
        first_delivery_dt,
        last_delivery_dt,
        contract_size,
        contract_cpn,
    )

    test_cases.header("BOND MATURITY", "COUPON", "CF")
    for bond in bonds:
        cf = bfut.conversion_factor(bond)
        test_cases.print(str(bond.maturity_dt), str(bond.cpn), cf)

    # Get the Invoice Prices
    futures_price = 125.265625

    test_cases.header("BOND MATURITY", "PRINCIPAL INVOICE PRICE")
    for bond in bonds:
        pip = bfut.principal_invoice(bond, futures_price)
        test_cases.print(str(bond.maturity_dt), pip)

    test_cases.header("BOND MATURITY", "TOTAL INVOICE AMOUNT")
    for bond in bonds:
        tia = bfut.total_invoice_amount(settle_dt, bond, futures_price)
        test_cases.print(str(bond.maturity_dt), tia)

    test_cases.header("BOND MATURITY", "IMPLIED REPO RATE")
    for bond, clean_price in zip(bonds, clean_prices):
        repo_rate = bfut.implied_repo_rate(bond, settle_dt, clean_price, futures_price)
        test_cases.print(str(bond.maturity_dt), repo_rate)

    ctd = bfut.ctd(bonds, prices, futures_price)
    test_cases.header("CTD MATURITY", "CTD COUPON")
    test_cases.print(str(ctd.maturity_dt), ctd.cpn)

    # Prepare list of dicts to collect results
    results = []
    repo_rate = 0.015

    if 1 == 0:
        for bond, clean_price in zip(bonds, clean_prices):

            cf = bfut.conversion_factor(bond)

            yld = bond.yield_to_maturity(settle_dt, clean_price, yield_convention)

            del_years = bfut.delivery_years(bond)

            pip = bfut.principal_invoice(bond, futures_price)

            tia = bfut.total_invoice_amount(settle_dt, bond, futures_price)

            gross_basis = bfut.gross_basis(bond, clean_price, futures_price)

            net_basis = bfut.net_basis(
                bond, settle_dt, clean_price, futures_price, repo_rate
            )

            irr = bfut.implied_repo_rate(bond, settle_dt, clean_price, futures_price)

            fut_dv01 = 0.0

            # Append dictionary for each bond
            results.append(
                {
                    "Coupon": bond.cpn * 100,
                    "Maturity": bond.maturity_dt,
                    "Clean": clean_price,
                    "del_yrs": del_years,
                    "TCF": cf,
                    "PIP": pip,
                    "TIA": tia,
                    "GROSS_BASIS": gross_basis,
                    "NET_BASIS": net_basis,
                    "FUT_DV01": fut_dv01,
                    "Yield (%)": yld * 100,
                    "IRR (%)": irr * 100,
                }
            )

        # Create data_frame
        df = pd.data_frame(results)
        df = df.sort_values(by="IRR (%)", ascending=False).reset_index(drop=True)

        formatted_df = df.copy()
        formatted_df["Coupon"] = formatted_df["Coupon"].map("{:.3f}".format)
        formatted_df["del_yrs"] = formatted_df["del_yrs"].map("{:.2f}".format)
        formatted_df["TCF"] = formatted_df["TCF"].map("{:.4f}".format)
        formatted_df["PIP"] = formatted_df["PIP"].map("{:.2f}".format)
        formatted_df["TIA"] = formatted_df["TIA"].map("{:.2f}".format)
        formatted_df["GROSS_BASIS"] = formatted_df["GROSS_BASIS"].map("{:.3f}".format)
        formatted_df["NET_BASIS"] = formatted_df["NET_BASIS"].map("{:.3f}".format)
        formatted_df["FUT_DV01"] = formatted_df["FUT_DV01"].map("{:.3f}".format)
        formatted_df["IRR (%)"] = formatted_df["IRR (%)"].map("{:.3f}".format)
        formatted_df["Yield (%)"] = formatted_df["Yield (%)"].map("{:.3f}".format)

        print(formatted_df.to_string(index=False))

    ##########################################################################


########################################################################################


def test_BondFutures_BBG_table():

    # https://i.sstatic.net/gg6eg.png

    freq = FrequencyTypes.SEMI_ANNUAL
    basis = DayCountTypes.ACT_ACT_ICMA
    issue_dt = Date(15, 2, 2004)
    yield_convention = YTMCalcType.US_TREASURY

    test_cases.banner("TABLE 3 EXAMPLE FROM CME")
    test_cases.banner("========================")
    settle_dt = Date(10, 10, 2017)

    bonds = []
    prices = []
    new_prices = []

    bond = Bond(issue_dt, Date(31, 7, 2024), 0.02125, freq, basis)
    bonds.append(bond)
    new_prices.append(98.5703)

    bond = Bond(issue_dt, Date(30, 6, 2024), 0.020, freq, basis)
    bonds.append(bond)
    new_prices.append(97.8516)

    bond = Bond(issue_dt, Date(31, 8, 2024), 0.01875, freq, basis)
    bonds.append(bond)
    new_prices.append(97.0469)

    bond = Bond(issue_dt, Date(30, 11, 2024), 0.02125, freq, basis)
    bonds.append(bond)
    new_prices.append(98.4297)

    bond = Bond(issue_dt, Date(31, 10, 2024), 0.0225, freq, basis)
    bonds.append(bond)
    new_prices.append(99.2578)

    bond = Bond(issue_dt, Date(15, 11, 2024), 0.0225, freq, basis)
    bonds.append(bond)
    new_prices.append(99.2188)

    bond = Bond(issue_dt, Date(30, 9, 2024), 0.02125, freq, basis)
    bonds.append(bond)
    new_prices.append(98.4844)

    bond = Bond(issue_dt, Date(15, 8, 2024), 0.02375, freq, basis)
    bonds.append(bond)
    new_prices.append(101.2266)

    # bonds.reverse()
    # prices.reverse()
    # new_prices.reverse()

    test_cases.header("BOND MATURITY", "COUPON", "PRICE")
    for bond, clean_price in zip(bonds, new_prices):
        test_cases.print(str(bond.maturity_dt), str(bond.cpn), clean_price)

    test_cases.header("BOND MATURITY", "COUPON", "YIELD")
    for bond, clean_price in zip(bonds, new_prices):
        yld = bond.yield_to_maturity(settle_dt, clean_price, yield_convention)
        test_cases.print(str(bond.maturity_dt), str(bond.cpn), yld)

    first_delivery_dt = Date(1, 12, 2017)
    last_delivery_dt = Date(29, 12, 2017)

    contract_size = 100000
    contract_cpn = 0.06

    bfut = BondFuture(
        "TYZ7",
        first_delivery_dt,
        last_delivery_dt,
        contract_size,
        contract_cpn,
    )

    test_cases.header("BOND MATURITY", "COUPON", "CF")
    for bond in bonds:
        cf = bfut.conversion_factor(bond)
        test_cases.print(str(bond.maturity_dt), str(bond.cpn), cf)

    # Get the Invoice Prices
    futures_price = 125.265625

    test_cases.header("BOND MATURITY", "PRINCIPAL INVOICE PRICE")
    for bond in bonds:
        pip = bfut.principal_invoice(bond, futures_price)
        test_cases.print(str(bond.maturity_dt), pip)

    test_cases.header("BOND MATURITY", "TOTAL INVOICE AMOUNT")
    for bond in bonds:
        tia = bfut.total_invoice_amount(settle_dt, bond, futures_price)
        test_cases.print(str(bond.maturity_dt), tia)

    ctd = bfut.ctd(bonds, new_prices, futures_price)
    test_cases.header("CTD MATURITY", "CTD COUPON")
    test_cases.print(str(ctd.maturity_dt), ctd.cpn)

    # Prepare list of dicts to collect results
    results = []
    repo_rate = 0.01499

    if 1 == 0:
        for bond, clean_price in zip(bonds, new_prices):

            mid_price = clean_price

            cf = bfut.conversion_factor(bond)

            yld = bond.yield_to_maturity(settle_dt, mid_price, yield_convention)

            del_years = bfut.delivery_years(bond)

            pip = bfut.principal_invoice(bond, futures_price)

            tia = bfut.total_invoice_amount(settle_dt, bond, futures_price)

            gross_basis = bfut.gross_basis(bond, mid_price, futures_price)

            irr = bfut.implied_repo_rate(bond, settle_dt, mid_price, futures_price)

            net_basis = bfut.net_basis(
                bond, settle_dt, clean_price, futures_price, repo_rate
            )

            # Append dictionary for each bond
            results.append(
                {
                    "Coupon": bond.cpn * 100,
                    "Maturity": bond.maturity_dt,
                    "Mid": mid_price,
                    "del_yrs": del_years,
                    "TCF": cf,
                    "PIP": pip,
                    "TIA": tia,
                    "GROSS_BASIS": gross_basis,
                    "NET_BASIS": net_basis,
                    "Yield (%)": yld * 100,
                    "IRR (%)": irr * 100,
                }
            )

        # Create data_frame
        df = pd.data_frame(results)
        #    df = df.sort_values(by="IRR (%)", ascending=False).reset_index(drop=True)

        formatted_df = df.copy()
        formatted_df["Coupon"] = formatted_df["Coupon"].map("{:.3f}".format)
        formatted_df["del_yrs"] = formatted_df["del_yrs"].map("{:.2f}".format)
        formatted_df["TCF"] = formatted_df["TCF"].map("{:.4f}".format)
        formatted_df["PIP"] = formatted_df["PIP"].map("{:.2f}".format)
        formatted_df["TIA"] = formatted_df["TIA"].map("{:.2f}".format)
        formatted_df["GROSS_BASIS"] = formatted_df["GROSS_BASIS"].map("{:.3f}".format)
        formatted_df["NET_BASIS"] = formatted_df["NET_BASIS"].map("{:.3f}".format)
        formatted_df["IRR (%)"] = formatted_df["IRR (%)"].map("{:.3f}".format)
        formatted_df["Yield (%)"] = formatted_df["Yield (%)"].map("{:.6f}".format)

        print(formatted_df.to_string(index=False))

    ##########################################################################


# test_BondFutures()
# test_BondFutures_CME_two_bond_examples()
test_BondFutures_CME_table()
test_BondFutures_BBG_table()
test_cases.compare_test_cases()

########################################################################################

