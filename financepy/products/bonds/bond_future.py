###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from ...utils.global_vars import g_days_in_year
from ...products.bonds.bond import Bond
from ...utils.date import Date
from ...utils.helpers import label_to_string, check_argument_types


# TODO: Examine other exchange conventions.
# TODO: Delivery option model
###############################################################################


class BondFuture:
    """Class for managing futures contracts on government bonds that follows
    CME conventions and related analytics."""

    def __init__(
        self,
        ticker_name: str,
        first_delivery_dt: Date,
        last_delivery_dt: Date,
        contract_size: int,
        cpn: float,
    ):

        check_argument_types(self.__init__, locals())

        self.ticker_name = ticker_name
        self.first_delivery_dt = first_delivery_dt  # This is the IMM date
        self.last_delivery_dt = last_delivery_dt
        self.contract_size = contract_size
        self.cpn = cpn

    ###########################################################################

    def conversion_factor(self, bond: Bond):
        """Determine the conversion factor for a specific bond using CME
        convention. To do this we need to know the contract standard coupon and
        must round the bond maturity (starting its life on the first delivery
        date) to the nearest 3 month multiple and then calculate the bond clean
        price."""

        # See
        # https://www.cmegroup.com//trading//interest-rates//us-treasury-futures-conversion-factor-lookup-tables.html
        # for a reference.

        t_mat = (bond.maturity_dt - self.first_delivery_dt) / g_days_in_year
        roundedt_matInMonths = int(t_mat * 4.0) * 3
        new_mat = self.first_delivery_dt.add_months(roundedt_matInMonths)
        ex_div_days = 0

        issue_dt = Date(new_mat.d, new_mat.m, 2000)

        new_bond = Bond(
            issue_dt,
            new_mat,
            bond.cpn,
            bond.freq_type,
            bond.dc_type,
            ex_div_days,
        )

        p = new_bond.clean_price_from_ytm(self.first_delivery_dt, self.cpn)

        # Convention is to round the conversion factor to 4dp
        p = round(p, 4)
        return p

    ###########################################################################

    def principal_invoice_price(self, bond: Bond, futures_price: float):
        """The principal invoice price as defined by the CME."""
        cf = self.conversion_factor(bond)
        pip = self.contract_size * (futures_price * cf) / 100.0
        pip = round(pip, 2)
        return pip

    ###########################################################################

    def total_invoice_amount(
        self, settle_dt: Date, bond: Bond, futures_price: float
    ):
        "The total invoice amount paid to take delivery of bond."

        if bond.accrued_int is None:
            bond._calculate_cpn_dts(settle_dt)

        pip = self.principal_invoice_price(bond, futures_price)
        accrued = bond.accrued_int * self.contract_size / 100.0
        tia = pip + accrued
        tia = round(tia, 2)
        return tia

    ###########################################################################

    def cheapest_to_deliver(
        self, bonds: list, bond_clean_prices: list, futures_price: float
    ):
        """Determination of CTD as deliverable bond with the lowest cost to buy
        versus what is received when the bond is delivered."""
        ctd_bond = None
        ctd_net = -self.contract_size * 100
        for bondCleanPrice, bond in zip(bond_clean_prices, bonds):
            receive_on_future = self.principal_invoice_price(
                bond, futures_price
            )
            pay_for_bond = self.contract_size * bondCleanPrice / 100.0
            net = receive_on_future - pay_for_bond
            if net > ctd_net:
                ctd_bond = bond
                ctd_net = net

        return ctd_bond

    ###########################################################################

    def delivery_gain_loss(
        self, bond: Bond, bond_clean_price: float, futures_price: float
    ):
        """Determination of what is received when the bond is delivered."""
        receive_on_future = self.principal_invoice_price(bond, futures_price)
        pay_for_bond = self.contract_size * bond_clean_price / 100.0
        net = receive_on_future - pay_for_bond
        return net, pay_for_bond, receive_on_future

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("TICKER NAME", self.ticker_name)
        s += label_to_string("FIRST DELIVERY DATE", self.first_delivery_dt)
        s += label_to_string("LAST DELIVERY DATE", self.last_delivery_dt)
        s += label_to_string("CONTRACT SIZE", self.contract_size)
        s += label_to_string("COUPON", self.cpn)
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


###############################################################################
