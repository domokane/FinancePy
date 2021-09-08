###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from ...utils.global_vars import gDaysInYear
from ...products.bonds.bond import Bond
from ...utils.date import Date
from ...utils.helpers import label_to_string, check_argument_types


# TODO: Examine other exchange conventions.
# TODO: Delivery option model
###############################################################################


class BondFuture:
    """ Class for managing futures contracts on government bonds that follows
    CME conventions and related analytics. """

    def __init__(self,
                 ticker_name: str,
                 first_delivery_date: Date,
                 last_delivery_date: Date,
                 contract_size: int,
                 coupon: float):

        check_argument_types(self.__init__, locals())

        self._ticker_name = ticker_name
        self._first_delivery_date = first_delivery_date  # This is the IMM date
        self._last_delivery_date = last_delivery_date
        self._contract_size = contract_size
        self._coupon = coupon

###############################################################################

    def conversion_factor(self,
                          bond: Bond):
        """ Determine the conversion factor for a specific bond using CME
        convention. To do this we need to know the contract standard coupon and
        must round the bond maturity (starting its life on the first delivery
        date) to the nearest 3 month multiple and then calculate the bond clean
        price. """

        # See
        # https://www.cmegroup.com//trading//interest-rates//us-treasury-futures-conversion-factor-lookup-tables.html
        # for a reference.

        tmat = (bond._maturity_date - self._first_delivery_date) / gDaysInYear
        roundedTmatInMonths = int(tmat * 4.0) * 3
        newMat = self._first_delivery_date.add_months(roundedTmatInMonths)
        face = 1.0

        issue_date = Date(newMat._d, newMat._m, 2000)

        newBond = Bond(issue_date,
                       newMat,
                       bond._coupon,
                       bond._freq_type,
                       bond._accrual_type,
                       face)

        p = newBond.clean_price_from_ytm(self._first_delivery_date,
                                         self._coupon)

        # Convention is to round the conversion factor to 4dp
        p = round(p, 4)
        return p

###############################################################################

    def principal_invoice_price(self,
                                bond: Bond,
                                futures_price: float):
        """ The principal invoice price as defined by the CME."""
        cf = self.conversion_factor(bond)
        pip = self._contract_size * (futures_price * cf) / 100.0
        pip = round(pip, 2)
        return pip

###############################################################################

    def total_invoice_amount(self,
                             settlement_date: Date,
                             bond: Bond,
                             futures_price: float):
        ' The total invoice amount paid to take delivery of bond. '

        if bond._accrued_interest is None:
            bond.calculate_flow_dates(settlement_date)

        accrued_interest = bond._accrued_interest

        pip = self.principal_invoice_price(bond, futures_price)
        accrued = accrued_interest * self._contract_size / 100.0
        tia = pip + accrued
        tia = round(tia, 2)
        return tia

###############################################################################

    def cheapest_to_deliver(self,
                            bonds: list,
                            bond_clean_prices: list,
                            futures_price: float):
        """ Determination of CTD as deliverable bond with lowest cost to buy
        versus what is received when the bond is delivered. """
        ctdBond = None
        ctdNet = -self._contract_size * 100
        for bondCleanPrice, bond in zip(bond_clean_prices, bonds):
            receiveOnFuture = self.principal_invoice_price(bond, futures_price)
            payForBond = self._contract_size * bondCleanPrice / 100.0
            net = receiveOnFuture - payForBond
            if net > ctdNet:
                ctdBond = bond
                ctdNet = net

        return ctdBond

###############################################################################

    def delivery_gain_loss(self,
                           bond: Bond,
                           bond_clean_price: float,
                           futures_price: float):
        """ Determination of what is received when the bond is delivered. """
        receiveOnFuture = self.principal_invoice_price(bond, futures_price)
        payForBond = self._contract_size * bond_clean_price / 100.0
        net = receiveOnFuture - payForBond
        return net, payForBond, receiveOnFuture

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("TICKER NAME", self._ticker_name)
        s += label_to_string("FIRST DELIVERY DATE", self._first_delivery_date)
        s += label_to_string("LAST DELIVERY DATE", self._last_delivery_date)
        s += label_to_string("CONTRACT SIZE", self._contract_size)
        s += label_to_string("COUPON", self._coupon)
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
