###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from ...utils.global_vars import gDaysInYear
from ...products.bonds.bond import Bond
from ...utils.date import Date
from ...utils.helpers import labelToString, check_argument_types


# TODO: Examine other exchange conventions.
# TODO: Delivery option model
###############################################################################


class BondFuture(object):
    """ Class for managing futures contracts on government bonds that follows
    CME conventions and related analytics. """

    def __init__(self,
                 tickerName: str,
                 firstDeliveryDate: Date,
                 lastDeliveryDate: Date,
                 contractSize: int,
                 coupon: float):

        check_argument_types(self.__init__, locals())

        self._tickerName = tickerName
        self._firstDeliveryDate = firstDeliveryDate  # This is the IMM date
        self._lastDeliveryDate = lastDeliveryDate
        self._contractSize = contractSize
        self._coupon = coupon

###############################################################################

    def conversionFactor(self,
                         bond: Bond):
        """ Determine the conversion factor for a specific bond using CME
        convention. To do this we need to know the contract standard coupon and
        must round the bond maturity (starting its life on the first delivery
        date) to the nearest 3 month multiple and then calculate the bond clean
        price. """

        # See
        # https://www.cmegroup.com//trading//interest-rates//us-treasury-futures-conversion-factor-lookup-tables.html
        # for a reference.

        tmat = (bond._maturity_date - self._firstDeliveryDate) / gDaysInYear
        roundedTmatInMonths = int(tmat * 4.0) * 3
        newMat = self._firstDeliveryDate.addMonths(roundedTmatInMonths)
        face = 1.0

        issue_date = Date(newMat._d, newMat._m, 2000)

        newBond = Bond(issue_date,
                          newMat,
                          bond._coupon,
                          bond._freq_type,
                          bond._accrual_type,
                          face)

        p = newBond.clean_price_from_ytm(self._firstDeliveryDate,
                                      self._coupon)

        # Convention is to round the conversion factor to 4dp
        p = round(p, 4)
        return p

###############################################################################

    def principalInvoicePrice(self,
                              bond: Bond,
                              futures_price: float):
        """ The principal invoice price as defined by the CME."""
        cf = self.conversionFactor(bond)
        pip = self._contractSize * (futures_price * cf) / 100.0
        pip = round(pip, 2)
        return pip

###############################################################################

    def totalInvoiceAmount(self,
                           settlement_date: Date,
                           bond: Bond,
                           futures_price: float):
        ' The total invoice amount paid to take delivery of bond. '

        if bond._accruedInterest is None:
            bond.calculate_flow_dates(settlement_date)

        accrued_interest= bond._accruedInterest

        pip = self.principalInvoicePrice(bond, futures_price)
        accrued = accrued_interest* self._contractSize / 100.0
        tia = pip + accrued
        tia = round(tia, 2)
        return tia

###############################################################################

    def cheapestToDeliver(self,
                          bonds: list,
                          bondCleanPrices: list,
                          futures_price: float):
        """ Determination of CTD as deliverable bond with lowest cost to buy
        versus what is received when the bond is delivered. """
        ctdBond = None
        ctdNet = -self._contractSize * 100
        for bondCleanPrice, bond in zip(bondCleanPrices, bonds):
            receiveOnFuture = self.principalInvoicePrice(bond, futures_price)
            payForBond = self._contractSize * bondCleanPrice / 100.0
            net = receiveOnFuture - payForBond
            if net > ctdNet:
                ctdBond = bond
                ctdNet = net

        return ctdBond

###############################################################################

    def deliveryGainLoss(self,
                         bond: Bond,
                         bondCleanPrice: float,
                         futures_price: float):
        """ Determination of what is received when the bond is delivered. """
        receiveOnFuture = self.principalInvoicePrice(bond, futures_price)
        payForBond = self._contractSize * bondCleanPrice / 100.0
        net = receiveOnFuture - payForBond
        return net, payForBond, receiveOnFuture

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("TICKER NAME", self._tickerName)
        s += labelToString("FIRST DELIVERY DATE", self._firstDeliveryDate)
        s += labelToString("LAST DELIVERY DATE", self._lastDeliveryDate)
        s += labelToString("CONTRACT SIZE", self._contractSize)
        s += labelToString("COUPON", self._coupon)
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
