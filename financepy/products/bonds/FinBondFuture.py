###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from ...finutils.FinGlobalVariables import gDaysInYear
from ...products.bonds.FinBond import FinBond
from ...finutils.FinDate import FinDate
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes


# TODO: Examine other exchange conventions.
# TODO: Delivery option model
###############################################################################


class FinBondFuture(object):
    ''' Class for managing futures contracts on government bonds that follows
    CME conventions and related analytics. '''

    def __init__(self,
                 tickerName: str,
                 firstDeliveryDate: FinDate,
                 lastDeliveryDate: FinDate,
                 contractSize: int,
                 coupon: float):

        checkArgumentTypes(self.__init__, locals())

        self._tickerName = tickerName
        self._firstDeliveryDate = firstDeliveryDate  # This is the IMM date
        self._lastDeliveryDate = lastDeliveryDate
        self._contractSize = contractSize
        self._coupon = coupon

###############################################################################

    def conversionFactor(self,
                         bond: FinBond):
        ''' Determine the conversion factor for a specific bond using CME
        convention. To do this we need to know the contract standard coupon and
        must round the bond maturity (starting its life on the first delivery
        date) to the nearest 3 month multiple and then calculate the bond clean
        price. '''

        # See
        # https://www.cmegroup.com//trading//interest-rates//us-treasury-futures-conversion-factor-lookup-tables.html
        # for a reference.

        tmat = (bond._maturityDate - self._firstDeliveryDate) / gDaysInYear
        roundedTmatInMonths = int(tmat * 4.0) * 3
        newMat = self._firstDeliveryDate.addMonths(roundedTmatInMonths)
        face = 1.0

        issueDate = FinDate(newMat._d, newMat._m, 2000)

        newBond = FinBond(issueDate,
                          newMat,
                          bond._coupon,
                          bond._freqType,
                          bond._accrualType,
                          face)

        p = newBond.cleanPriceFromYTM(self._firstDeliveryDate,
                                      self._coupon)

        # Convention is to round the conversion factor to 4dp
        p = round(p, 4)
        return p

###############################################################################

    def principalInvoicePrice(self,
                              bond: FinBond,
                              futuresPrice: float):
        ' The principal invoice price as defined by the CME.'''
        cf = self.conversionFactor(bond)
        pip = self._contractSize * (futuresPrice * cf) / 100.0
        pip = round(pip, 2)
        return pip

###############################################################################

    def totalInvoiceAmount(self,
                           settlementDate: FinDate,
                           bond: FinBond,
                           futuresPrice: float):
        ' The total invoice amount paid to take delivery of bond. '

        if bond._accruedInterest is None:
            bond.calculateFlowDates(settlementDate)

        accd = bond._accruedInterest

        pip = self.principalInvoicePrice(bond, futuresPrice)
        accrued = accd * self._contractSize / 100.0
        tia = pip + accrued
        tia = round(tia, 2)
        return tia

###############################################################################

    def cheapestToDeliver(self,
                          bonds: list,
                          bondCleanPrices: list,
                          futuresPrice: float):
        ''' Determination of CTD as deliverable bond with lowest cost to buy
        versus what is received when the bond is delivered. '''
        ctdBond = None
        ctdNet = -self._contractSize * 100
        for bondCleanPrice, bond in zip(bondCleanPrices, bonds):
            receiveOnFuture = self.principalInvoicePrice(bond, futuresPrice)
            payForBond = self._contractSize * bondCleanPrice / 100.0
            net = receiveOnFuture - payForBond
            if net > ctdNet:
                ctdBond = bond
                ctdNet = net

        return ctdBond

###############################################################################

    def deliveryGainLoss(self,
                         bond: FinBond,
                         bondCleanPrice: float,
                         futuresPrice: float):
        ''' Determination of what is received when the bond is delivered. '''
        receiveOnFuture = self.principalInvoicePrice(bond, futuresPrice)
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
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
