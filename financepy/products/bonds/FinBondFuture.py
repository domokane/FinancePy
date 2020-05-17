# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2019

@author: Dominic O'Kane
"""

from ...finutils.FinGlobalVariables import gDaysInYear
from ...products.bonds.FinBond import FinBond

# TODO: Examine other exchange conventions.
# TODO: Delivery option model
##########################################################################


class FinBondFuture(object):
    ''' Class for managing futures contracts on government bonds that follows
    CME conventions and related analytics. '''

    def __init__(self,
                 tickerName,
                 firstDeliveryDate,
                 lastDeliveryDate,
                 contractSize,
                 coupon):

        self._tickerName = tickerName
        self._firstDeliveryDate = firstDeliveryDate
        self._lastDeliveryDate = lastDeliveryDate
        self._contractSize = contractSize
        self._coupon = coupon

##########################################################################

    def conversionFactor(self, bond):
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
        newBond = FinBond(newMat,
                          bond._coupon,
                          bond._frequencyType,
                          bond._accrualType,
                          face)

        p = newBond.cleanPriceFromYield(self._firstDeliveryDate,
                                        self._coupon)

        # Convention is to round the conversion factor to 4dp
        p = round(p, 4)
        return p

##########################################################################

    def principalInvoicePrice(self,
                              bond,
                              futuresPrice):
        ' The principal invoice price as defined by the CME.'''
        cf = self.conversionFactor(bond)
        pip = self._contractSize * (futuresPrice * cf) / 100.0
        pip = round(pip, 2)
        return pip

##########################################################################

    def totalInvoiceAmount(self,
                           settlementDate,
                           bond,
                           futuresPrice):
        ' The total invoice amount paid to take delivery of bond. '

        if bond._accruedInterest is None:
            bond.calculateFlowDates(settlementDate)

        accd = bond._accruedInterest

        pip = self.principalInvoicePrice(bond, futuresPrice)
        accrued = accd * self._contractSize / 100.0
        tia = pip + accrued
        tia = round(tia, 2)
        return tia

##########################################################################

    def cheapestToDeliver(self,
                          bonds,
                          bondCleanPrices,
                          futuresPrice):
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

##########################################################################

    def deliveryGainLoss(self,
                         bond,
                         bondCleanPrice,
                         futuresPrice):
        ''' Determination of what is received when the bond is delivered. '''
        receiveOnFuture = self.principalInvoicePrice(bond, futuresPrice)
        payForBond = self._contractSize * bondCleanPrice / 100.0
        net = receiveOnFuture - payForBond
        return net, payForBond, receiveOnFuture

##########################################################################
