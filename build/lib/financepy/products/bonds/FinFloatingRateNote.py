# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

from scipy import optimize

from ...finutils.FinDate import FinDate
from ...finutils.FinError import FinError
from ...finutils.FinFrequency import FinFrequency, FinFrequencyTypes
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes

#TODO: Need to complete and verify the risk sensitity calculations.

###############################################################################

def f(dm, *args):
    ''' Function used to do solve root search in discount margin calculation. '''
    self = args[0]
    settlementDate = args[1]
    nextCoupon = args[2]
    futureLibor = args[3]
    fullPrice = args[4]

    px = self.fullPriceFromDiscountMargin(settlementDate,
                                          nextCoupon,
                                          futureLibor,
                                          dm)

    objFn = px - fullPrice
    return objFn

################################################################################
################################################################################

class FinFloatingRateNote(object):

    ''' Class for managing floating rate notes that pay a floating index plus a 
    quoted margin.'''

    def __init__(self, 
                 maturityDate,
                 quotedMargin,
                 frequencyType, 
                 accrualType,
                 face = 100.0, 
                 redemption = 1.0):

        ''' Create FinFloatingRateNote object. '''

        if frequencyType not in FinFrequencyTypes:
           raise FinError("Invalid Frequency:" + str(frequencyType))
           return

        if accrualType not in FinDayCountTypes:
           raise FinError("Unknown Bond Accrued Convention type " + str(accrualType))

        self._maturityDate = maturityDate
        self._quotedMargin = quotedMargin
        self._frequencyType = frequencyType
        self._accrualType = accrualType
        self._flowDates = []
        self._frequency = FinFrequency(frequencyType)
        self._face = face
        self._redemption = redemption

        ''' I do not determine cashflow dates as I do not want to require
        users to supply the issue date and without that I do not know how 
        far to go back in the cashflow date schedule. '''

        self._settlementDate = FinDate(1900,1,1)

################################################################################

    def calculateFlowDates(self, settlementDate):
        ''' Determine the bond cashflow payment dates. '''

        # No need to generate flows if settlement date has not changed
        if settlementDate == self._settlementDate:
            return

        self._settlementDate = settlementDate
        self._flowDates = []

        nextDate = self._maturityDate
        months = int(12 / self._frequency)

        while nextDate._excelDate > settlementDate._excelDate:
            self._flowDates.append(nextDate)
            nextDate = nextDate.addMonths(-months)

        self._flowDates.reverse()
        return

###############################################################################

    def fullPriceFromDiscountMargin(self, 
                                    settlementDate,
                                    nextCoupon,
                                    futureLibor,
                                    dm ):
        ''' Calculate the full price of the bond from its discount margin and #
        making assumptions about the future Libor rates. '''

        self.calculateFlowDates(settlementDate)
        months = int(12 / self._frequency)

        ncd = self._flowDates[0]
        dc = FinDate.datediff(settlementDate,ncd)
        pcd = FinDate.addMonths(ncd,-months)
        dbc = FinDate.datediff(pcd,ncd)

        alpha = float(dc) / float(dbc)
        f = self._frequency
        q = self._quotedMargin
        y = futureLibor + dm
        c = futureLibor + q
        numFlows = len(self._flowDates)

        dayCount = FinDayCount(self._accrualType)
        yearFrac = dayCount.yearFrac(pcd,ncd)
        pv = (nextCoupon/f)/pow(1.0+y*yearFrac,alpha)

        for n in range(1,numFlows):
            pcd = self._flowDates[n-1]
            ncd = self._flowDates[n]
            yearFrac = dayCount.yearFrac(pcd,ncd)
            pv = pv + c*yearFrac/pow(1.0+y*yearFrac,alpha+n)

        pv += self._redemption/pow(1.0+y*yearFrac,alpha+n)
        return pv * self._face

###############################################################################

    def dollarDuration(self, 
                       settlementDate, 
                       nextCoupon,
                       futureLibor,
                       dm):
        ''' Calculate the risk or dP/dy of the bond by bumping. '''

        if dm > 10.0:
            raise FinError("Discount margin exceeds 100000bp")

        self.calculateFlowDates(settlementDate)
        dy = 0.0001

        p0 = self.fullPriceFromDiscountMargin(settlementDate,
                                              nextCoupon,
                                              futureLibor+dy,
                                              dm)

        p2 = self.fullPriceFromDiscountMargin(settlementDate,
                                              nextCoupon,
                                              futureLibor-dy,
                                              dm)
        durn = -(p2-p0)/dy/2.0
        durn *= self._frequency
        return durn

###############################################################################

    def macauleyDuration(self,
                         settlementDate,
                         nextCoupon,
                         futureLibor,
                         dm ):
        ''' Calculate the Macauley duration of the bond on a settlement date 
        given its yield to maturity. '''

        if dm > 10.0:
            raise FinError("Discount margin exceeds 100000bp")

        dd = self.dollarDuration(settlementDate,
                                 nextCoupon,
                                 futureLibor,
                                 dm )

        fp = self.fullPriceFromDiscountMargin(settlementDate,
                                              nextCoupon,
                                              futureLibor,
                                              dm )

        lastLiborReset = nextCoupon - self._quotedMargin
        md = dd * (1.0 + (lastLiborReset + dm)/self._frequency) / fp
        return md

###############################################################################

    def modifiedDuration(self,
                         settlementDate,
                         nextCoupon,
                         futureLibor,
                         dm ):

        ''' Calculate the modified duration of the bondon a settlement date 
        given its yield to maturity. '''

        if dm > 10.0:
            raise FinError("Discount margin exceeds 100000bp")

        dd = self.dollarDuration(settlementDate,
                                 nextCoupon,
                                 futureLibor,
                                 dm)

        fp = self.fullPriceFromDiscountMargin(settlementDate,
                                              nextCoupon,
                                              futureLibor,
                                              dm)
        md = dd/fp
        return md

###############################################################################

    def convexityFromDiscountMargin(self, 
                                    settlementDate,
                                    nextCoupon,
                                    futureLibor,
                                    dm ):
        ''' Calculate the bond convexity from the yield to maturity. '''

        if dm > 10.0:
            raise FinError("Discount margin exceeds 100000bp")

        self.calculateFlowDates(settlementDate)
        dy = 0.0001

        p0 = self.fullPriceFromDiscountMargin(settlementDate,
                                              nextCoupon,
                                              futureLibor,
                                              dm - dy)
              
        p1 = self.fullPriceFromDiscountMargin(settlementDate,
                                              nextCoupon,
                                              futureLibor,
                                              dm )

        p2 = self.fullPriceFromDiscountMargin(settlementDate,
                                              nextCoupon,
                                              futureLibor,
                                              dm + dy)

        conv = ((p2+p0)-2.0*p1)/dy/dy/p1/self._face
        return conv

###############################################################################

    def cleanPriceFromDiscountMargin(self, 
                                     settlementDate,
                                     nextCoupon,
                                     futureLibor,
                                     dm ):
        ''' Calculate the bond clean price from the yield. '''

        if dm > 10.0:
            raise FinError("Discount margin exceeds 100000bp")

        fullPrice = self.fullPriceFromDiscountMargin(settlementDate,
                                                     nextCoupon, 
                                                     futureLibor, 
                                                     dm)

        accrued = self.accruedInterest(settlementDate,nextCoupon)
        cleanPrice = fullPrice - accrued
        return cleanPrice

###############################################################################

    def fullPriceFromDiscountCurve(self, 
                                   settlementDate,
                                   indexCurve,
                                   discountCurve):

        ''' Calculate the bond price using some discount curve to present-value 
        the bond's cashflows. THIS IS NOT COMPLETE. '''
        self.calculateFlowDates(settlementDate)

        pv = 0.0

        for dt in self._flowDates:
            df = discountCurve.df(dt)
            flow = self._coupon / self._frequency
            pv = pv + flow * df

        pv = pv + df 
        return pv * self._face

################################################################################

    def discountMargin(self, 
                       settlementDate,
                       nextCoupon, 
                       futureLibor,
                       cleanPrice):

        ''' Calculate the bond's yield to maturity by solving the price 
        yield relationship using a one-dimensional root solver. '''
        self.calculateFlowDates(settlementDate)

        accrued = self.accruedInterest(settlementDate,nextCoupon)
        fullPrice = cleanPrice + accrued
        print("Clean",cleanPrice)
        print("Accrued",accrued)
        print("Full",fullPrice)
        
        argtuple = (self,settlementDate,nextCoupon,futureLibor,fullPrice)

        dm = optimize.newton(f,
                             x0=0.01, # initial value of 10%
                             fprime=None,
                             args=argtuple,
                             tol=1e-8,
                             maxiter=50,
                             fprime2=None)
        return dm

################################################################################

    def accruedDays(self, settlementDate):
        ''' Calculate number of days from previous coupon date to settlement.'''
        self.calculateFlowDates(settlementDate)

        if len(self._flowDates) <= 2:
            raise FinError("Accrued interest - not enough flow dates.")

        return settlementDate - self.pcd(settlementDate)

################################################################################

    def pcd(self, settlementDate):
        ''' Determine the previous coupon date before the settlement date. '''
        self.calculateFlowDates(settlementDate)

        if len(self._flowDates) <= 2:
            raise FinError("Accrued interest - not enough flow dates.")

        months = int(12/self._frequency)
        ncd = self._flowDates[0]
        pcd = FinDate.addMonths(ncd,-months)
        return pcd

################################################################################

    def accruedInterest(self, 
                        settlementDate, 
                        nextCoupon):
        ''' Calculate the amount of coupon that has accrued between the 
        previous coupon date and the settlement date. '''

        self.calculateFlowDates(settlementDate)
        pcd = self.pcd(settlementDate)
        dayCount = FinDayCount(self._accrualType)
        accruedPeriod = dayCount.yearFrac(pcd,settlementDate)
        accrued = accruedPeriod * nextCoupon
        return accrued * self._face

################################################################################
