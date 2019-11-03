# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""
# TODO
# INCORPORATE OTHER YIELD MEASURES AND SMALL TWEAKS
# OTHER SPREADS SUCH AS I-SPREAD

from ...finutils.FinDate import FinDate
from ...finutils.FinError import FinError
from ...finutils.FinFrequency import FinFrequency, FinFrequencyTypes
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinSchedule import FinSchedule
from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinBusDayConventionTypes, FinDateGenRuleTypes

from math import pow

from scipy import optimize

from enum import Enum

# NOT SURE I NEED ALL OF THESE


class FinBondAccruedTypes(Enum):
    THIRTY_360 = 1
    THIRTY_360_BOND = 2
    ACT_360 = 3
    ACT_365 = 4
    ACT_ACT = 5

###############################################################################


def f(y, *args):
    ''' Function used to do solve root search in price to yield calculation. '''
    self = args[0]
    settlementDate = args[1]
    price = args[2]
    px = self.fullPriceFromYield(settlementDate, y)
    objFn = px - price
    return objFn

###############################################################################


def g(oas, *args):
    ''' Function used to do solve root search in price to yield calculation. '''
    self = args[0]
    settlementDate = args[1]
    price = args[2]
    discountCurve = args[3]
    px = self.fullPriceFromOAS(settlementDate, discountCurve, oas)
    objFn = px - price
    return objFn

###############################################################################


class FinBond(object):

    ''' Class for managing fixed coupon bonds and performing related analytics. '''

    def __init__(self,
                 maturityDate,
                 coupon,
                 frequencyType,
                 accrualType,
                 face=100.0,
                 redemption=1.0):
        ''' Create FinBond object by providing Maturity Date, Frequency,
        coupon and the accrual convention type. '''

        if frequencyType not in FinFrequencyTypes:
            raise FinError("Invalid Frequency:" + str(frequencyType))
            return

        if accrualType not in FinBondAccruedTypes:
            raise FinError(
                "Unknown Bond Accrued Convention type " +
                str(accrualType))

        self._maturityDate = maturityDate
        self._coupon = coupon
        self._frequencyType = frequencyType
        self._accrualType = accrualType
        self._flowDates = []
        self._frequency = FinFrequency(frequencyType)
        self._face = face
        self._redemption = redemption

        ''' I do not determine cashflow dates as I do not want to require
        users to supply the issue date and without that I do not know how
        far to go back in the cashflow date schedule. '''

        self._settlementDate = FinDate(1900, 1, 1)

##########################################################################

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

    def fullPriceFromYield(self, settlementDate, ytm):
        ''' Calculate the full price of the bond from its yield to maturity. '''

        self.calculateFlowDates(settlementDate)

        alpha = self.initialPeriodFraction(settlementDate)
        f = self._frequencyType.value
        c = self._coupon

        numFlows = len(self._flowDates)

        pv = 0.0
        for n in range(0, numFlows):
            pv = pv + (c / f) / pow(1.0 + ytm / f, alpha + n)

        pv += self._redemption / pow(1.0 + ytm / f, alpha + n)
        return pv * self._face

###############################################################################

    def dollarDuration(self, settlementDate, ytm):
        ''' Calculate the risk or dP/dy of the bond by bumping. '''

        self.calculateFlowDates(settlementDate)
        dy = 0.0001
        p0 = self.fullPriceFromYield(settlementDate, ytm - dy)
        p2 = self.fullPriceFromYield(settlementDate, ytm + dy)
        durn = -(p2 - p0) / dy / 2.0
        return durn

###############################################################################

    def macauleyDuration(self, settlementDate, ytm):
        ''' Calculate the Macauley duration of the bond on a settlement date
        given its yield to maturity. '''

        dd = self.dollarDuration(settlementDate, ytm)
        fp = self.fullPriceFromYield(settlementDate, ytm)
        md = dd * (1.0 + ytm / self._frequency) / fp
        return md

###############################################################################

    def modifiedDuration(self, settlementDate, ytm):
        ''' Calculate the modified duration of the bondon a settlement date
        given its yield to maturity. '''

        dd = self.dollarDuration(settlementDate, ytm)
        fp = self.fullPriceFromYield(settlementDate, ytm)
        md = dd / fp
        return md

###############################################################################

    def convexityFromYield(self, settlementDate, ytm):
        ''' Calculate the bond convexity from the yield to maturity. '''

        self.calculateFlowDates(settlementDate)
        dy = 0.0001
        p0 = self.fullPriceFromYield(settlementDate, ytm - dy)
        p1 = self.fullPriceFromYield(settlementDate, ytm)
        p2 = self.fullPriceFromYield(settlementDate, ytm + dy)
        conv = ((p2 + p0) - 2.0 * p1) / dy / dy / p1 / self._face
        return conv

###############################################################################

    def cleanPriceFromYield(self, settlementDate, ytm):
        ''' Calculate the bond clean price from the yield to maturity. '''

        fullPrice = self.fullPriceFromYield(settlementDate, ytm)
        accrued = self.accruedInterest(settlementDate)
        cleanPrice = fullPrice - accrued
        return cleanPrice

###############################################################################

    def fullPriceFromDiscountCurve(self, settlementDate, discountCurve):
        ''' Calculate the bond price using some discount curve to present-value
        the bond's cashflows. '''

        self.calculateFlowDates(settlementDate)
        pv = 0.0

        for dt in self._flowDates:
            df = discountCurve.df(dt)
            flow = self._coupon / self._frequency
            pv = pv + flow * df

        pv = pv + self._redemption * df
        return pv * self._face

##########################################################################

    def currentYield(self, cleanPrice):
        ''' Calculate the current yield of the bond which is the
        coupon divided by the clean price (not the full price)'''

        y = self._coupon * self._face / cleanPrice
        return y

##########################################################################

    def yieldToMaturity(self, settlementDate, cleanPrice):
        ''' Calculate the bond's yield to maturity by solving the price
        yield relationship using a one-dimensional root solver. '''

        self.calculateFlowDates(settlementDate)
        accrued = self.accruedInterest(settlementDate)
        fullPrice = cleanPrice + accrued
        ymin = -0.10

        pMax = self.fullPriceFromYield(settlementDate, ymin)

        if fullPrice > pMax:
            raise FinError("Entered price is too high.")

        argtuple = (self, settlementDate, fullPrice)

        ytm = optimize.newton(f,
                              x0=0.10,  # initial value of 10%
                              fprime=None,
                              args=argtuple,
                              tol=1e-8,
                              maxiter=50,
                              fprime2=None)

        return ytm

##########################################################################

    def accruedDays(self, settlementDate):
        ''' Calculate number of days from previous coupon date to settlement.'''
        self.calculateFlowDates(settlementDate)

        if len(self._flowDates) <= 2:
            raise FinError("Accrued interest - not enough flow dates.")

        return settlementDate - self.pcd(settlementDate)

##########################################################################

    def pcd(self, settlementDate):
        ''' Determine the previous coupon date before the settlement date. '''
        self.calculateFlowDates(settlementDate)

        if len(self._flowDates) <= 2:
            raise FinError("Accrued interest - not enough flow dates.")

        months = int(12 / self._frequency)
        ncd = self._flowDates[0]
        pcd = FinDate.addMonths(ncd, -months)

        return pcd

##########################################################################

    def accruedInterest(self, settlementDate):
        ''' Calculate the amount of coupon that has accrued between the
        previous coupon date and the settlement date. '''

        self.calculateFlowDates(settlementDate)

        if len(self._flowDates) <= 2:
            raise FinError("Accrued interest - not enough flow dates.")

        months = int(12 / self._frequency)
        ncd = self._flowDates[0]
        pcd = FinDate.addMonths(ncd, -months)

        f = self._frequency
        dt1 = pcd
        dt2 = settlementDate

        if dt1 > dt2:
            raise FinError("First coupon must precede settlement date.")

        d1 = dt1._d
        d2 = dt2._d
        m1 = dt1._m
        m2 = dt2._m
        y1 = dt1._y
        y2 = dt2._y

        if self._accrualType == FinBondAccruedTypes.THIRTY_360:
            dayDiff = 360 * (y2 - y1) + 30 * (m2 - m1) + (d2 - d1)
            daysInPeriod = 360.0 / f
            accFactor = dayDiff / daysInPeriod
        elif self._accrualType == FinBondAccruedTypes.THIRTY_360_BOND:
            d1 = min(d1, 30)
            if d1 == 31 or d1 == 30:
                d2 = min(d2, 30)
            dayDiff = 360 * (y2 - y1) + 30 * (m2 - m1) + (d2 - d1)
            daysInPeriod = 360.0 / f
            accFactor = dayDiff / daysInPeriod
        elif self._accrualType == FinBondAccruedTypes.ACT_360:
            accFactor = (dt2 - dt1) / 360.0
            daysInPeriod = 360.0 / f
            accFactor = (dt2 - dt1) / daysInPeriod
        elif self._accrualType == FinBondAccruedTypes.ACT_365:
            daysInPeriod = 365.0 / f
            accFactor = (dt2 - dt1) / daysInPeriod
        elif self._accrualType == FinBondAccruedTypes.ACT_ACT:
            daysInPeriod = ncd - pcd
            accFactor = (settlementDate - pcd) / daysInPeriod
        else:
            raise FinError(str(self._type) + " is not one of BondAccrualTypes")

        flow = self._coupon / f
        accrued = accFactor * flow * self._face
        return accrued

##########################################################################

    def initialPeriodFraction(self, settlementDate):
        ''' Calculate the amount of coupon that has accrued between the
        previous coupon date and the settlement date. '''

        self.calculateFlowDates(settlementDate)

        if len(self._flowDates) <= 2:
            raise FinError("Accrued interest - not enough flow dates.")

        months = int(12 / self._frequency)
        ncd = self._flowDates[0]

        dt1 = settlementDate
        dt2 = ncd

        d1 = dt1._d
        d2 = dt2._d
        m1 = dt1._m
        m2 = dt2._m
        y1 = dt1._y
        y2 = dt2._y
        accFactor = 0.0

        if self._accrualType == FinBondAccruedTypes.THIRTY_360:
            dayDiff = 360 * (y2 - y1) + 30 * (m2 - m1) + (d2 - d1)
            accFactor = dayDiff / 360.0
        elif self._accrualType == FinBondAccruedTypes.THIRTY_360_BOND:
            d1 = min(d1, 30)
            if d1 == 31 or d1 == 30:
                d2 = min(d2, 30)
            dayDiff = 360 * (y2 - y1) + 30 * (m2 - m1) + (d2 - d1)
            accFactor = dayDiff / 360.0
        elif self._accrualType == FinBondAccruedTypes.ACT_360:
            accFactor = (dt2 - dt1) / 360.0
        elif self._accrualType == FinBondAccruedTypes.ACT_365:
            accFactor = (dt2 - dt1) / 365.0
        elif self._accrualType == FinBondAccruedTypes.ACT_ACT:
            pcd = FinDate.addMonths(ncd, -months)
            numer = ncd - settlementDate
            denom = ncd - pcd
            accFactor = numer / denom
        else:
            raise FinError(str(self._type) + " is not one of BondAccrualTypes")

        return accFactor

##########################################################################

    def assetSwapSpread(
            self,
            settlementDate,
            cleanPrice,
            discountCurve,
            swapFloatDayCountConventionType=FinDayCountTypes.ACT_360,
            swapFloatFrequencyType=FinFrequencyTypes.SEMI_ANNUAL,
            swapFloatCalendarType=FinCalendarTypes.WEEKEND,
            swapFloatBusDayAdjustRuleType=FinBusDayConventionTypes.FOLLOWING,
            swapFloatDateGenRuleType=FinDateGenRuleTypes.BACKWARD):
        ''' Calculate the par asset swap spread of the bond. The discount curve
        is a Libor quality curve that is passed in. The price is the clean price.
        TODO - Check how first coupon on floating leg is sized. '''

        self.calculateFlowDates(settlementDate)

        accrued = self.accruedInterest(settlementDate)
        bondPrice = cleanPrice + accrued

        # Calculate the price of the bond discounted on the Libor curve
        pvLibor = 0.0
        prevDate = self._flowDates[0]

        for dt in self._flowDates:

            df = discountCurve.df(dt)
            pvLibor += df * self._coupon / self._frequency

        pvLibor += df * self._redemption

        # Calculate the PV01 of the floating leg of the asset swap
        # I assume here that the coupon starts accruing on the settlement date
        prevDate = self.pcd(settlementDate)
        schedule = FinSchedule(settlementDate,
                               self._maturityDate,
                               swapFloatFrequencyType,
                               swapFloatCalendarType,
                               swapFloatBusDayAdjustRuleType,
                               swapFloatDateGenRuleType)

        dayCount = FinDayCount(swapFloatDayCountConventionType)

        prevDate = schedule._adjustedDates[0]
        pv01 = 0.0
        for dt in schedule._adjustedDates[1:]:
            df = discountCurve.df(dt)
            yearFrac = dayCount.yearFrac(prevDate, dt)
            pv01 = pv01 + yearFrac * df
            prevDate = dt

        asw = (pvLibor - bondPrice) / pv01
        return asw

##########################################################################

    def fullPriceFromOAS(self,
                         settlementDate,
                         discountCurve,
                         oas):
        ''' Calculate the full price of the bond from its OAS given the bond
        settlement date, a discount curve and the oas as a number. '''

        self.calculateFlowDates(settlementDate)
        f = self._frequency
        c = self._coupon

        pv = 0.0
        for dt in self._flowDates:
            t = (dt - settlementDate) / gDaysInYear
            df = discountCurve.df(dt)
            # determine the Libor implied zero rate
            r = f * (pow(df, -1.0 / t / f) - 1.0)
            # determine the OAS adjusted zero rate
            df_adjusted = pow(1.0 + (r + oas) / f, -t * f)
            pv = pv + (c / f) * df_adjusted

        pv = pv + self._redemption * df_adjusted
        return pv * self._face

##########################################################################

    def optionAdjustedSpread(self,
                             settlementDate,
                             cleanPrice,
                             discountCurve):
        ''' Return OAS for bullet bond given settlement date, clean bond price
        and the discount relative to which the spread is to be computed. '''

        self.calculateFlowDates(settlementDate)
        accrued = self.accruedInterest(settlementDate)
        fullPrice = cleanPrice + accrued

        argtuple = (self, settlementDate, fullPrice, discountCurve)

        oas = optimize.newton(g,
                              x0=0.10,  # initial value of 10%
                              fprime=None,
                              args=argtuple,
                              tol=1e-8,
                              maxiter=50,
                              fprime2=None)

        return oas

##########################################################################

    def printFlows(self, settlementDate):
        ''' Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. '''
        pcd = self.pcd(settlementDate)
        print("PCD", pcd)
        datediff = settlementDate - pcd
        print("NOW", settlementDate, datediff)
        self.calculateFlowDates(settlementDate)
        for dt in self._flowDates:
            print("NCD", dt)

##########################################################################

    def priceFromSurvivalCurve(self,
                               discountCurve,
                               survivalCurve,
                               recoveryRate):
        ''' Calculate discounted present value of flows assuming default model.
        This has not been completed. '''
        pass

##########################################################################

    def print(self):
        ''' Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. '''
        print("Maturity Date:", self._maturityDate)
        print("Coupon:", self._coupon)
        print("Frequency:", self._frequencyType)
        print("Accrual:", self._accrualType)

##########################################################################
