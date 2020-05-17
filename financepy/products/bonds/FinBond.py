# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

#  - ROUNDING CONVENTIONS FOR ACCRUED
#  - CHECK OAS CALCULATION
#  - Check how first coupon on floating leg is sized on asset swaps. '''

# https://www.dmo.gov.uk/media/15004/convention_changes.pdf
###############################################################################
# Conventions
#  GILTS - SEMI ANNUAL ACT/ACT
#  US TREASURIES
###############################################################################

import numpy as np
from ...finutils.FinDate import FinDate
from ...finutils.FinError import FinError
from ...finutils.FinFrequency import FinFrequency, FinFrequencyTypes
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinSchedule import FinSchedule
from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinDayAdjustTypes
from ...finutils.FinCalendar import FinDateGenRuleTypes

from scipy import optimize

# References https://www.dmo.gov.uk/media/15011/yldeqns_v1.pdf
# DO TRUE YIELD
# JAPANESE SIMPLE YIELD

from enum import Enum


class FinYieldConventions(Enum):
    UK_DMO = 1,
    US_STREET = 2,
    US_TREASURY = 3

###############################################################################


def f(y, *args):
    ''' Function used to do root search in price to yield calculation. '''
    bond = args[0]
    settlementDate = args[1]
    price = args[2]
    convention = args[3]
    px = bond.fullPriceFromYield(settlementDate, y, convention)
    objFn = px - price
    return objFn

###############################################################################


def g(oas, *args):
    ''' Function used to do root search in price to OAS calculation. '''
    bond = args[0]
    settlementDate = args[1]
    price = args[2]
    discountCurve = args[3]
    px = bond.fullPriceFromOAS(settlementDate, discountCurve, oas)
    objFn = px - price
    return objFn

###############################################################################


class FinBond(object):
    ''' Class for fixed coupon bonds and performing related analytics. These
    are bullet bonds which means they have regular coupon payments of a known
    size that are paid on known dates plus a payment of par at maturity.'''

    def __init__(self,
                 maturityDate,
                 coupon,
                 frequencyType,
                 accrualType,
                 face=100.0):
        ''' Create FinBond object by providing Maturity Date, Frequency,
        coupon and the accrual convention type. '''

        if frequencyType not in FinFrequencyTypes:
            raise FinError("Invalid Frequency:" + str(frequencyType))
            return

        if accrualType not in FinDayCountTypes:
            raise FinError("Unknown Bond Accrued Convention type " +
                           str(accrualType))

        self._maturityDate = maturityDate
        self._coupon = coupon
        self._frequencyType = frequencyType
        self._accrualType = accrualType
        self._frequency = FinFrequency(frequencyType)
        self._face = face  # This is the position size
        self._par = 100.0  # This is how price is quoted

        ''' I do not determine cashflow dates as I do not want to require
        users to supply the issue date and without that I do not know how
        far to go back in the cashflow date schedule. '''

        self._flowDates = []
        self._settlementDate = FinDate(1, 1, 1900)
        self._accruedInterest = None
        self._accruedDays = 0.0
        self._alpha = 0.0

##########################################################################

    def calculateFlowDates(self, settlementDate):
        ''' Determine the bond cashflow payment dates. '''

        # No need to generate flows if settlement date has not changed
        if settlementDate == self._settlementDate:
            return

        self._settlementDate = settlementDate
        calendarType = FinCalendarTypes.NONE
        busDayRuleType = FinDayAdjustTypes.NONE
        dateGenRuleType = FinDateGenRuleTypes.BACKWARD

        self._flowDates = FinSchedule(settlementDate,
                                      self._maturityDate,
                                      self._frequencyType,
                                      calendarType,
                                      busDayRuleType,
                                      dateGenRuleType).generate()

        self._pcd = self._flowDates[0]
        self._ncd = self._flowDates[1]
        self.calcAccruedInterest(settlementDate)

###############################################################################

    def fullPriceFromYield(self, settlementDate, y,
                           convention=FinYieldConventions.UK_DMO):
        ''' Calculate the full price of bond from its yield to maturity. This
        function is vectorised with respect to the yield input. '''

        if convention not in FinYieldConventions:
            raise FinError("Yield convention unknown." + str(convention))

        y = np.array(y)  # VECTORIZED
        y = y + 0.000000000012345  # SNEAKY LOW-COST TRICK TO AVOID y=0

        self.calculateFlowDates(settlementDate)
        f = self._frequencyType.value
        c = self._coupon
        v = 1.0 / (1.0 + y/f)

        # n is the number of flows after the next coupon - we remove 2 because
        # the first element is the previous coupon date and then the ncd
        n = len(self._flowDates) - 2

        if convention == FinYieldConventions.UK_DMO:
            if n == 0:
                fp = (v**(self._alpha))*(1.0+c/f)
            else:
                term1 = c/f
                term2 = c*v/f
                term3 = c*v*v*(1.0-v**(n-1))/f/(1.0-v)
                term4 = v**n
                fp = (v**(self._alpha))*(term1 + term2 + term3 + term4)
        elif convention == FinYieldConventions.US_TREASURY:
            if n == 0:
                fp = (v**(self._alpha))*(1.0+c/f)
            else:
                term1 = c/f
                term2 = c*v/f
                term3 = c*v*v*(1.0-v**(n-1))/f/(1.0-v)
                term4 = v**n
                vw = 1.0 / (1.0 + self._alpha * y/f)
                fp = (vw)*(term1 + term2 + term3 + term4)
        elif convention == FinYieldConventions.US_STREET:
            vw = 1.0 / (1.0 + self._alpha * y/f)
            if n == 0:
                vw = 1.0 / (1.0 + self._alpha * y/f)
                fp = vw*(1.0+c/f)
            else:
                term1 = c/f
                term2 = c*v/f
                term3 = c*v*v*(1.0-v**(n-1))/f/(1.0-v)
                term4 = v**n
                fp = (v**(self._alpha))*(term1 + term2 + term3 + term4)
        else:
            raise ValueError("Unknown yield convention")

        return fp * self._par

###############################################################################

    def principal(self,
                  settlementDate,
                  y,
                  convention):
        ''' Calculate the principal value of the bond based on the face
        amount from its discount margin and making assumptions about the
        future Libor rates. '''

        fullPrice = self.fullPriceFromYield(settlementDate,
                                            y, convention)

        principal = fullPrice * self._face / self._par - self._accrued
        return principal

###############################################################################

    def dollarDuration(self,
                       settlementDate,
                       ytm,
                       convention=FinYieldConventions.UK_DMO):
        ''' Calculate the risk or dP/dy of the bond by bumping. '''

        self.calculateFlowDates(settlementDate)
        dy = 0.0001
        p0 = self.fullPriceFromYield(settlementDate, ytm - dy, convention)
        p2 = self.fullPriceFromYield(settlementDate, ytm + dy, convention)
        durn = -(p2 - p0) / dy / 2.0
        return durn

###############################################################################

    def macauleyDuration(self,
                         settlementDate,
                         ytm,
                         convention=FinYieldConventions.UK_DMO):
        ''' Calculate the Macauley duration of the bond on a settlement date
        given its yield to maturity. '''

        dd = self.dollarDuration(settlementDate, ytm, convention)
        fp = self.fullPriceFromYield(settlementDate, ytm, convention)
        md = dd * (1.0 + ytm / self._frequency) / fp
        return md

###############################################################################

    def modifiedDuration(self,
                         settlementDate,
                         ytm,
                         convention=FinYieldConventions.UK_DMO):
        ''' Calculate the modified duration of the bondon a settlement date
        given its yield to maturity. '''

        dd = self.dollarDuration(settlementDate, ytm, convention)
        fp = self.fullPriceFromYield(settlementDate, ytm, convention)
        md = dd / fp
        return md

###############################################################################

    def convexityFromYield(self,
                           settlementDate,
                           ytm,
                           convention=FinYieldConventions.UK_DMO):
        ''' Calculate the bond convexity from the yield to maturity. This
        function is vectorised with respect to the yield input. '''

        self.calculateFlowDates(settlementDate)
        dy = 0.0001
        p0 = self.fullPriceFromYield(settlementDate, ytm - dy, convention)
        p1 = self.fullPriceFromYield(settlementDate, ytm, convention)
        p2 = self.fullPriceFromYield(settlementDate, ytm + dy, convention)
        conv = ((p2 + p0) - 2.0 * p1) / dy / dy / p1 / self._par
        return conv

###############################################################################

    def cleanPriceFromYield(self, settlementDate, ytm,
                            convention=FinYieldConventions.UK_DMO):
        ''' Calculate the bond clean price from the yield to maturity. This
        function is vectorised with respect to the yield input. '''

        fullPrice = self.fullPriceFromYield(settlementDate, ytm, convention)
        cleanPrice = fullPrice - self._accruedInterest * self._par / self._face
        return cleanPrice

###############################################################################

    def cleanValueFromDiscountCurve(self, settlementDate, discountCurve):
        ''' Calculate the clean bond value using some discount curve to
        present-value the bond's cashflows back to the curve anchor date and
        not to the settlement date. '''

        fullPrice = self.valueBondUsingDiscountCurve(settlementDate,
                                                     discountCurve)

        accrued = self._accruedInterest * self._par / self._face
        cleanPrice = fullPrice - accrued
        return cleanPrice

##############################################################################

    def valueBondUsingDiscountCurve(self, settlementDate, discountCurve,
                                    verbose=False):
        ''' Calculate the bond *value* using some discount curve to PV the
        bond's cashflows to the curve anchor date. The anchor of the discount
        curve should be on the valuation date and so be 0-3 days before the
        settlement of the bond. This is not the same as the full price which
        is only the correct price on the settlement date of the bond which may
        be in the future.'''

        if discountCurve._curveDate > settlementDate:
            raise FinError("Discount curve date is after bond settlement date")

        self.calculateFlowDates(settlementDate)
        pv = 0.0

        for dt in self._flowDates[1:]:
            df = discountCurve.df(dt)
            flow = self._coupon / self._frequency
            pv = pv + flow * df

            if verbose is True:
                print(dt, flow, df, pv)

        pv = pv + df

        if verbose is True:
            print(dt, 1.0, df, pv)

        return pv * self._par

###############################################################################

    def currentYield(self, cleanPrice):
        ''' Calculate the current yield of the bond which is the
        coupon divided by the clean price (not the full price)'''

        y = self._coupon * self._par / cleanPrice
        return y

###############################################################################

    def yieldToMaturity(self,
                        settlementDate,
                        cleanPrice,
                        convention=FinYieldConventions.US_TREASURY):
        ''' Calculate the bond's yield to maturity by solving the price
        yield relationship using a one-dimensional root solver. '''

        if type(cleanPrice) is float or type(cleanPrice) is np.float64:
            cleanPrices = np.array([cleanPrice])
        elif type(cleanPrice) is list or type(cleanPrice) is np.ndarray:
            cleanPrices = np.array(cleanPrice)
        else:
            raise FinError("Unknown type for cleanPrice "
                           + str(type(cleanPrice)))

        self.calculateFlowDates(settlementDate)
        fullPrices = (cleanPrices + self._accruedInterest*self._par/self._face)
        ytms = []

        for fullPrice in fullPrices:

            argtuple = (self, settlementDate, fullPrice, convention)

            ytm = optimize.newton(f,
                                  x0=0.10,  # guess initial value of 10%
                                  fprime=None,
                                  args=argtuple,
                                  tol=1e-8,
                                  maxiter=50,
                                  fprime2=None)

            ytms.append(ytm)

        if len(ytms) == 1:
            return ytms[0]
        else:
            return np.array(ytms)

###############################################################################

    def calcAccruedInterest(self, settlementDate):
        ''' Calculate the amount of coupon that has accrued between the
        previous coupon date and the settlement date. '''

        if settlementDate != self._settlementDate:
            self.calculateFlowDates(settlementDate)

        if len(self._flowDates) == 0:
            raise FinError("Accrued interest - not enough flow dates.")

        dc = FinDayCount(self._accrualType)

        if self._accrualType == FinDayCountTypes.ACT_ACT_ICMA:
            accFactor = dc.yearFrac(self._pcd, settlementDate, self._ncd)
            alpha = 1.0 - accFactor
            accFactor = accFactor/self._frequency
        else:
            accFactor = dc.yearFrac(self._pcd, settlementDate)
            alpha = 1.0 - accFactor

        self._accruedInterest = accFactor * self._face * self._coupon
        self._alpha = alpha
        self._accruedDays = settlementDate - self._pcd

        return self._accruedInterest

###############################################################################

    def assetSwapSpread(
            self,
            settlementDate,
            cleanPrice,
            discountCurve,
            swapFloatDayCountConventionType=FinDayCountTypes.ACT_360,
            swapFloatFrequencyType=FinFrequencyTypes.SEMI_ANNUAL,
            swapFloatCalendarType=FinCalendarTypes.WEEKEND,
            swapFloatBusDayAdjustRuleType=FinDayAdjustTypes.FOLLOWING,
            swapFloatDateGenRuleType=FinDateGenRuleTypes.BACKWARD):
        ''' Calculate the par asset swap spread of the bond. The discount curve
        is a Libor curve that is passed in. This function is vectorised with
        respect to the clean price. '''

        cleanPrice = np.array(cleanPrice)
        self.calculateFlowDates(settlementDate)

        bondPrice = (cleanPrice + self._accruedInterest*self._par / self._face)
        # Calculate the price of the bond discounted on the Libor curve
        pvLibor = 0.0
        prevDate = self._pcd

        for dt in self._flowDates[1:]:
            df = discountCurve.df(dt)
            pvLibor += df * self._coupon / self._frequency
        pvLibor += df

        # Calculate the PV01 of the floating leg of the asset swap
        # I assume here that the coupon starts accruing on the settlement date
        prevDate = self._pcd
        schedule = FinSchedule(settlementDate,
                               self._maturityDate,
                               swapFloatFrequencyType,
                               swapFloatCalendarType,
                               swapFloatBusDayAdjustRuleType,
                               swapFloatDateGenRuleType)

        dayCount = FinDayCount(swapFloatDayCountConventionType)

        prevDate = self._pcd
        pv01 = 0.0
        for dt in schedule._adjustedDates[1:]:
            df = discountCurve.df(dt)
            yearFrac = dayCount.yearFrac(prevDate, dt)
            pv01 = pv01 + yearFrac * df
            prevDate = dt

        asw = (pvLibor - bondPrice/self._par) / pv01
        return asw

###############################################################################

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
        for dt in self._flowDates[1:]:
            t = (dt - settlementDate) / gDaysInYear
            df = discountCurve.df(dt)
            # determine the Libor implied zero rate
            r = f * (np.power(df, -1.0 / t / f) - 1.0)
            # determine the OAS adjusted zero rate
            df_adjusted = np.power(1.0 + (r + oas)/f, -t * f)
            pv = pv + (c / f) * df_adjusted

        pv = pv + df_adjusted
        pv *= self._par
        return pv

###############################################################################

    def optionAdjustedSpread(self,
                             settlementDate,
                             cleanPrice,
                             discountCurve):
        ''' Return OAS for bullet bond given settlement date, clean bond price
        and the discount relative to which the spread is to be computed. '''

        if type(cleanPrice) is float or type(cleanPrice) is np.float64:
            cleanPrices = np.array([cleanPrice])
        elif type(cleanPrice) is list or type(cleanPrice) is np.ndarray:
            cleanPrices = np.array(cleanPrice)
        else:
            raise FinError("Unknown type for cleanPrice "
                           + str(type(cleanPrice)))

        self.calculateFlowDates(settlementDate)
        fullPrices = (cleanPrices + self._accruedInterest*self._par/self._face)
        oass = []

        for fullPrice in fullPrices:

            argtuple = (self, settlementDate, fullPrice, discountCurve)

            oas = optimize.newton(g,
                                  x0=0.01,  # initial value of 1%
                                  fprime=None,
                                  args=argtuple,
                                  tol=1e-8,
                                  maxiter=50,
                                  fprime2=None)

            oass.append(oas)

        if len(oass) == 1:
            return oass[0]
        else:
            return np.array(oass)

###############################################################################

    def printFlows(self, settlementDate):
        ''' Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. '''
        self.calculateFlowDates(settlementDate)

        for dt in self._flowDates[1:-1]:
            flow = self._face * self._coupon/self._frequency
            print(dt, ",", flow)

        print(self._flowDates[-1], ",", self._face + flow)

###############################################################################

    def priceFromSurvivalCurve(self,
                               discountCurve,
                               survivalCurve,
                               recoveryRate):
        ''' Calculate discounted present value of flows assuming default model.
        This has not been completed. '''
        pass

###############################################################################

    def print(self):
        ''' Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. '''
        print("MATURITY DATE:", self._maturityDate)
        print("COUPON:", self._coupon)
        print("FREQUENCY:", self._frequencyType)
        print("ACCRUED TYPE:", self._accrualType)
        print("FACE:", self._face)

###############################################################################
