# -*- coding: utf-8 -*-
"""
Created on Sun Feb 07 14:23:13 2016

@author: Dominic O'Kane
"""

from ...finutils.FinDayCount import FinDayCount
from ...finutils.FinDayCount import FinDayCountTypes
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinBusDayConventionTypes, FinDateGenRuleTypes
from ...finutils.FinSchedule import FinSchedule
from ...finutils.FinMath import ONE_MILLION

################################################################################

class FinLiborSwap(object):
    ''' Class for managing an interest rate swap contract. '''

    def __init__(self, 
                 startDate, 
                 maturityDate,
                 fixedCoupon, 
                 fixedFreqType, 
                 fixedDayCountType,
                 notional = ONE_MILLION,
                 floatSpread = 0.0, 
                 floatFreqType = FinFrequencyTypes.QUARTERLY, 
                 floatDayCountType = FinDayCountTypes.THIRTY_360,
                 firstFixing = 0.0,
                 payFixedFlag=True,
                 calendarType = FinCalendarTypes.WEEKEND,
                 busDayAdjustType = FinBusDayConventionTypes.FOLLOWING,
                 dateGenRuleType=FinDateGenRuleTypes.BACKWARD):

        ''' Create an interest rate swap contract. '''
        if startDate > maturityDate:
            raise ValueError("Start date after maturity date")

        if fixedDayCountType not in FinDayCountTypes:
            raise ValueError("Unknown Fixed Day Count Rule type " + str(fixedDayCountType))

        if floatDayCountType not in FinDayCountTypes:
            raise ValueError("Unknown Float Day Count Rule type " + str(floatDayCountType))

        if fixedFreqType not in FinFrequencyTypes:
            raise ValueError("Unknown Fixed Frequency type " + str(fixedFreqType))

        if floatFreqType not in FinFrequencyTypes:
            raise ValueError("Unknown Float Frequency type " + str(fixedFreqType))

        if calendarType not in FinCalendarTypes:
            raise ValueError("Unknown Calendar type " + str(calendarType))

        if busDayAdjustType not in FinBusDayConventionTypes:
            raise ValueError("Unknown Business Day Adjust type " + str(busDayAdjustType))

        if dateGenRuleType not in FinDateGenRuleTypes:
            raise ValueError("Unknown Date Gen Rule type " + str(dateGenRuleType))

        self._maturityDate = maturityDate
        self._payFixedLeg = payFixedFlag
        self._notional = notional
        self._startDate = startDate

        self._fixedCoupon = fixedCoupon
        self._floatSpread = floatSpread

        self._fixedFrequencyType = fixedFreqType
        self._floatFrequencyType = floatFreqType

        self._fixedDayCountType = fixedDayCountType
        self._floatDayCountType = floatDayCountType

        self._payFixedFlag = payFixedFlag

        self._calendarType = calendarType
        self._busDayAdjustType = busDayAdjustType
        self._dateGenRuleType = dateGenRuleType

        self._fixedFlows = []
        self._floatingFlows = []

        self._adjustedFixedDates = []
        self._adjustedFloatDates = []

################################################################################

    def value(self, valueDate, discountCurve):
        ''' Value the interest rate swap on a value date given a single Libor 
        discount curve. '''
        principal = 1.0

        fixedLegValue = self.fixedLegValue(valueDate,
                                            discountCurve, 
                                            principal)

        floatLegValue = self.floatLegValue(valueDate,
                                            discountCurve, 
                                            discountCurve, 
                                            principal)

        value = fixedLegValue - floatLegValue

        if self._payFixedLeg is True:
            value = value * (-1.0)

        return value * self._notional

################################################################################

    def generateFixedLegFlows(self, valueDate):

        self._adjustedFixedDates = FinSchedule(self._startDate,
                                               self._maturityDate,
                                               self._fixedFrequencyType,
                                               self._calendarType,
                                               self._busDayAdjustType,
                                               self._dateGenRuleType).generate()

        dayCounter = FinDayCount(self._fixedDayCountType)
        self._fixedFlows = []
        d1 = valueDate
        for dt in self._adjustedFixedDates:
            flow = dayCounter.yearFrac(d1,dt) * self._fixedCoupon
            self._fixedFlows.append(flow)
            d1 = dt

################################################################################

    def fixedLegValue(self, valueDate, discountCurve, principal = 0.0):

        self.generateFixedLegFlows(self._startDate)
        
        pv = 0.0
        df = 1.0

        for dt, flow in zip(self._adjustedFixedDates, self._fixedFlows):
            df = discountCurve.df(dt)
            pv += df * flow

        pv = pv + principal * df
        z0 = discountCurve.df(valueDate)
        pv = pv / z0
        return pv

################################################################################

    def pv01(self, valueDate, discountCurve, principal = 0.0):

        if self._fixedFlows == []:
            self.generateFixedLegFlows(self._startDate)

        pv = 0.0
        df = 1.0

        for dt, flow in zip(self._adjustedFixedDates, self._fixedFlows):
            df = discountCurve.df(dt)
            pv += df * flow

        z0 = discountCurve.df(valueDate)
        pv01 = pv / z0 / self._fixedCoupon
        return pv01

################################################################################

    def parCoupon(self, valueDate, discountCurve, principal = 0.0):

        if self._fixedFlows == []:
            self.generateFixedLegFlows(self._startDate)

        pv01 = self.pv01(valueDate,discountCurve,principal)
        df0 = discountCurve.df(self._startDate)
        dfT = discountCurve.df(self._maturityDate)        
        cpn = (df0 - dfT) / (abs(pv01) + 1e-20)
        return cpn

################################################################################

    def printFixedLeg(self):
        ''' Prints the fixed leg amounts. '''
        print("START DATE:",self._startDate)
        print("MATURITY DATE:",self._maturityDate)
        print("COUPON (%):", self._fixedCoupon*100)
        print("FIXED LEG FREQUENCY:",str(self._fixedFrequencyType))
        print("FIXED LEG DAY COUNT:",str(self._fixedDayCountType))

        dayCounter = FinDayCount(self._fixedDayCountType)
        numFlows = len(self._adjustedFixedDates)

        print("PAYMENT DATE","AMOUNT")
        for i in range(1,numFlows):
            dt1 = self._adjustedFixedDates[i-1]
            dt2 = self._adjustedFixedDates[i]
            acc = dayCounter.yearFrac(dt1,dt2)
            print(dt2,acc*self._coupon)

################################################################################

    def generateFloatingLegFlows(self, valueDate, indexCurve):
        ''' Generate the payment amounts on floating leg implied by index curve '''
        self._adjustedFloatingDates = FinSchedule(self._startDate,
                                          self._maturityDate,
                                          self._floatFrequencyType,
                                          self._calendarType,
                                          self._busDayAdjustType,
                                          self._dateGenRuleType).generate()

        basis = FinDayCount(self._floatDayCountType)

        dt1 = valueDate
        df1 = indexCurve.df(dt1)

        self._floatingFlows = []

        for dt2 in self._adjustedFloatingDates[1:]:

            alpha = basis.yearFrac(dt1,dt2)
            df2 = indexCurve.df(dt2)
            flow =  (df1/df2-1.0) / alpha
            dt1 = dt2
            df1 = df2

            self._floatingFlows.append(flow)

################################################################################

    def floatLegValue(self, 
                          valueDate, 
                          discountCurve, 
                          indexCurve, 
                          principal = 0.0 ):
        ''' Value the floating leg with payments from an index curve and 
        discounting based on a supplied discount curve. '''
        basis = FinDayCount(self._floatDayCountType)

        if self._floatingFlows == []:
            self.generateFloatingLegFlows(valueDate,indexCurve)

        dt1 = valueDate
        pv = 0.0

        for dt2, flow in zip(self._adjustedFloatingDates,self._floatingFlows):
            df = discountCurve.df(dt2)
            alpha = basis.yearFrac(dt1,dt2)
            pv += df * flow * alpha
            dt1 = dt2

        pv = pv + df * principal
        return pv

################################################################################
