# -*- coding: utf-8 -*-
"""
Created on Sun Feb 07 14:23:13 2016

@author: Dominic O'Kane
"""

from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinBusDayConventionTypes, FinDateGenRuleTypes
from ...finutils.FinSchedule import FinSchedule
from ...finutils.FinMath import ONE_MILLION

##########################################################################


class FinLiborSwap(object):
    ''' Class for managing an interest rate swap contract. '''

    def __init__(self,
                 startDate,
                 maturityDate,
                 fixedCoupon,
                 fixedFreqType,
                 fixedDayCountType,
                 notional=ONE_MILLION,
                 floatSpread=0.0,
                 floatFreqType=FinFrequencyTypes.QUARTERLY,
                 floatDayCountType=FinDayCountTypes.THIRTY_360,
                 payFixedFlag=True,
                 calendarType=FinCalendarTypes.WEEKEND,
                 busDayAdjustType=FinBusDayConventionTypes.FOLLOWING,
                 dateGenRuleType=FinDateGenRuleTypes.BACKWARD):
        ''' Create an interest rate swap contract. '''
        if startDate > maturityDate:
            raise ValueError("Start date after maturity date")

        if fixedDayCountType not in FinDayCountTypes:
            raise ValueError(
                "Unknown Fixed Day Count Rule type " +
                str(fixedDayCountType))

        if floatDayCountType not in FinDayCountTypes:
            raise ValueError(
                "Unknown Float Day Count Rule type " +
                str(floatDayCountType))

        if fixedFreqType not in FinFrequencyTypes:
            raise ValueError(
                "Unknown Fixed Frequency type " +
                str(fixedFreqType))

        if floatFreqType not in FinFrequencyTypes:
            raise ValueError(
                "Unknown Float Frequency type " +
                str(fixedFreqType))

        if calendarType not in FinCalendarTypes:
            raise ValueError("Unknown Calendar type " + str(calendarType))

        if busDayAdjustType not in FinBusDayConventionTypes:
            raise ValueError(
                "Unknown Business Day Adjust type " +
                str(busDayAdjustType))

        if dateGenRuleType not in FinDateGenRuleTypes:
            raise ValueError(
                "Unknown Date Gen Rule type " +
                str(dateGenRuleType))

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

        # These are generated immediately as they are for the entire
        # life of the swap. Given a valuation date we can determine
        # which cash flows are in the future and value the swap
        self.generateFixedLegPaymentDates()
        self.generateFloatLegPaymentDates()

##########################################################################

    def value(self,
              valuationDate,
              discountCurve,
              indexCurve,
              firstFixingRate,
              principal=0.0):
        ''' Value the interest rate swap on a value date given a single Libor
        discount curve. '''

        fixedLegValue = self.fixedLegValue(valuationDate,
                                           discountCurve,
                                           principal)

        floatLegValue = self.floatLegValue(valuationDate,
                                           discountCurve,
                                           indexCurve,
                                           firstFixingRate,
                                           principal)

        value = fixedLegValue - floatLegValue

        if self._payFixedLeg is True:
            value = value * (-1.0)

        return value * self._notional

##########################################################################

    def generateFixedLegPaymentDates(self):
        ''' Generate the fixed leg payment dates all the way back to
        the start date of the swap which may precede the valuation date'''
        self._adjustedFixedDates = FinSchedule(
            self._startDate,
            self._maturityDate,
            self._fixedFrequencyType,
            self._calendarType,
            self._busDayAdjustType,
            self._dateGenRuleType).generate()

##########################################################################

    def generateFloatLegPaymentDates(self):
        ''' Generate the floating leg payment dates all the way back to
        the start date of the swap which may precede the valuation date'''
        self._adjustedFloatDates = FinSchedule(
            self._startDate,
            self._maturityDate,
            self._floatFrequencyType,
            self._calendarType,
            self._busDayAdjustType,
            self._dateGenRuleType).generate()

##########################################################################

    def pv01(self, valuationDate, discountCurve):
        ''' Calculate the value of 1 basis point coupon on the fixed leg. '''

        pv = self.fixedLegValue(valuationDate, discountCurve)
        pv01 = pv / self._fixedCoupon
        return pv01

##########################################################################

    def parCoupon(self, valuationDate, discountCurve):
        ''' Calculate the fixed leg coupon that makes the swap worth zero. '''

        pv01 = self.pv01(valuationDate, discountCurve)
        df0 = discountCurve.df(valuationDate)
        dfT = discountCurve.df(self._maturityDate)
        cpn = (df0 - dfT) / pv01
        return cpn

##########################################################################

    def fixedLegValue(self, valuationDate, discountCurve, principal=0.0):

        self._fixedYearFracs = []
        self._fixedFlows = []
        self._fixedDfs = []
        self._fixedFlowPVs = []

        basis = FinDayCount(self._fixedDayCountType)

        ''' The swap may have started in the past but we can only value
        payments that have occurred after the valuation date. '''
        startIndex = 0
        while self._adjustedFixedDates[startIndex] < valuationDate:
            startIndex += 1

        ''' If the swap has yet to settle then we do not include the
        start date of the swap as a coupon payment date. '''
        if valuationDate <= self._startDate:
            startIndex = 1

        self._fixedStartIndex = startIndex

        ''' Now PV fixed leg flows. '''
        pv = 0.0
        prevDt = self._adjustedFixedDates[startIndex - 1]

        for nextDt in self._adjustedFixedDates[startIndex:]:
            alpha = basis.yearFrac(prevDt, nextDt)
            df_discount = discountCurve.df(nextDt)
            flow = self._fixedCoupon * alpha
            pv += flow * df_discount
            prevDt = nextDt

            self._fixedYearFracs.append(alpha)
            self._fixedFlows.append(flow)
            self._fixedDfs.append(df_discount)
            self._fixedFlowPVs.append(flow * df_discount)

        pv = pv + principal * df_discount
        self._fixedFlowPVs[-1] += principal * df_discount
        self._fixedFlows[-1] += principal

        z0 = discountCurve.df(valuationDate)
        pv = pv / z0
        return pv

##########################################################################

    def floatLegValue(self,
                      valuationDate,
                      discountCurve,
                      indexCurve,
                      firstFixingRate=None,
                      principal=0.0):
        ''' Value the floating leg with payments from an index curve and
        discounting based on a supplied discount curve. '''

        self._floatYearFracs = []
        self._floatFlows = []
        self._floatFlowPVs = []
        self._floatDfs = []

        basis = FinDayCount(self._floatDayCountType)
        pv = 0.0

        ''' The swap may have started in the past but we can only value
        payments that have occurred after the start date. '''
        startIndex = 0
        while self._adjustedFloatDates[startIndex] < valuationDate:
            startIndex += 1

        ''' If the swap has yet to settle then we do not include the
        start date of the swap as a coupon payment date. '''
        if valuationDate <= self._startDate:
            startIndex = 1

        self._floatStartIndex = startIndex

        ''' The first floating payment is usually already fixed so is
        not implied by the index curve. '''
        prevDt = self._adjustedFloatDates[startIndex - 1]
        nextDt = self._adjustedFloatDates[startIndex]
        alpha = basis.yearFrac(prevDt, nextDt)
        df1_index = 1.0  # Cannot be previous date as that has past
        df2_index = indexCurve.df(nextDt)

        if firstFixingRate is None:
            libor = (df1_index / df2_index - 1.0) / alpha
            flow = libor * alpha
        else:
            flow = firstFixingRate * alpha

        df_discount = discountCurve.df(nextDt)
        pv += flow * df_discount

        self._floatYearFracs.append(alpha)
        self._floatFlows.append(flow)
        self._floatDfs.append(df_discount)
        self._floatFlowPVs.append(flow * df_discount)

        prevDt = nextDt
        df1_index = indexCurve.df(prevDt)

        for nextDt in self._adjustedFloatDates[startIndex + 1:]:
            alpha = basis.yearFrac(prevDt, nextDt)
            df2_index = indexCurve.df(nextDt)
            flow = (df1_index / df2_index - 1.0)  # The accrual factors cancel
            df_discount = discountCurve.df(nextDt)
            pv += flow * df_discount
            df1_index = df2_index
            prevDt = nextDt

            self._floatFlows.append(flow)
            self._floatYearFracs.append(alpha)
            self._floatDfs.append(df_discount)
            self._floatFlowPVs.append(flow * df_discount)

        pv = pv + principal * df_discount
        self._floatFlows[-1] += principal
        self._floatFlowPVs[-1] += principal * df_discount
        z0 = discountCurve.df(valuationDate)
        pv = pv / z0

        return pv

##########################################################################

    def printFixedLeg(self, valuationDate):
        ''' Prints the fixed leg amounts. '''

        print("START DATE:", self._startDate)
        print("MATURITY DATE:", self._maturityDate)
        print("COUPON (%):", self._fixedCoupon * 100)
        print("FIXED LEG FREQUENCY:", str(self._fixedFrequencyType))
        print("FIXED LEG DAY COUNT:", str(self._fixedDayCountType))
        print("VALUATION DATE", valuationDate)

        if self._fixedFlows is None:
            print("Fixed Flows not calculated.")
            return

        print("PAYMENT DATE", "AMOUNT")
        numFlows = len(self._adjustedFixedDates)
        totalPV = 0.0

        for i in range(self._fixedStartIndex, numFlows):
            paymentDate = self._adjustedFixedDates[i]
            iFlow = i - self._fixedStartIndex
            flow = self._fixedFlows[iFlow] * self._notional
            alpha = self._fixedYearFracs[iFlow]
            df = self._fixedDfs[iFlow]
            flowPV = self._fixedFlowPVs[iFlow] * self._notional
            totalPV += flowPV
            print("%s %10.7f %12.0f %12.6f %12.0f %12.0f" %
                  (paymentDate, alpha, flow, df, flowPV, totalPV))

##########################################################################

    def printFloatLeg(self, valuationDate):
        ''' Prints the floating leg amounts. '''

        print("START DATE:", self._startDate)
        print("MATURITY DATE:", self._maturityDate)
        print("SPREAD COUPON (%):", self._floatSpread * 100)
        print("FLOAT LEG FREQUENCY:", str(self._floatFrequencyType))
        print("FLOAT LEG DAY COUNT:", str(self._floatDayCountType))
        print("VALUATION DATE", valuationDate)

        if self._floatFlows is None:
            print("Floating Flows not calculated.")
            return

        print("PAYMENT DATE", "AMOUNT")
        numFlows = len(self._adjustedFloatDates)
        totalPV = 0.0

        for i in range(self._floatStartIndex, numFlows):
            paymentDate = self._adjustedFloatDates[i]
            iFlow = i - self._floatStartIndex
            flow = self._floatFlows[iFlow] * self._notional
            alpha = self._floatYearFracs[iFlow]
            df = self._floatDfs[iFlow]
            flowPV = self._floatFlowPVs[iFlow] * self._notional
            totalPV += flowPV
            print("%s %10.7f %12.0f %12.6f %12.0f %12.0f" %
                  (paymentDate, alpha, flow, df, flowPV, totalPV))

##########################################################################
