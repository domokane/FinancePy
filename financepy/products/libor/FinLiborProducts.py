# -*- coding: utf-8 -*-
"""
Created on Sun Feb 07 14:23:13 2016

@author: Dominic O'Kane
"""

from ...finutils.FinDate import FinDate
from ...finutils.FinSchedule import schedule

###############################################################################

class FinLiborSwap(object):

    def __init__(self, startDate, endDate,
                 fixedCoupon, fixedFreq, fixedBasis,
                 floatSpread, floatFreq, floatBasis,
                 firstFixing=None,
                 payFixedFlag=True,
                 calendarName="WEEKEND",
                 businessDateAdjust="FOLLOWING",
                 dateGenRule="BACKWARD"):

        self.payFixedLeg = payFixedFlag
                
        self.fixedLeg = FinLiborSwapFixedLeg(startDate,
                                             endDate,
                                             fixedCoupon,
                                             fixedFreq,
                                             fixedBasis,
                                             calendarName,
                                             businessDateAdjust,
                                             dateGenRule)

        self.floatLeg = FinLiborSwapFloatLeg(startDate,
                                             endDate,
                                             floatSpread,
                                             floatFreq,
                                             floatBasis,
                                             firstFixing,
                                             calendarName,
                                             businessDateAdjust,
                                             dateGenRule)

        self.payFixedFlag = payFixedFlag

###############################################################################

    def value(self, valueDate, discountCurve, indexCurve):

        fixedLegValue = self.fixedLeg.value(valueDate,
                                            discountCurve)

        floatLegValue = self.floatLeg.value(valueDate,
                                            discountCurve,
                                            indexCurve)

        value = fixedLegValue - floatLegValue

        if self.payFixedLeg is True:
            value = value * (-1.0)

        return value

###############################################################################


    def dump(self):
        self.fixedLeg.dump()
#        self.floatLeg.dump()
                
###############################################################################
###############################################################################
###############################################################################


class FinLiborSwapFixedLeg(object):

    def __init__(self,
                 startDate,
                 maturityDate,
                 coupon,
                 freq,
                 basis,
                 calendarName="WEEKEND",
                 businessDateAdjust="MODIFIED_FOLLOWING",
                 dateGenRule="BACKWARD"):

        self.startDate = startDate
        self.maturityDate = maturityDate
        self.coupon = coupon
        self.freq = freq
        self.basis = basis

        self.flows = []

        self.schedule = schedule(self.startDate,
                                 self.maturityDate,
                                 self.freq,
                                 calendarName,
                                 businessDateAdjust,
                                 dateGenRule)

        self.generateFlows(self.basis)
        
###############################################################################

    def value(self, valueDate, discountCurve):
        
        df = discountCurve.df(valueDate)

        pv = 0.0

        for flow in self.flows[1:]:
            df = discountCurve.df(flow.date)
            pv += df * flow.amount

        return pv

###############################################################################

    def generateFlows(self, fixedBasis):

        if len(self.schedule) < 1:
            print("Error: Schedule has not been computed")
            return

        prevCpnDate = self.schedule[0]

        for dt in self.schedule[1:]:

            accruedPeriod = yearfrac(prevCpnDate, dt, fixedBasis)
            amount = accruedPeriod * self.coupon
            flow = FinFlow(dt, amount)
            self.flows.append(flow)
            prevCpnDate = dt

        numFlows = len(self.flows)

        self.flows[numFlows-1].amount += 1.0
                
        self.dump()

########################################################################

    def dump(self):

#        print "###################################################"

#        print "Swap Fixed Leg DUMP"
        
        if len(self.flows) < 1:
            print("Flows not calculated")
            return
        
        for flow in self.flows:

            flow.dump()

#        print "###################################################"
            
###############################################################################
########################################################################
########################################################################


class FinLiborSwapFloatLeg(object):

    def __init__(self,
                 startDate,
                 endDate,
                 floatSpread,
                 floatFreq,
                 floatBasis,
                 firstFixing,
                 calendarName,
                 businessDateAdjust,
                 dateGenRule):

        self.startDate = startDate
        self.endDate = endDate
        self.floatSpread = floatSpread
        self.freq = floatFreq
        self.basis = floatBasis
        self.firstFixing = firstFixing

        self.schedule = schedule(startDate,
                                  endDate,
                                  floatFreq,
                                  calendarName,
                                  businessDateAdjust,
                                  dateGenRule)
        
        self.flows = []
        
######################################################################

    def value(self, valueDate, discountCurve, indexCurve):

        self.generateFlows(indexCurve)

        z0 = discountCurve.df(valueDate)
        pv = 0.0

        for flow in self.flows[1:]:
            df = discountCurve.df(flow.date)
            pv += df * flow.amount

        self.dump()
        
        pv = pv / z0
        return pv

######################################################################

    def generateFlows(self, indexCurve):

        if len(self.schedule) <= 2:
            print("Schedule too short")
            return

        prevCpnDate = self.schedule[0]

        df1 = 1.0
        
        for dt in self.schedule[1:]:

            accruedPeriod = yearfrac(prevCpnDate, dt, self.basis)

            df2 = indexCurve.df(dt)

            liborFwd = (df1/df2-1.0) / accruedPeriod
            amount = accruedPeriod * liborFwd 
            flow = FinFlow(dt, amount)
            self.flows.append(flow)
            
            prevCpnDate = dt
 
            df1 = df2

        numFlows = len(self.flows)

        self.flows[numFlows-1].amount += 1.0

###############################################################################

    def dump(self):

#        print "###################################################")
        
        if len(self.flows) < 2:
            print("Flows not calculated")
            return
        
        for flow in self.flows:

            flow.dump()
                        
########################################################################
########################################################################

if __name__ == '__main__':

    print("==================================================================")    
    print("SEMI-ANNUAL FREQUENCY")
    print("==================================================================")    

    d1 = FinDate(20,6,2018)
    d2 = FinDate(20,6,2028)
    frequency = FinFrequency(2)
    calendar = FinCalendar("WEEKEND")
    businessDateAdjust = "FOLLOWING"
    dateGenRule = "BACKWARD"
    