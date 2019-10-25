# -*- coding: utf-8 -*-
"""
Created on Sun Aug 07 14:23:13 2019

@author: Dominic O'Kane
"""

from math import log, sqrt

from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinBusDayConventionTypes, FinDateGenRuleTypes
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinMath import ONE_MILLION, N
from ...finutils.FinError import FinError
from ...finutils.FinSchedule import FinSchedule
from ...models.FinSABRModel import blackVolFromSABR

################################################################################

from enum import Enum

class FinLiborCapFloorType(Enum):
    CAP = 1
    FLOOR = 2

class FinLiborCapFloorModelTypes(Enum):
    BLACK = 1
    SHIFTED_BLACK = 2
    SABR = 3

################################################################################
    
class FinLiborCapFloor():

    def __init__(self, 
                 startDate,
                 maturityDate,
                 optionType,
                 strikeRate,
                 lastFixing = None,
                 frequencyType = FinFrequencyTypes.QUARTERLY,
                 dayCountType = FinDayCountTypes.THIRTY_E_360_ISDA,
                 notional = ONE_MILLION, 
                 calendarType = FinCalendarTypes.WEEKEND,
                 busDayAdjustType = FinBusDayConventionTypes.FOLLOWING,
                 dateGenRuleType = FinDateGenRuleTypes.BACKWARD):

        if startDate > maturityDate: 
            raise FinError("Start date must be before maturity date")

        if optionType not in FinLiborCapFloorType:
            raise FinError("Unknown Libor Cap Floor type " + str(optionType))

        if dayCountType not in FinDayCountTypes:
            raise FinError("Unknown Cap Floor DayCountRule type " + str(dayCountType))

        if frequencyType not in FinFrequencyTypes:
            raise FinError("Unknown CapFloor Frequency type " + str(frequencyType))

        if calendarType not in FinCalendarTypes:
            raise FinError("Unknown Calendar type " + str(calendarType))

        if busDayAdjustType not in FinBusDayConventionTypes:
            raise FinError("Unknown Business Day Adjust type " + str(busDayAdjustType))

        if dateGenRuleType not in FinDateGenRuleTypes:
            raise FinError("Unknown Date Gen Rule type " + str(dateGenRuleType))

        self._startDate = startDate
        self._maturityDate = maturityDate
        self._optionType = optionType
        self._strikeRate = strikeRate
        self._lastFixing = lastFixing
        self._frequencyType = frequencyType
        self._dayCountType = dayCountType
        self._notional = notional
        self._calendarType = calendarType
        self._busDayAdjustType = busDayAdjustType
        self._dateGenRuleType = dateGenRuleType
        
################################################################################
        
    def value(self, 
              valuationDate,
              liborCurve, 
              modelType, 
              modelParams):

        self._capFloorDates = FinSchedule(self._startDate,
                                          self._maturityDate,
                                          self._frequencyType,
                                          self._calendarType,
                                          self._busDayAdjustType,
                                          self._dateGenRuleType).generate()

        dayCounter = FinDayCount(self._dayCountType)
        numOptions = len(self._capFloorDates)
        strikeRate = self._strikeRate

#        for dt in self._capFloorDates:
#            print(dt)
            
        if strikeRate <= 0.0:
            raise FinError("Strike <= 0.0")

        if numOptions <= 1:
            raise FinError("Number of options in capfloor equals 1")
            
        ########################################################################
        
        capFloorValue = 0.0
        
        # Value the first caplet or floorlet with known payoff
        
        if self._lastFixing is None:

            # Assume that the fixing is set today and that there is 
            # potentially some intrinsic value depending on strike
            startDate = self._startDate
            endDate = self._capFloorDates[1]
            fwdRate = liborCurve.fwd(startDate,endDate,self._dayCountType)
            alpha = dayCounter.yearFrac(startDate,endDate)
            df = liborCurve.df(endDate)
                 
            if self._optionType == FinLiborCapFloorType.CAP:
                capFloorLetValue = df * alpha * max(fwdRate - strikeRate,0)
            elif self._optionType == FinLiborCapFloorType.FLOOR:
                capFloorLetValue = df * alpha * max(strikeRate - fwdRate,0)

        else:
            
            startDate = self._startDate
            endDate = self._capFloorDates[1]
            fwdRate = self._lastFixing
            df = liborCurve.df(endDate)
            alpha = dayCounter.yearFrac(startDate,endDate)
                 
            if self._optionType == FinLiborCapFloorType.CAP:
                capFloorLetValue = df * alpha * max(fwdRate - strikeRate,0)
            elif self._optionType == FinLiborCapFloorType.FLOOR:
                capFloorLetValue = df * alpha * max(strikeRate - fwdRate,0)
    
        capFloorValue += capFloorLetValue
               
        for i in range(2,numOptions):
            
            startDate = self._capFloorDates[i-1]
            endDate = self._capFloorDates[i]

            capFloorLetValue = self.valueCapletFloorlet(valuationDate,
                                                   startDate, 
                                                   endDate, 
                                                   liborCurve, 
                                                   modelType, 
                                                   modelParams)
    
            capFloorValue +=capFloorLetValue 
            
        capFloorValue = capFloorValue * self._notional
        return capFloorValue


################################################################################

    def valueCapletFloorlet(self,
                            valuationDate,
                            startDate, 
                            endDate, 
                            liborCurve, 
                            modelType, 
                            modelParams):

        df = liborCurve.df(endDate)
        t = (startDate - valuationDate)/gDaysInYear
        f = liborCurve.fwd(startDate, endDate, self._dayCountType)
        k = self._strikeRate
        
        if modelType == FinLiborCapFloorModelTypes.BLACK:

            v = modelParams["volatility"]

            if v <= 0:
                raise FinError("Black Volatility must be positive")
            
            d1 = (log(f/k) + v*v*t/2.0)/v/sqrt(t)
            d2 = d1 - v * sqrt(t)

            if self._optionType == FinLiborCapFloorType.CAP:
                capFloorLetValue = df * (f * N(+d1) - k * N(+d2))
            elif self._optionType == FinLiborCapFloorType.FLOOR:
                capFloorLetValue = df * (k * N(-d2) - f * N(-d1))

        elif modelType == FinLiborCapFloorModelTypes.SHIFTED_BLACK:

            v = modelParams["volatility"]
            h = modelParams["shift"]

            if v <= 0:
                raise FinError("Black Volatility must be positive")
            
            d1 = (log((f-h)/(k-h)) + v*v*t/2.0)/v/sqrt(t)
            d2 = d1 - v * sqrt(t)

            if self._optionType == FinLiborCapFloorType.CAP:
                capFloorLetValue = df * ((f-h) * N(+d1) - (k-h) * N(+d2))
            elif self._optionType == FinLiborCapFloorType.FLOOR:
                capFloorLetValue = df * ((k-h) * N(-d2) - (f-h) * N(-d1))

        elif modelType == FinLiborCapFloorModelTypes.SABR:
            
            alpha = modelParams["alpha"]
            beta = modelParams["beta"]
            rho = modelParams["rho"]
            nu = modelParams["nu"]
            
            v = blackVolFromSABR(alpha,beta,rho,nu,f,k,t)
                        
            d1 = (log((f)/(k)) + v*v*t/2.0)/v/sqrt(t)
            d2 = d1 - v * sqrt(t)

            if self._optionType == FinLiborCapFloorType.CAP:
                capFloorLetValue = df * ((f) * N(+d1) - (k) * N(+d2))
            elif self._optionType == FinLiborCapFloorType.FLOOR:
                capFloorLetValue = df * ((k) * N(-d2) - (f) * N(-d1))
                
        else:
            raise FinError("Unknown model type " + str(modelType))

        return capFloorLetValue
    
################################################################################

    def print(self):
        print("CAP FLOOR START DATE:",self._startDate)
        print("CAP FLOOR MATURITY DATE:",self._maturityDate)
        print("CAP FLOOR STRIKE COUPON:", self._strikeRate*100)
        print("CAP FLOOR OPTION TYPE:",str(self._optionType))
        print("CAP FLOOR FREQUENCY:",str(self._frequencyType))
        print("CAP FLOOR DAY COUNT:",str(self._dayCountType))
        
################################################################################
