##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO: Extend to allow term structure of volatility
# TODO: Extend to allow two fixed legs in underlying swap
# TODO: Cash settled swaptions

''' This module implements the LMM in the spot measure. It combines both model
and product specific code - I am not sure if it is better to separate these. At
the moment this seems to work ok. '''

import numpy as np

from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinBusDayAdjustTypes
from ...finutils.FinCalendar import FinDateGenRuleTypes
from ...finutils.FinDayCount import FinDayCountTypes
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinError import FinError
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes
from ...finutils.FinDate import FinDate
from ...products.libor.FinLiborSwap import FinLiborSwap

from ...models.FinModelBlack import FinModelBlack

from ...finutils.FinDate import FinDate
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinCalendar import FinCalendarTypes,  FinDateGenRuleTypes
from ...finutils.FinCalendar import FinBusDayAdjustTypes
from ...finutils.FinSchedule import FinSchedule

from ...models.FinModelRatesLMM import *

from ...finutils.FinOptionTypes import FinOptionTypes

##########################################################################

from enum import Enum


class FinLiborSwaptionTypes(Enum):
    PAYER = 1
    RECEIVER = 2

##########################################################################


class FinLiborLMMProducts():
    ''' This is the class for pricing Libor products using the LMM. '''

    def __init__(self,
                 settlementDate: FinDate,
                 exerciseDate: FinDate,
                 maturityDate: FinDate,
                 floatFrequencyType: FinFrequencyTypes = FinFrequencyTypes.QUARTERLY,
                 floatDayCountType: FinDayCountTypes = FinDayCountTypes.THIRTY_360,
                 calendarType: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                 busDayAdjustType: FinBusDayAdjustTypes = FinBusDayAdjustTypes.FOLLOWING,
                 dateGenRuleType: FinDateGenRuleTypes = FinDateGenRuleTypes.BACKWARD):
        ''' Create a European-style swaption by defining the exercise date of
        the swaption, and all of the details of the underlying interest rate
        swap including the fixed coupon and the details of the fixed and the
        floating leg payment schedules. '''

        checkArgumentTypes(self.__init__, locals())

        if settlementDate > exerciseDate:
            raise FinError("Settlement date must be before expiry date")

        if exerciseDate > maturityDate:
            raise FinError("Exercise date must be before swap maturity date")

        ''' Set up the grid for the Libor rates that are to be simulated. These
        must be consistent with the floating rate leg of the product that is to
        be priced. '''

        self._gridDates = FinSchedule(settlementDate,
                                      maturityDate,
                                      floatFrequencyType,
                                      calendarType,
                                      busDayAdjustType,
                                      dateGenRuleType).generate()

        self._floatDayCountType = floatDayCountType

        basis = FinDayCount(self._fixedDayCountType)

        prevDt = self._gridDates[0]
        for nextDt in self._gridDates[1:]:
            tau = basis.yearFrac(prevDt, nextDt)
            self._accrualFactors.append(tau)
            prevDt = nextDt

        self._fwds = None

##########################################################################

    def simulate(self,
                 discountCurve,
                 volCurves,
                 correlationMatrix: np.ndarray,
                 modelType: FinRateModelLMMModelTypes,
                 numPaths: int = 1000,
                 numeraireIndex: int = 0,
                 useSobol: bool = True,
                 seed: int = 42):
        ''' Run the simulation to generate and store all of the Libor forward
        rate paths. '''

        if numPaths < 2 or numPaths > 1000000:
            raise FinError("NumPaths must be between 2 and 1 million")

        if isinstance(modelType, FinRateModelLMMModelTypes) is False:
            raise FinError("Model type must be type FinRateModelLMMModelTypes")

        if discountCurve.curveDate != self._startDate:
            raise FinError("Curve anchor date not the same as LMM start date.")

        self._numPaths = numPaths
        self._volCurves = volCurves
        self._correlationMatrix = correlationMatrix
        self._modelType = modelType
        self._numeraireIndex = numeraireIndex
        self._useSobol = useSobol

        self._numForwards = len(self._gridDates) - 1
        self._forwardCurve = []

        for i in range(1, self._numForwards):
            startDate = self._gridDates[i-1]
            endDate = self._gridDates[i]
            fwdRate = discountCurve.forwardRate(startDate, endDate,
                                                self._floatDayCountType)
            self._forwardCurve.append(fwdRate)

        if modelType == FinLMMModelTypes.LMM_ONE_FACTOR:

            self._fwds = LMMSimulateFwds1F(numForwards,
                                           numPaths,
                                           numeraireIndex,
                                           self._forwardCurve,
                                           gammas,
                                           self._accrualFactors,
                                           useSobol,
                                           seed)

        elif modelType == FinLMMModelTypes.LMM_HW_M_FACTOR:
            pass
        elif modelType == FinLMMModelTypes.LMM_FULL_N_FACTOR:
            pass        
        else:
            raise FinError("Unknown LMM Model Type")


###############################################################################

     def valueSwaption(self, 
                       settlementDate: FinDate,
                       exerciseDate: FinDate,
                       maturityDate: FinDate,
                       swaptionType: FinLiborSwaptionTypes,
                       fixedCoupon: float,
                       fixedFrequencyType: FinFrequencyTypes,
                       fixedDayCountType: FinDayCountTypes,
                       notional: float = ONE_MILLION,
                       floatFrequencyType: FinFrequencyTypes = FinFrequencyTypes.QUARTERLY,
                       floatDayCountType: FinDayCountTypes = FinDayCountTypes.THIRTY_360,
                       calendarType: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                       busDayAdjustType: FinBusDayAdjustTypes = FinBusDayAdjustTypes.FOLLOWING,
                       dateGenRuleType: FinDateGenRuleTypes = FinDateGenRuleTypes.BACKWARD):

        swaptionFloatDates = FinSchedule(settlementDate,
                                         maturityDate,
                                         floatFrequencyType,
                                         calendarType,
                                         busDayAdjustType,
                                         dateGenRuleType).generate()

        for swaptionDt in swaptionFloatDates:
            foundDt = False
            for gridDt in self._gridDates;
                if swaptionDt == gridDt:
                    foundDt = True
                    break
            if foundDt is False:
                raise FinError("Swaption float leg not on grid.")

        swaptionFixedDates = FinSchedule(settlementDate,
                                         maturityDate,
                                         fixedFrequencyType,
                                         calendarType,
                                         busDayAdjustType,
                                         dateGenRuleType).generate()

        for swaptionDt in swaptionFixedDates:
            foundDt = False
            for gridDt in self._gridDates;
                if swaptionDt == gridDt:
                    foundDt = True
                    break
            if foundDt is False:
                raise FinError("Swaption fixed leg not on grid.")
        
        a = 0
        b = 0
        
        for gridDt in self._gridDates;
            if gridDt == exerciseDate:
                break
            else:
                a += 1

        for gridDt in self._gridDates;
            if gridDt == maturityDate:
                break
            else:
                b += 1

         v = LMMSwaptionPricer(strike, a, b, numPaths, 
                               fwd0, fwds, taus, isPayer):

        return v    

###############################################################################


     def valueCapFloor(self, 
                       exerciseDate
                       maturityDate, 
                       couponStrike):
###############################################################################

    def __repr__(self):
        ''' Function to allow us to print the swaption details. '''

        s = labelToString("SETTLEMENT DATE", self._settlementDate)
        s += labelToString("EXERCISE DATE", self._exerciseDate)
        s += labelToString("SWAPTION TYPE", str(self._swaptionType))
        s += labelToString("MATURITY DATE", self._maturityDate)
        s += labelToString("SWAP NOTIONAL", self._notional)
        s += labelToString("FIXED COUPON", self._fixedCoupon * 100)
        s += labelToString("FIXED FREQUENCY", str(self._fixedFrequencyType))
        s += labelToString("FIXED DAY COUNT", str(self._fixedDayCountType))
        s += labelToString("FLOAT FREQUENCY", str(self._floatFrequencyType))
        s += labelToString("FLOAT DAY COUNT", str(self._floatDayCountType))

        if self._pv01 is not None:
            s += labelToString("PV01", self._pv01)
            s += labelToString("FWD SWAP RATE", self._fwdSwapRate*100)
            s += labelToString("FWD DF TO EXPIRY", self._forwardDf, "")

        return s

###############################################################################

    def print(self):
        ''' Alternative print method. '''

        print(self)

###############################################################################
