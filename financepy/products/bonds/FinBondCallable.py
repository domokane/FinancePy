# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

from ...finutils.FinGlobalVariables import gDaysInYear
from ...models.FinHullWhiteRateModel import FinHullWhiteRateModel
from ...models.FinBlackKarasinskiRateModel import FinBlackKarasinskiRateModel

from enum import Enum
import numpy as np

###############################################################################


class FinBondModelTypes(Enum):
    BLACK = 1
    HO_LEE = 2
    HULL_WHITE = 3
    BLACK_KARASINSKI = 4

###############################################################################


class FinBondOptionTypes(Enum):
    EUROPEAN_CALL = 1
    EUROPEAN_PUT = 2
    AMERICAN_CALL = 3
    AMERICAN_PUT = 4


###############################################################################


class FinBondCallable():
    ''' Class for options on fixed coupon bonds. '''

    def __init__(self,
                 bond,
                 callDates,
                 callPrices,
                 putDates,
                 putPrices,
                 face):

        self._callDates = callDates
        self._callPrices = callPrices
        self._putDates = putDates
        self._putPrices = putPrices
        self._bond = bond
        self._face = face

###############################################################################

    def value(self,
              valueDate,
              discountCurve,
              model):
        ''' Value the bond option using the specified model. '''

        texp = (self._expiryDate - valueDate) / gDaysInYear

        dfTimes = discountCurve._times
        dfValues = discountCurve._values

        # We need all of the flows in case the option is American
        self._bond.calculateFlowDates(valueDate)
        cpn = self._bond._coupon/self._bond._frequency
        cpnTimes = []
        cpnAmounts = []

        for flowDate in self._bond._flowDates:
            cpnTime = (flowDate - valueDate) / gDaysInYear
            cpnTimes.append(cpnTime)
            cpnAmounts.append(cpn)

        cpnTimes = np.array(cpnTimes)
        cpnAmounts = np.array(cpnAmounts)

        callTimes = []
        for dt in self._callDates:
            callTime = (dt - valueDate) / gDaysInYear
            callTimes.append(callTime)
        callTimes = np.array(callTimes)
        callPrices = np.array(self._callPrices)

        putTimes = []
        for dt in self._putDates:
            putTime = (dt - valueDate) / gDaysInYear
            putTimes.append(putTime)
        putTimes = np.array(putTimes)
        putPrices = np.array(self._putPrices)

        if type(model) == FinHullWhiteRateModel:

            numTimeSteps = 100
            model.buildTree(texp, numTimeSteps, dfTimes, dfValues)

            v = model.callablePuttableBond(cpnTimes, cpnAmounts,
                                           callTimes, callPrices,
                                           putTimes, putPrices)

            return v[0]

        elif type(model) == FinBlackKarasinskiRateModel:

            maturityDate = self._bond._maturityDate
            numTimeSteps = 100
            tmat = (maturityDate - valueDate) / gDaysInYear
            model.buildTree(tmat, numTimeSteps, dfTimes, dfValues)

            v = model.callablePuttableBond(cpnTimes, cpnAmounts,
                                           callTimes, callPrices,
                                           putTimes, putPrices)

            return v


        return 999.0

###############################################################################
