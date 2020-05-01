# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

from ...finutils.FinGlobalVariables import gDaysInYear
from ...models.FinModelRatesHullWhite import FinModelRatesHullWhite
from ...models.FinModelRatesBlackKarasinski import FinModelRatesBlackKarasinski

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


class FinBondOption():
    ''' Class for options on fixed coupon bonds. '''

    def __init__(self,
                 bond,
                 expiryDate,
                 strikePrice,
                 face,
                 optionType):

        self._expiryDate = expiryDate
        self._strikePrice = strikePrice
        self._bond = bond
        self._optionType = optionType
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

        if type(model) == FinModelRatesHullWhite:

            if self._optionType == FinBondOptionTypes.EUROPEAN_CALL:

                v = model.bondOption_Jamshidian(texp, self._strikePrice,
                                                cpnTimes, cpnAmounts,
                                                dfTimes, dfValues)

                return v[0]

            elif self._optionType == FinBondOptionTypes.EUROPEAN_PUT:

                v = model.bondOption_Jamshidian(texp, self._strikePrice,
                                                cpnTimes, cpnAmounts,
                                                dfTimes, dfValues)

                return v[1]

            elif self._optionType == FinBondOptionTypes.AMERICAN_CALL:

                numTimeSteps = 100
                model.buildTree(texp, numTimeSteps, dfTimes, dfValues)
                americanExercise = True

                v = model.americanBondOption_Tree(texp, self._strikePrice,
                                                  self._face,
                                                  cpnTimes, cpnAmounts,
                                                  americanExercise)

                return v[0]

            elif self._optionType == FinBondOptionTypes.AMERICAN_PUT:

                numTimeSteps = 100
                model.buildTree(texp, numTimeSteps, dfTimes, dfValues)
                americanExercise = True

                v = model.americanBondOption_Tree(texp, self._strikePrice,
                                                  self._face,
                                                  cpnTimes, cpnAmounts,
                                                  americanExercise)

                return v[1]

        elif type(model) == FinModelRatesBlackKarasinski:

            maturityDate = self._bond._maturityDate
            numTimeSteps = 100
            tmat = (maturityDate - valueDate) / gDaysInYear
            model.buildTree(tmat, numTimeSteps, dfTimes, dfValues)

            if self._optionType == FinBondOptionTypes.EUROPEAN_CALL:

                v = model.bondOption(texp, self._strikePrice, self._face,
                                     cpnTimes, cpnAmounts, False)

                return v[0]

            elif self._optionType == FinBondOptionTypes.EUROPEAN_PUT:

                v = model.bondOption(texp, self._strikePrice, self._face,
                                     cpnTimes, cpnAmounts, False)

                return v[1]

            elif self._optionType == FinBondOptionTypes.AMERICAN_CALL:

                v = model.bondOption(texp, self._strikePrice, self._face,
                                     cpnTimes, cpnAmounts, True)

                return v[0]

            elif self._optionType == FinBondOptionTypes.AMERICAN_PUT:

                v = model.bondOption(texp, self._strikePrice, self._face,
                                     cpnTimes, cpnAmounts, True)

                return v[1]

        return 999.0

###############################################################################
