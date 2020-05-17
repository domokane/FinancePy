# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

from ...finutils.FinGlobalVariables import gDaysInYear
from ...models.FinModelRatesHW import FinModelRatesHW
from ...models.FinModelRatesBK import FinModelRatesBK
from ...finutils.FinError import FinError

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
        cpnTimes = [0.0]
        cpnAmounts = [0.0]

        # The first flow is always the previous coupon date
        for flowDate in self._bond._flowDates[1:]:
            cpnTime = (flowDate - valueDate) / gDaysInYear
            cpnTimes.append(cpnTime)
            cpnAmounts.append(cpn)

        cpnTimes = np.array(cpnTimes)
        cpnAmounts = np.array(cpnAmounts)

        if np.any(cpnTimes < 0.0):
            raise FinError("No coupon times can be before the value date.")

        if isinstance(model, FinModelRatesHW):

            if self._optionType == FinBondOptionTypes.EUROPEAN_CALL \
                    and model._useJamshidian is True:

                v = model.europeanBondOption_Jamshidian(texp,
                                                        self._strikePrice,
                                                        self._face,
                                                        cpnTimes,
                                                        cpnAmounts,
                                                        dfTimes, dfValues)

                return v['call']

            elif self._optionType == FinBondOptionTypes.EUROPEAN_PUT  \
                    and model._useJamshidian is True:

                v = model.europeanBondOption_Jamshidian(texp,
                                                        self._strikePrice,
                                                        self._face,
                                                        cpnTimes,
                                                        cpnAmounts,
                                                        dfTimes, dfValues)

                return v['put']

            elif self._optionType == FinBondOptionTypes.EUROPEAN_CALL  \
                    and model._useJamshidian is False:

                model.buildTree(texp, dfTimes, dfValues)
                americanExercise = False

                v1 = model.americanBondOption_Tree(texp, self._strikePrice,
                                                   self._face,
                                                   cpnTimes, cpnAmounts,
                                                   americanExercise)

                model._numTimeSteps += 1
                model.buildTree(texp, dfTimes, dfValues)
                v2 = model.americanBondOption_Tree(texp, self._strikePrice,
                                                   self._face,
                                                   cpnTimes, cpnAmounts,
                                                   americanExercise)

                v = (v1['call'] + v2['call'])/2.0
                return v

            elif self._optionType == FinBondOptionTypes.EUROPEAN_PUT  \
                    and model._useJamshidian is False:

                americanExercise = False
                model.buildTree(texp, dfTimes, dfValues)
                v1 = model.americanBondOption_Tree(texp, self._strikePrice,
                                                   self._face,
                                                   cpnTimes, cpnAmounts,
                                                   americanExercise)

                model._numTimeSteps += 1
                model.buildTree(texp, dfTimes, dfValues)
                v2 = model.americanBondOption_Tree(texp, self._strikePrice,
                                                   self._face,
                                                   cpnTimes, cpnAmounts,
                                                   americanExercise)

                model._numTimeSteps -= 1
                v = (v1['put'] + v2['put'])/2.0
                return v

            elif self._optionType == FinBondOptionTypes.AMERICAN_CALL:

                americanExercise = True
                model.buildTree(texp, dfTimes, dfValues)
                v1 = model.americanBondOption_Tree(texp, self._strikePrice,
                                                   self._face,
                                                   cpnTimes, cpnAmounts,
                                                   americanExercise)

                model._numTimeSteps += 1
                model.buildTree(texp, dfTimes, dfValues)
                v2 = model.americanBondOption_Tree(texp, self._strikePrice,
                                                   self._face,
                                                   cpnTimes, cpnAmounts,
                                                   americanExercise)

                model._numTimeSteps -= 1
                v = (v1['call'] + v2['call'])/2.0
                return v

            elif self._optionType == FinBondOptionTypes.AMERICAN_PUT:

                americanExercise = True
                model.buildTree(texp, dfTimes, dfValues)
                v1 = model.americanBondOption_Tree(texp, self._strikePrice,
                                                   self._face,
                                                   cpnTimes, cpnAmounts,
                                                   americanExercise)

                model._numTimeSteps += 1
                model.buildTree(texp, dfTimes, dfValues)
                v2 = model.americanBondOption_Tree(texp, self._strikePrice,
                                                   self._face,
                                                   cpnTimes, cpnAmounts,
                                                   americanExercise)

                model._numTimeSteps -= 1
                v = (v1['put'] + v2['put'])/2.0
                return v

        elif type(model) == FinModelRatesBK:

            maturityDate = self._bond._maturityDate
            tmat = (maturityDate - valueDate) / gDaysInYear

            if self._optionType == FinBondOptionTypes.EUROPEAN_CALL:

                model.buildTree(tmat, dfTimes, dfValues)
                v1 = model.bondOption(texp, self._strikePrice, self._face,
                                      cpnTimes, cpnAmounts, False)
                model._numTimeSteps += 1
                model.buildTree(tmat, dfTimes, dfValues)
                v2 = model.bondOption(texp, self._strikePrice, self._face,
                                      cpnTimes, cpnAmounts, True)

                v = (v1['call'] + v2['call'])/2.0
                return v

            elif self._optionType == FinBondOptionTypes.EUROPEAN_PUT:

                model.buildTree(tmat, dfTimes, dfValues)
                v1 = model.bondOption(texp, self._strikePrice, self._face,
                                      cpnTimes, cpnAmounts, False)

                model._numTimeSteps += 1
                model.buildTree(tmat, dfTimes, dfValues)
                v2 = model.bondOption(texp, self._strikePrice, self._face,
                                      cpnTimes, cpnAmounts, True)

                v = (v1['put'] + v2['put'])/2.0
                return v

            elif self._optionType == FinBondOptionTypes.AMERICAN_CALL:

                model.buildTree(tmat, dfTimes, dfValues)
                v1 = model.bondOption(texp, self._strikePrice, self._face,
                                      cpnTimes, cpnAmounts, True)

                model._numTimeSteps += 1
                model.buildTree(tmat, dfTimes, dfValues)
                v2 = model.bondOption(texp, self._strikePrice, self._face,
                                      cpnTimes, cpnAmounts, True)

                v = (v1['call'] + v2['call'])/2.0
                return v

            elif self._optionType == FinBondOptionTypes.AMERICAN_PUT:

                model.buildTree(tmat, dfTimes, dfValues)
                v1 = model.bondOption(texp, self._strikePrice, self._face,
                                      cpnTimes, cpnAmounts, True)
                model._numTimeSteps += 1
                model.buildTree(tmat, dfTimes, dfValues)
                v2 = model.bondOption(texp, self._strikePrice, self._face,
                                      cpnTimes, cpnAmounts, True)

                v = (v1['put'] + v2['put'])/2.0
                return v

        else:
            raise FinError("Unknown model and option combination")

        return 999.0

###############################################################################
