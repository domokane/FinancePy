##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...finutils.FinGlobalVariables import gDaysInYear
from ...models.FinModelRatesHW import FinModelRatesHW
from ...models.FinModelRatesBK import FinModelRatesBK
from ...finutils.FinError import FinError
from ...finutils.FinDate import FinDate
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes
from ...market.curves.FinDiscountCurve import FinDiscountCurve
from .FinBond import FinBond
from ...finutils.FinOptionTypes import FinOptionTypes, FinOptionExerciseTypes
from enum import Enum
import numpy as np


###############################################################################


class FinBondModelTypes(Enum):
    BLACK = 1
    HO_LEE = 2
    HULL_WHITE = 3
    BLACK_KARASINSKI = 4

###############################################################################


class FinBondOption():
    ''' Class for options on fixed coupon bonds. '''

    def __init__(self,
                 bond: FinBond,
                 expiryDate: FinDate,
                 strikePrice: float,
                 face: float,
                 optionType: FinOptionTypes):

        checkArgumentTypes(self.__init__, locals())

        self._expiryDate = expiryDate
        self._strikePrice = strikePrice
        self._bond = bond
        self._optionType = optionType
        self._face = face

###############################################################################

    def value(self,
              valueDate: FinDate,
              discountCurve: FinDiscountCurve,
              model):
        ''' Value a bond option (option on a bond) using the specified model
        which include Hull-White Tree, Black-Karasinski Tree. '''

        texp = (self._expiryDate - valueDate) / gDaysInYear

        dfTimes = discountCurve._times
        dfValues = discountCurve._dfValues

        # We need all of the flows in case the option is American
        self._bond._calculateFlowDates(valueDate)
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

            if self._optionType == FinOptionTypes.EUROPEAN_CALL \
                    and model._useJamshidian is True:

                v = model.europeanBondOption_Jamshidian(texp,
                                                        self._strikePrice,
                                                        self._face,
                                                        cpnTimes,
                                                        cpnAmounts,
                                                        dfTimes, dfValues)

                return v['call']

            elif self._optionType == FinOptionTypes.EUROPEAN_PUT  \
                    and model._useJamshidian is True:

                v = model.europeanBondOption_Jamshidian(texp,
                                                        self._strikePrice,
                                                        self._face,
                                                        cpnTimes,
                                                        cpnAmounts,
                                                        dfTimes, dfValues)

                return v['put']

            elif self._optionType == FinOptionTypes.EUROPEAN_CALL  \
                    and model._useJamshidian is False:

                model.buildTree(texp, dfTimes, dfValues)
                exerciseType = FinOptionTypes.EUROPEAN_CALL
                v1 = model.americanBondOption_Tree(texp,
                                                   self._strikePrice,
                                                   self._face,
                                                   cpnTimes,
                                                   cpnAmounts,
                                                   exerciseType)

                model._numTimeSteps += 1
                model.buildTree(texp, dfTimes, dfValues)
                v2 = model.americanBondOption_Tree(texp,
                                                   self._strikePrice,
                                                   self._face,
                                                   cpnTimes,
                                                   cpnAmounts,
                                                   exerciseType)

                v = (v1['call'] + v2['call'])/2.0
                return v

            elif self._optionType == FinOptionTypes.EUROPEAN_PUT  \
                    and model._useJamshidian is False:

                exerciseType = FinOptionExerciseTypes.AMERICAN
                model.buildTree(texp, dfTimes, dfValues)
                v1 = model.americanBondOption_Tree(texp, self._strikePrice,
                                                   self._face,
                                                   cpnTimes, cpnAmounts,
                                                   exerciseType)

                model._numTimeSteps += 1
                model.buildTree(texp, dfTimes, dfValues)
                v2 = model.americanBondOption_Tree(texp, self._strikePrice,
                                                   self._face,
                                                   cpnTimes, cpnAmounts,
                                                   exerciseType)

                model._numTimeSteps -= 1
                v = (v1['put'] + v2['put'])/2.0
                return v

            elif self._optionType == FinOptionTypes.AMERICAN_CALL:

                exerciseType = FinOptionExerciseTypes.AMERICAN

                model.buildTree(texp, dfTimes, dfValues)
                v1 = model.americanBondOption_Tree(texp, self._strikePrice,
                                                   self._face,
                                                   cpnTimes, cpnAmounts,
                                                   exerciseType)

                model._numTimeSteps += 1
                model.buildTree(texp, dfTimes, dfValues)
                v2 = model.americanBondOption_Tree(texp, self._strikePrice,
                                                   self._face,
                                                   cpnTimes, cpnAmounts,
                                                   exerciseType)

                model._numTimeSteps -= 1
                v = (v1['call'] + v2['call'])/2.0
                return v

            elif self._optionType == FinOptionTypes.AMERICAN_PUT:

                exerciseType = FinOptionExerciseTypes.AMERICAN
                model.buildTree(texp, dfTimes, dfValues)
                v1 = model.americanBondOption_Tree(texp, self._strikePrice,
                                                   self._face,
                                                   cpnTimes, cpnAmounts,
                                                   exerciseType)

                model._numTimeSteps += 1
                model.buildTree(texp, dfTimes, dfValues)
                v2 = model.americanBondOption_Tree(texp, self._strikePrice,
                                                   self._face,
                                                   cpnTimes, cpnAmounts,
                                                   exerciseType)

                model._numTimeSteps -= 1
                v = (v1['put'] + v2['put'])/2.0
                return v

        elif type(model) == FinModelRatesBK:

            maturityDate = self._bond._maturityDate
            tmat = (maturityDate - valueDate) / gDaysInYear

            if self._optionType == FinOptionTypes.EUROPEAN_CALL:

                exerciseType = FinOptionExerciseTypes.EUROPEAN
                model.buildTree(tmat, dfTimes, dfValues)
                v1 = model.bondOption(texp, self._strikePrice, self._face,
                                      cpnTimes, cpnAmounts, exerciseType)
                model._numTimeSteps += 1
                model.buildTree(tmat, dfTimes, dfValues)
                v2 = model.bondOption(texp, self._strikePrice, self._face,
                                      cpnTimes, cpnAmounts, exerciseType)

                v = (v1['call'] + v2['call'])/2.0
                return v

            elif self._optionType == FinOptionTypes.EUROPEAN_PUT:

                exerciseType = FinOptionExerciseTypes.EUROPEAN
                model.buildTree(tmat, dfTimes, dfValues)
                v1 = model.bondOption(texp, self._strikePrice, self._face,
                                      cpnTimes, cpnAmounts, exerciseType)

                model._numTimeSteps += 1
                model.buildTree(tmat, dfTimes, dfValues)
                v2 = model.bondOption(texp, self._strikePrice, self._face,
                                      cpnTimes, cpnAmounts, exerciseType)

                v = (v1['put'] + v2['put'])/2.0
                return v

            elif self._optionType == FinOptionTypes.AMERICAN_CALL:

                exerciseType = FinOptionExerciseTypes.AMERICAN

                model.buildTree(tmat, dfTimes, dfValues)
                v1 = model.bondOption(texp, self._strikePrice, self._face,
                                      cpnTimes, cpnAmounts, exerciseType)

                model._numTimeSteps += 1
                model.buildTree(tmat, dfTimes, dfValues)
                v2 = model.bondOption(texp, self._strikePrice, self._face,
                                      cpnTimes, cpnAmounts, exerciseType)

                v = (v1['call'] + v2['call'])/2.0
                return v

            elif self._optionType == FinOptionTypes.AMERICAN_PUT:

                exerciseType = FinOptionExerciseTypes.AMERICAN

                model.buildTree(tmat, dfTimes, dfValues)
                v1 = model.bondOption(texp, self._strikePrice, self._face,
                                      cpnTimes, cpnAmounts, exerciseType)
                model._numTimeSteps += 1
                model.buildTree(tmat, dfTimes, dfValues)
                v2 = model.bondOption(texp, self._strikePrice, self._face,
                                      cpnTimes, cpnAmounts, exerciseType)

                v = (v1['put'] + v2['put'])/2.0
                return v

        else:
            raise FinError("Unknown model and option combination")

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("EXPIRY DATE", self._expiryDate)
        s += labelToString("STRIKE", self._strikePrice)
        s += labelToString("OPTION TYPE", self._optionType)
        s += labelToString("FACE", self._face, "")
        s += "Underlying Bond\n"
        s += str(self._bond)
        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
