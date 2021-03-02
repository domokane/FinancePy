##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...utils.FinGlobalVariables import gDaysInYear
from ...utils.FinError import FinError
from ...utils.Date import Date
from ...utils.FinHelperFunctions import labelToString, checkArgumentTypes
from ...market.curves.FinDiscountCurve import FinDiscountCurve
from ...utils.FinGlobalTypes import FinOptionTypes, FinExerciseTypes
from ...products.bonds.FinBond import FinBond

from enum import Enum
import numpy as np

###############################################################################
# TODO: Add BDT model to valuation
# TODO: Reorganise code - too much duplication
###############################################################################


class FinBondModelTypes(Enum):
    BLACK = 1
    HO_LEE = 2
    HULL_WHITE = 3
    BLACK_KARASINSKI = 4

###############################################################################


class FinBondOption():
    """ Class for options on fixed coupon bonds. These are options to either
    buy or sell a bond on or before a specific future expiry date at a strike
    price that is set on trade date. A European option only allows the bond to
    be exercised into on a specific expiry date. An American option allows the
    option holder to exercise early, potentially allowing earlier coupons to
    be received. """

    def __init__(self,
                 bond: FinBond,
                 expiry_date: Date,
                 strikePrice: float,
                 face_amount: float,
                 optionType: FinOptionTypes):

        checkArgumentTypes(self.__init__, locals())

        self._expiry_date = expiry_date
        self._strikePrice = strikePrice
        self._bond = bond
        self._optionType = optionType
        self._face_amount = face_amount

###############################################################################

    def value(self,
              valuation_date: Date,
              discount_curve: FinDiscountCurve,
              model):
        """ Value a bond option (option on a bond) using a specified model
        which include the Hull-White, Black-Karasinski and Black-Derman-Toy
        model which are all implemented as short rate tree models. """

        texp = (self._expiry_date - valuation_date) / gDaysInYear
        tmat = (self._bond._maturity_date - valuation_date) / gDaysInYear

        dfTimes = discount_curve._times
        dfValues = discount_curve._dfs

        # We need all of the flows in case the option is American and some occur before expiry
        flow_dates = self._bond._flow_dates
        flow_amounts = self._bond._flow_amounts

        coupon_times = []
        coupon_flows = []

        numFlows = len(self._bond._flow_dates)

        # Want the first flow to be the previous coupon date
        # This is needed to calculate accrued correctly        
        for i in range(1, numFlows):
            pcd = flow_dates[i-1]
            ncd = flow_dates[i]
            if pcd < valuation_date and ncd > valuation_date:
                flow_time = (pcd - valuation_date) / gDaysInYear
                coupon_times.append(flow_time)
                coupon_flows.append(flow_amounts[i])
                break

        for i in range(1, numFlows):
            if flow_dates[i] == valuation_date:
                coupon_times.append(0.0)
                coupon_flows.append(flow_amounts[i])
                
        # Now calculate the remaining coupons
        for i in range(1, numFlows):
            ncd = flow_dates[i]
            if ncd > valuation_date:
                flow_time = (ncd - valuation_date) / gDaysInYear
                coupon_times.append(flow_time)
                coupon_flows.append(flow_amounts[i])

        ##################################################################

        coupon_times = np.array(coupon_times)
        coupon_flows = np.array(coupon_flows)

        exerciseType = FinExerciseTypes.AMERICAN

        if self._optionType == FinOptionTypes.EUROPEAN_CALL \
            or self._optionType == FinOptionTypes.EUROPEAN_PUT:
                exerciseType = FinExerciseTypes.EUROPEAN                

        # This is wasteful if the model is Jamshidian but how to do neat design ?
        model.buildTree(tmat, dfTimes, dfValues)

        v = model.bondOption(texp, self._strikePrice, self._face_amount,
                             coupon_times, coupon_flows, exerciseType)

        if self._optionType == FinOptionTypes.EUROPEAN_CALL \
            or self._optionType == FinOptionTypes.AMERICAN_CALL:
                return v['call']    
        elif self._optionType == FinOptionTypes.EUROPEAN_PUT \
            or self._optionType == FinOptionTypes.AMERICAN_PUT:
                return v['put']
        else:
            print(self._optionType)
            raise FinError("Unknown option type.")

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("EXPIRY DATE", self._expiry_date)
        s += labelToString("STRIKE", self._strikePrice)
        s += labelToString("OPTION TYPE", self._optionType)
        s += labelToString("FACE AMOUNT", self._face_amount, "")
        s += "Underlying Bond\n"
        s += str(self._bond)
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
