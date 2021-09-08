##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...utils.global_vars import gDaysInYear
from ...utils.error import FinError
from ...utils.date import Date
from ...utils.helpers import label_to_string, check_argument_types
from ...market.curves.discount_curve import DiscountCurve
from ...utils.global_types import OptionTypes, FinExerciseTypes
from ...products.bonds.bond import Bond

from enum import Enum
import numpy as np

###############################################################################
# TODO: Add BDT model to valuation
# TODO: Reorganise code - too much duplication
###############################################################################


class BondModelTypes(Enum):
    BLACK = 1
    HO_LEE = 2
    HULL_WHITE = 3
    BLACK_KARASINSKI = 4

###############################################################################


class BondOption():
    """ Class for options on fixed coupon bonds. These are options to either
    buy or sell a bond on or before a specific future expiry date at a strike
    price that is set on trade date. A European option only allows the bond to
    be exercised into on a specific expiry date. An American option allows the
    option holder to exercise early, potentially allowing earlier coupons to
    be received. """

    def __init__(self,
                 bond: Bond,
                 expiry_date: Date,
                 strike_price: float,
                 face_amount: float,
                 option_type: OptionTypes):

        check_argument_types(self.__init__, locals())

        self._expiry_date = expiry_date
        self._strike_price = strike_price
        self._bond = bond
        self._option_type = option_type
        self._face_amount = face_amount

###############################################################################

    def value(self,
              valuation_date: Date,
              discount_curve: DiscountCurve,
              model):
        """ Value a bond option (option on a bond) using a specified model
        which include the Hull-White, Black-Karasinski and Black-Derman-Toy
        model which are all implemented as short rate tree models. """

        texp = (self._expiry_date - valuation_date) / gDaysInYear
        tmat = (self._bond._maturity_date - valuation_date) / gDaysInYear

        df_times = discount_curve._times
        df_values = discount_curve._dfs

        # We need all of the flows in case the option is American
        # and some occur before expiry
        flow_dates = self._bond._flow_dates
        flow_amounts = self._bond._flow_amounts

        coupon_times = []
        coupon_flows = []

        num_flows = len(self._bond._flow_dates)

        # Want the first flow to be the previous coupon date
        # This is needed to calculate accrued correctly
        for i in range(1, num_flows):
            pcd = flow_dates[i-1]
            ncd = flow_dates[i]
            if pcd < valuation_date and ncd > valuation_date:
                flow_time = (pcd - valuation_date) / gDaysInYear
                coupon_times.append(flow_time)
                coupon_flows.append(flow_amounts[i])
                break

        for i in range(1, num_flows):
            if flow_dates[i] == valuation_date:
                coupon_times.append(0.0)
                coupon_flows.append(flow_amounts[i])

        # Now calculate the remaining coupons
        for i in range(1, num_flows):
            ncd = flow_dates[i]
            if ncd > valuation_date:
                flow_time = (ncd - valuation_date) / gDaysInYear
                coupon_times.append(flow_time)
                coupon_flows.append(flow_amounts[i])

        ##################################################################

        coupon_times = np.array(coupon_times)
        coupon_flows = np.array(coupon_flows)

        exercise_type = FinExerciseTypes.AMERICAN

        if self._option_type == OptionTypes.EUROPEAN_CALL \
                or self._option_type == OptionTypes.EUROPEAN_PUT:
            exercise_type = FinExerciseTypes.EUROPEAN

        # This is wasteful if model is Jamshidian but how to do neat design
        model.build_tree(tmat, df_times, df_values)

        v = model.bond_option(texp, self._strike_price, self._face_amount,
                              coupon_times, coupon_flows, exercise_type)

        if self._option_type == OptionTypes.EUROPEAN_CALL \
                or self._option_type == OptionTypes.AMERICAN_CALL:
            return v['call']
        elif self._option_type == OptionTypes.EUROPEAN_PUT \
                or self._option_type == OptionTypes.AMERICAN_PUT:
            return v['put']
        else:
            print(self._option_type)
            raise FinError("Unknown option type.")

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self._expiry_date)
        s += label_to_string("STRIKE", self._strike_price)
        s += label_to_string("OPTION TYPE", self._option_type)
        s += label_to_string("FACE AMOUNT", self._face_amount, "")
        s += "Underlying Bond\n"
        s += str(self._bond)
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
