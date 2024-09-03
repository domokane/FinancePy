##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...utils.global_vars import g_days_in_year
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
                 expiry_dt: Date,
                 strike_price: float,
                 option_type: OptionTypes):

        check_argument_types(self.__init__, locals())

        self.expiry_dt = expiry_dt
        self.strike_price = strike_price
        self.bond = bond
        self.option_type = option_type
        self.par = 100.0

###############################################################################

    def value(self,
              value_dt: Date,
              discount_curve: DiscountCurve,
              model):
        """ Value a bond option (option on a bond) using a specified model
        which include the Hull-White, Black-Karasinski and Black-Derman-Toy
        model which are all implemented as short rate tree models. """

        t_exp = (self.expiry_dt - value_dt) / g_days_in_year
        t_mat = (self.bond.maturity_dt - value_dt) / g_days_in_year

        df_times = discount_curve._times
        df_values = discount_curve._dfs

        # We need all the flows in case the option is American
        # and some occur before expiry
        flow_dts = self.bond.cpn_dts
        flow_amounts = self.bond.flow_amounts

        cpn_times = []
        cpn_flows = []

        num_flows = len(self.bond.cpn_dts)

        # Want the first flow to be the previous coupon date
        # This is needed to calculate accrued correctly
        for i in range(1, num_flows):
            pcd = flow_dts[i-1]
            ncd = flow_dts[i]
            if pcd < value_dt and ncd > value_dt:
                flow_time = (pcd - value_dt) / g_days_in_year
                cpn_times.append(flow_time)
                cpn_flows.append(flow_amounts[i])
                break

        for i in range(1, num_flows):
            if flow_dts[i] == value_dt:
                cpn_times.append(0.0)
                cpn_flows.append(flow_amounts[i])

        # Now calculate the remaining coupons
        for i in range(1, num_flows):
            ncd = flow_dts[i]
            if ncd > value_dt:
                flow_time = (ncd - value_dt) / g_days_in_year
                cpn_times.append(flow_time)
                cpn_flows.append(flow_amounts[i])

        ##################################################################

        cpn_times = np.array(cpn_times)
        cpn_flows = np.array(cpn_flows)

        exercise_type = FinExerciseTypes.AMERICAN

        if self.option_type == OptionTypes.EUROPEAN_CALL \
                or self.option_type == OptionTypes.EUROPEAN_PUT:
            exercise_type = FinExerciseTypes.EUROPEAN

        # This is wasteful if model is Jamshidian but how to do neat design
        model.build_tree(t_mat, df_times, df_values)

        v = model.bond_option(t_exp, self.strike_price, self.par,
                              cpn_times, cpn_flows, exercise_type)

        if self.option_type == OptionTypes.EUROPEAN_CALL \
                or self.option_type == OptionTypes.AMERICAN_CALL:
            return v['call']
        elif self.option_type == OptionTypes.EUROPEAN_PUT \
                or self.option_type == OptionTypes.AMERICAN_PUT:
            return v['put']
        else:
            print(self.option_type)
            raise FinError("Unknown option type.")

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("STRIKE", self.strike_price)
        s += label_to_string("OPTION TYPE", self.option_type)
        s += "Underlying Bond\n"
        s += str(self.bond)
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
