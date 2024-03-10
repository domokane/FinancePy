##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from enum import Enum

from typing import List
import numpy as np

from ...utils.global_vars import gDaysInYear
from ...models.hw_tree import HWTree
from ...models.bk_tree import BKTree
from ...utils.error import FinError
from ...utils.frequency import FrequencyTypes
from ...utils.day_count import DayCountTypes
from ...products.bonds.bond import Bond

from ...utils.date import Date
from ...utils.helpers import label_to_string, check_argument_types
from ...market.curves.discount_curve import DiscountCurve


###############################################################################
# TODO: Make it possible to specify start and end of American Callable/Puttable
###############################################################################


class BondModelTypes(Enum):
    BLACK = 1
    HO_LEE = 2
    HULL_WHITE = 3
    BLACK_KARASINSKI = 4

###############################################################################


class BondOptionTypes(Enum):
    EUROPEAN_CALL = 1
    EUROPEAN_PUT = 2
    AMERICAN_CALL = 3
    AMERICAN_PUT = 4


###############################################################################


class BondEmbeddedOption:
    """ Class for fixed coupon bonds with embedded call or put optionality. """

    def __init__(self,
                 issue_dt: Date,
                 maturity_dt: Date,  # Date
                 coupon: float,  # Annualised coupon - 0.03 = 3.00%
                 freq_type: FrequencyTypes,
                 dc_type: DayCountTypes,
                 call_dts: List[Date],
                 call_prices: List[float],
                 put_dts: List[Date],
                 put_prices: List[float]):
        """ Create a BondEmbeddedOption object with a maturity date, coupon
        and all the bond inputs. """

        check_argument_types(self.__init__, locals())

        self._issue_dt = issue_dt
        self._maturity_dt = maturity_dt
        self._cpn = coupon
        self._freq_type = freq_type
        self._dc_type = dc_type

        self._ex_div_days = 0

        self._bond = Bond(issue_dt,
                          maturity_dt,
                          coupon,
                          freq_type,
                          dc_type,
                          self._ex_div_days)

        # Validate call and put schedules
        for dt in call_dts:
            if dt > self._maturity_dt:
                raise FinError("Call date after bond maturity date")

        if len(call_dts) > 0:
            dtprev = call_dts[0]
            for dt in call_dts[1:]:
                if dt <= dtprev:
                    raise FinError("Call dates not increasing")
                else:
                    dtprev = dt

        for dt in put_dts:
            if dt > self._maturity_dt:
                raise FinError("Put date after bond maturity date")

        if len(put_dts) > 0:
            dtprev = put_dts[0]
            for dt in put_dts[1:]:
                if dt <= dtprev:
                    raise FinError("Put dates not increasing")
                else:
                    dtprev = dt

        for px in call_prices:
            if px < 0.0:
                raise FinError("Call price must be positive.")

        for px in put_prices:
            if px < 0.0:
                raise FinError("Put price must be positive.")

        if len(call_dts) != len(call_prices):
            raise FinError("Number of call dates and call prices not the same")

        if len(put_dts) != len(put_prices):
            raise FinError("Number of put dates and put prices not the same")

        self._call_dts = call_dts
        self._call_prices = call_prices
        self._put_dts = put_dts
        self._put_prices = put_prices
        self._par = 100.0
        self._bond._calculate_cpn_dts()

###############################################################################

    def value(self,
              settle_dt: Date,
              discount_curve: DiscountCurve,
              model):
        """ Value the bond that settles on the specified date that can have
        both embedded call and put options. This is done using the specified
        model and a discount curve. """

        # Generate bond coupon flow schedule
        cpn = self._bond.cpn/self._bond.freq

        cpn_times = []
        cpn_amounts = []

        for flow_dt in self._bond.cpn_dts[1:]:
            if flow_dt > settle_dt:
                cpn_time = (flow_dt - settle_dt) / gDaysInYear
                cpn_times.append(cpn_time)
                cpn_amounts.append(cpn)

        cpn_times = np.array(cpn_times)
        cpn_amounts = np.array(cpn_amounts)

        # Generate bond call times and prices
        call_times = []
        for dt in self._call_dts:
            if dt > settle_dt:
                call_time = (dt - settle_dt) / gDaysInYear
                call_times.append(call_time)
        call_times = np.array(call_times)
        call_prices = np.array(self._call_prices)

        # Generate bond put times and prices
        put_times = []
        for dt in self._put_dts:
            if dt > settle_dt:
                put_time = (dt - settle_dt) / gDaysInYear
                put_times.append(put_time)
        put_times = np.array(put_times)
        put_prices = np.array(self._put_prices)

        maturity_dt = self._bond.maturity_dt
        t_mat = (maturity_dt - settle_dt) / gDaysInYear
        df_times = discount_curve._times
        df_values = discount_curve._dfs

        face_amount = self._par

        if isinstance(model, HWTree):

            """ We need to build the tree out to the bond maturity date. To be
            more precise we only need to go out the the last option date but
            we can do that refinement at a later date. """

            model.build_tree(t_mat, df_times, df_values)
            v1 = model.callable_puttable_bond_tree(cpn_times, cpn_amounts,
                                                   call_times, call_prices,
                                                   put_times, put_prices,
                                                   face_amount)
            model._num_time_steps += 1
            model.build_tree(t_mat, df_times, df_values)
            v2 = model.callable_puttable_bond_tree(cpn_times, cpn_amounts,
                                                   call_times, call_prices,
                                                   put_times, put_prices,
                                                   face_amount)
            model._num_time_steps -= 1

            v_bondwithoption = (v1['bondwithoption'] + v2['bondwithoption'])/2
            v_bondpure = (v1['bondpure'] + v2['bondpure'])/2

            return {'bondwithoption': v_bondwithoption, 'bondpure': v_bondpure}

        elif isinstance(model, BKTree):

            """ Because we not have a closed form bond price we need to build
            the tree out to the bond maturity which is after option expiry. """

            model.build_tree(t_mat, df_times, df_values)
            v1 = model.callable_puttable_bond_tree(cpn_times, cpn_amounts,
                                                   call_times, call_prices,
                                                   put_times, put_prices,
                                                   face_amount)
            model._num_time_steps += 1
            model.build_tree(t_mat, df_times, df_values)
            v2 = model.callable_puttable_bond_tree(cpn_times, cpn_amounts,
                                                   call_times, call_prices,
                                                   put_times, put_prices,
                                                   face_amount)
            model._num_time_steps -= 1

            v_bondwithoption = (v1['bondwithoption'] + v2['bondwithoption'])/2
            v_bondpure = (v1['bondpure'] + v2['bondpure'])/2

            return {'bondwithoption': v_bondwithoption, 'bondpure': v_bondpure}
        else:
            raise FinError("Unknown model type")

###############################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("ISSUE DATE", self._issue_dt)
        s += label_to_string("MATURITY DATE", self._maturity_dt)
        s += label_to_string("COUPON", self._cpn)
        s += label_to_string("FREQUENCY", self._freq_type)
        s += label_to_string("DAY COUNT TYPE", self._dc_type)
        s += label_to_string("EX-DIV DAYS", self._ex_div_days)

        s += label_to_string("NUM CALL DATES", len(self._call_dts))
        for i in range(0, len(self._call_dts)):
            s += "%12s %12.6f\n" % (self._call_dts[i], self._call_prices[i])

        s += label_to_string("NUM PUT DATES", len(self._put_dts))
        for i in range(0, len(self._put_dts)):
            s += "%12s %12.6f\n" % (self._put_dts[i], self._put_prices[i])

        return s

###############################################################################

    def _print(self):
        print(self)

###############################################################################
