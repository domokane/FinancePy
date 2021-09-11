##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

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

from enum import Enum
import numpy as np
from typing import List


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
                 issue_date: Date,
                 maturity_date: Date,  # Date
                 coupon: float,  # Annualised coupon - 0.03 = 3.00%
                 freq_type: FrequencyTypes,
                 accrual_type: DayCountTypes,
                 call_dates: List[Date],
                 call_prices: List[float],
                 put_dates: List[Date],
                 put_prices: List[float],
                 face_amount: float = 100.0):
        """ Create a BondEmbeddedOption object with a maturity date, coupon
        and all of the bond inputs. """

        check_argument_types(self.__init__, locals())

        self._issue_date = issue_date
        self._maturity_date = maturity_date
        self._coupon = coupon
        self._freq_type = freq_type
        self._accrual_type = accrual_type

        self._bond = Bond(issue_date,
                          maturity_date,
                          coupon,
                          freq_type,
                          accrual_type,
                          face_amount)

        # Validate call and put schedules
        for dt in call_dates:
            if dt > self._maturity_date:
                raise FinError("Call date after bond maturity date")

        if len(call_dates) > 0:
            dtprev = call_dates[0]
            for dt in call_dates[1:]:
                if dt <= dtprev:
                    raise FinError("Call dates not increasing")
                else:
                    dtprev = dt

        for dt in put_dates:
            if dt > self._maturity_date:
                raise FinError("Put date after bond maturity date")

        if len(put_dates) > 0:
            dtprev = put_dates[0]
            for dt in put_dates[1:]:
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

        if len(call_dates) != len(call_prices):
            raise FinError("Number of call dates and call prices not the same")

        if len(put_dates) != len(put_prices):
            raise FinError("Number of put dates and put prices not the same")

        self._call_dates = call_dates
        self._call_prices = call_prices
        self._put_dates = put_dates
        self._put_prices = put_prices
        self._face_amount = face_amount
        self._bond._calculate_flow_dates()

###############################################################################

    def value(self,
              settlement_date: Date,
              discount_curve: DiscountCurve,
              model):
        """ Value the bond that settles on the specified date that can have
        both embedded call and put options. This is done using the specified
        model and a discount curve. """

        # Generate bond coupon flow schedule
        cpn = self._bond._coupon/self._bond._frequency
        cpn_times = []
        cpn_amounts = []

        for flow_date in self._bond._flow_dates[1:]:
            if flow_date > settlement_date:
                cpn_time = (flow_date - settlement_date) / gDaysInYear
                cpn_times.append(cpn_time)
                cpn_amounts.append(cpn)

        cpn_times = np.array(cpn_times)
        cpn_amounts = np.array(cpn_amounts)

        # Generate bond call times and prices
        call_times = []
        for dt in self._call_dates:
            if dt > settlement_date:
                call_time = (dt - settlement_date) / gDaysInYear
                call_times.append(call_time)
        call_times = np.array(call_times)
        call_prices = np.array(self._call_prices)

        # Generate bond put times and prices
        put_times = []
        for dt in self._put_dates:
            if dt > settlement_date:
                put_time = (dt - settlement_date) / gDaysInYear
                put_times.append(put_time)
        put_times = np.array(put_times)
        put_prices = np.array(self._put_prices)

        maturity_date = self._bond._maturity_date
        tmat = (maturity_date - settlement_date) / gDaysInYear
        df_times = discount_curve._times
        df_values = discount_curve._dfs

        face_amount = self._bond._face_amount

        if isinstance(model, HWTree):

            """ We need to build the tree out to the bond maturity date. To be
            more precise we only need to go out the the last option date but
            we can do that refinement at a later date. """

            model.build_tree(tmat, df_times, df_values)
            v1 = model.callable_puttable_bond_tree(cpn_times, cpn_amounts,
                                                   call_times, call_prices,
                                                   put_times, put_prices,
                                                   face_amount)
            model._num_time_steps += 1
            model.build_tree(tmat, df_times, df_values)
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

            model.build_tree(tmat, df_times, df_values)
            v1 = model.callable_puttable_bond_tree(cpn_times, cpn_amounts,
                                                   call_times, call_prices,
                                                   put_times, put_prices,
                                                   face_amount)
            model._num_time_steps += 1
            model.build_tree(tmat, df_times, df_values)
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
        s += label_to_string("ISSUE DATE", self._issue_date)
        s += label_to_string("MATURITY DATE", self._maturity_date)
        s += label_to_string("COUPON", self._coupon)
        s += label_to_string("FREQUENCY", self._freq_type)
        s += label_to_string("ACCRUAL TYPE", self._accrual_type)
        s += label_to_string("FACE AMOUNT", self._face_amount)

        s += label_to_string("NUM CALL DATES", len(self._call_dates))
        for i in range(0, len(self._call_dates)):
            s += "%12s %12.6f\n" % (self._call_dates[i], self._call_prices[i])

        s += label_to_string("NUM PUT DATES", len(self._put_dates))
        for i in range(0, len(self._put_dates)):
            s += "%12s %12.6f\n" % (self._put_dates[i], self._put_prices[i])

        return s

###############################################################################

    def _print(self):
        print(self)

###############################################################################
