##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...utils.FinGlobalVariables import gDaysInYear
from ...models.FinModelRatesHW import FinModelRatesHW
from ...models.FinModelRatesBK import FinModelRatesBK
from ...utils.FinError import FinError
from ...utils.Frequency import FinFrequencyTypes
from ...utils.DayCount import FinDayCountTypes
from ...products.bonds.FinBond import FinBond

from ...utils.Date import Date
from ...utils.FinHelperFunctions import labelToString, checkArgumentTypes
from ...market.curves.FinDiscountCurve import FinDiscountCurve

from enum import Enum
import numpy as np
from typing import List


###############################################################################
# TODO: Make it possible to specify start and end of American Callable/Puttable
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


class FinBondEmbeddedOption(object):
    """ Class for fixed coupon bonds with embedded call or put optionality. """

    def __init__(self,
                 issue_date: Date,
                 maturity_date: Date,  # FinDate
                 coupon: float,  # Annualised coupon - 0.03 = 3.00%
                 freq_type: FinFrequencyTypes,
                 accrual_type: FinDayCountTypes,
                 callDates: List[Date],
                 call_prices: List[float],
                 putDates: List[Date],
                 put_prices: List[float],
                 face_amount: float = 100.0):
        """ Create a FinBondEmbeddedOption object with a maturity date, coupon
        and all of the bond inputs. """

        checkArgumentTypes(self.__init__, locals())

        self._issue_date = issue_date
        self._maturity_date = maturity_date
        self._coupon = coupon
        self._freq_type = freq_type
        self._accrual_type = accrual_type

        self._bond = FinBond(issue_date,
                             maturity_date,
                             coupon,
                             freq_type,
                             accrual_type,
                             face_amount)

        # Validate call and put schedules
        for dt in callDates:
            if dt > self._maturity_date:
                raise FinError("Call date after bond maturity date")

        if len(callDates) > 0:
            dtprev = callDates[0]
            for dt in callDates[1:]:
                if dt <= dtprev:
                    raise FinError("Call dates not increasing")
                else:
                    dtprev = dt

        for dt in putDates:
            if dt > self._maturity_date:
                raise FinError("Put date after bond maturity date")

        if len(putDates) > 0:
            dtprev = putDates[0]
            for dt in putDates[1:]:
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

        if len(callDates) != len(call_prices):
            raise FinError("Number of call dates and call prices not the same")

        if len(putDates) != len(put_prices):
            raise FinError("Number of put dates and put prices not the same")

        self._callDates = callDates
        self._call_prices = call_prices
        self._putDates = putDates
        self._put_prices = put_prices
        self._face_amount = face_amount
        self._bond._calculateFlowDates()

###############################################################################

    def value(self,
              settlement_date: Date,
              discount_curve: FinDiscountCurve,
              model):
        """ Value the bond that settles on the specified date that can have
        both embedded call and put options. This is done using the specified
        model and a discount curve. """

        # Generate bond coupon flow schedule
        cpn = self._bond._coupon/self._bond._frequency
        cpnTimes = []
        cpnAmounts = []

        for flowDate in self._bond._flow_dates[1:]:
            if flowDate > settlement_date:
                cpnTime = (flowDate - settlement_date) / gDaysInYear
                cpnTimes.append(cpnTime)
                cpnAmounts.append(cpn)

        cpnTimes = np.array(cpnTimes)
        cpnAmounts = np.array(cpnAmounts)

        # Generate bond call times and prices
        call_times = []
        for dt in self._callDates:
            if dt > settlement_date:
                call_time = (dt - settlement_date) / gDaysInYear
                call_times.append(call_time)
        call_times = np.array(call_times)
        call_prices = np.array(self._call_prices)

        # Generate bond put times and prices
        put_times = []
        for dt in self._putDates:
            if dt > settlement_date:
                put_time = (dt - settlement_date) / gDaysInYear
                put_times.append(put_time)
        put_times = np.array(put_times)
        put_prices = np.array(self._put_prices)

        maturity_date = self._bond._maturity_date
        tmat = (maturity_date - settlement_date) / gDaysInYear
        dfTimes = discount_curve._times
        dfValues = discount_curve._dfs

        face_amount = self._bond._face_amount

        if isinstance(model, FinModelRatesHW):

            """ We need to build the tree out to the bond maturity date. To be
            more precise we only need to go out the the last option date but
            we can do that refinement at a later date. """

            model.buildTree(tmat, dfTimes, dfValues)
            v1 = model.callablePuttableBond_Tree(cpnTimes, cpnAmounts,
                                                 call_times, call_prices,
                                                 put_times, put_prices,
                                                 face_amount)
            model._numTimeSteps += 1
            model.buildTree(tmat, dfTimes, dfValues)
            v2 = model.callablePuttableBond_Tree(cpnTimes, cpnAmounts,
                                                 call_times, call_prices,
                                                 put_times, put_prices,
                                                 face_amount)
            model._numTimeSteps -= 1

            v_bondwithoption = (v1['bondwithoption'] + v2['bondwithoption'])/2
            v_bondpure = (v1['bondpure'] + v2['bondpure'])/2

            return {'bondwithoption': v_bondwithoption, 'bondpure': v_bondpure}

        elif isinstance(model, FinModelRatesBK):

            """ Because we not have a closed form bond price we need to build
            the tree out to the bond maturity which is after option expiry. """

            model.buildTree(tmat, dfTimes, dfValues)
            v1 = model.callablePuttableBond_Tree(cpnTimes, cpnAmounts,
                                                 call_times, call_prices,
                                                 put_times, put_prices,
                                                 face_amount)
            model._numTimeSteps += 1
            model.buildTree(tmat, dfTimes, dfValues)
            v2 = model.callablePuttableBond_Tree(cpnTimes, cpnAmounts,
                                                 call_times, call_prices,
                                                 put_times, put_prices,
                                                 face_amount)
            model._numTimeSteps -= 1

            v_bondwithoption = (v1['bondwithoption'] + v2['bondwithoption'])/2
            v_bondpure = (v1['bondpure'] + v2['bondpure'])/2

            return {'bondwithoption': v_bondwithoption, 'bondpure': v_bondpure}
        else:
            raise FinError("Unknown model type")

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("ISSUE DATE", self._issue_date)
        s += labelToString("MATURITY DATE", self._maturity_date)
        s += labelToString("COUPON", self._coupon)
        s += labelToString("FREQUENCY", self._freq_type)
        s += labelToString("ACCRUAL TYPE", self._accrual_type)
        s += labelToString("FACE AMOUNT", self._face_amount)

        s += labelToString("NUM CALL DATES", len(self._callDates))
        for i in range(0, len(self._callDates)):
            s += "%12s %12.6f\n" % (self._callDates[i], self._call_prices[i])

        s += labelToString("NUM PUT DATES", len(self._putDates))
        for i in range(0, len(self._putDates)):
            s += "%12s %12.6f\n" % (self._putDates[i], self._put_prices[i])

        return s

###############################################################################

    def _print(self):
        print(self)

###############################################################################
