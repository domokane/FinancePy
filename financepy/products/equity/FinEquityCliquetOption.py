##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import numpy as np

from ...utils.frequency import FrequencyTypes
from ...utils.global_variables import gDaysInYear
from ...utils.FinError import FinError
from ...utils.FinGlobalTypes import FinOptionTypes

from ...utils.helper_functions import labelToString, check_argument_types
from ...utils.date import Date
from ...utils.day_count import DayCountTypes
from ...utils.calendar import BusDayAdjustTypes
from ...utils.calendar import CalendarTypes,  DateGenRuleTypes
from ...utils.schedule import Schedule
from ...products.equity.FinEquityOption import FinEquityOption
from ...market.curves.discount_curve_flat import DiscountCurve

from ...models.black_scholes import bsValue, FinModelBlackScholes
from ...models.FinModel import FinModel

###############################################################################
# TODO: Do we need to day count adjust option payoffs ?
# TODO: Monte Carlo pricer
###############################################################################

class FinEquityCliquetOption(FinEquityOption):
    """ A FinEquityCliquetOption is a series of options which start and stop at
    successive times with each subsequent option resetting its strike to be ATM
    at the start of its life. This is also known as a reset option."""

    def __init__(self,
                 start_date: Date,
                 finalExpiryDate: Date,
                 optionType: FinOptionTypes,
                 freq_type: FrequencyTypes,
                 day_count_type: DayCountTypes = DayCountTypes.THIRTY_E_360,
                 calendar_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD):
        """ Create the FinEquityCliquetOption by passing in the start date
        and the end date and whether it is a call or a put. Some additional
        data is needed in order to calculate the individual payments. """

        check_argument_types(self.__init__, locals())

        if optionType != FinOptionTypes.EUROPEAN_CALL and \
           optionType != FinOptionTypes.EUROPEAN_PUT:
            raise FinError("Unknown Option Type" + str(optionType))

        if finalExpiryDate < start_date:
            raise FinError("Expiry date precedes start date")

        self._start_date = start_date
        self._finalExpiryDate = finalExpiryDate
        self._optionType = optionType
        self._freq_type = freq_type
        self._day_count_type = day_count_type
        self._calendar_type = calendar_type
        self._bus_day_adjust_type = bus_day_adjust_type
        self._date_gen_rule_type = date_gen_rule_type

        self._expiry_dates = Schedule(self._start_date,
                                     self._finalExpiryDate,
                                     self._freq_type,
                                     self._calendar_type,
                                     self._bus_day_adjust_type,
                                     self._date_gen_rule_type)._generate()

###############################################################################

    def value(self,
              valuation_date: Date,
              stock_price: float,
              discount_curve: DiscountCurve,
              dividendCurve: DiscountCurve,
              model:FinModel):
        """ Value the cliquet option as a sequence of options using the Black-
        Scholes model. """

        if valuation_date > self._finalExpiryDate:
            raise FinError("Value date after final expiry date.")

        s = stock_price
        v_cliquet = 0.0

        self._v_options = []
        self._dfs = []
        self._actualDates = []

        CALL = FinOptionTypes.EUROPEAN_CALL
        PUT = FinOptionTypes.EUROPEAN_PUT

        if isinstance(model, FinModelBlackScholes):

            v = model._volatility
            v = max(v, 1e-6)
            tprev = 0.0

            for dt in self._expiry_dates:

                if dt > valuation_date:

                    df = discount_curve.df(dt)
                    texp = (dt - valuation_date) / gDaysInYear
                    r = -np.log(df) / texp

                    # option life
                    tau = texp - tprev

                    # The deflator is out to the option reset time
                    dq = dividendCurve._df(tprev)

                    # The option dividend is over the option life
                    dqMat = dividendCurve._df(texp)

                    q = -np.log(dqMat/dq)/tau

                    if self._optionType == CALL:
                        v_fwd_opt = s * dq * bsValue(1.0, tau, 1.0, r, q, v, CALL.value)
                        v_cliquet += v_fwd_opt
                    elif self._optionType == PUT:
                        v_fwd_opt = s * dq * bsValue(1.0, tau, 1.0, r, q, v, PUT.value)
                        v_cliquet += v_fwd_opt
                    else:
                        raise FinError("Unknown option type")

#                    print(dt, r, df, q, v_fwd_opt, v_cliquet)

                    self._dfs.append(df)
                    self._v_options.append(v)
                    self._actualDates.append(dt)
                    tprev = texp
        else:
            raise FinError("Unknown Model Type")

        return v_cliquet

###############################################################################

    def printFlows(self):
        numOptions = len(self._v_options)
        for i in range(0, numOptions):
            print(self._actualDates[i], self._dfs[i], self._v_options[i])

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("START DATE", self._start_date)
        s += labelToString("FINAL EXPIRY DATE", self._finalExpiryDate)
        s += labelToString("OPTION TYPE", self._optionType)
        s += labelToString("FREQUENCY TYPE", self._freq_type)
        s += labelToString("DAY COUNT TYPE", self._day_count_type)
        s += labelToString("CALENDAR TYPE", self._calendar_type)
        s += labelToString("BUS DAY ADJUST TYPE", self._bus_day_adjust_type)
        s += labelToString("DATE GEN RULE TYPE", self._date_gen_rule_type, "")
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
