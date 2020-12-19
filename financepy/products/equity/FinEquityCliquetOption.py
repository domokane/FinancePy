##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import numpy as np

from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinError import FinError
from ...finutils.FinGlobalTypes import FinOptionTypes

from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes
from ...finutils.FinDate import FinDate
from ...finutils.FinDayCount import FinDayCountTypes
from ...finutils.FinCalendar import FinBusDayAdjustTypes
from ...finutils.FinCalendar import FinCalendarTypes,  FinDateGenRuleTypes
from ...finutils.FinSchedule import FinSchedule
from ...products.equity.FinEquityOption import FinEquityOption
from ...market.curves.FinDiscountCurveFlat import FinDiscountCurve

from ...models.FinModelBlackScholes import bsValue, FinModelBlackScholes
from ...models.FinModel import FinModel


###############################################################################
# TODO: Do we need to day count adjust option payoffs ?
# TODO: Monte Carlo pricer
###############################################################################


###############################################################################


class FinEquityCliquetOption(FinEquityOption):
    ''' A FinEquityCliquetOption is a series of options which start and stop at
    successive times with each subsequent option resetting its strike to be ATM
    at the start of its life. This is also known as a reset option.'''

    def __init__(self,
                 startDate: FinDate,
                 finalExpiryDate: FinDate,
                 optionType: FinOptionTypes,
                 freqType: FinFrequencyTypes,
                 dayCountType: FinDayCountTypes = FinDayCountTypes.THIRTY_E_360,
                 calendarType: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                 busDayAdjustType: FinBusDayAdjustTypes = FinBusDayAdjustTypes.FOLLOWING,
                 dateGenRuleType: FinDateGenRuleTypes = FinDateGenRuleTypes.BACKWARD):
        ''' Create the FinEquityCliquetOption by passing in the start date
        and the end date and whether it is a call or a put. Some additional
        data is needed in order to calculate the individual payments. '''

        checkArgumentTypes(self.__init__, locals())

        if optionType != FinOptionTypes.EUROPEAN_CALL and \
           optionType != FinOptionTypes.EUROPEAN_PUT:
            raise FinError("Unknown Option Type" + str(optionType))

        if finalExpiryDate < startDate:
            raise FinError("Expiry date precedes start date")

        self._startDate = startDate
        self._finalExpiryDate = finalExpiryDate
        self._optionType = optionType
        self._freqType = freqType
        self._dayCountType = dayCountType
        self._calendarType = calendarType
        self._busDayAdjustType = busDayAdjustType
        self._dateGenRuleType = dateGenRuleType

        self._expiryDates = FinSchedule(self._startDate,
                                        self._finalExpiryDate,
                                        self._freqType,
                                        self._calendarType,
                                        self._busDayAdjustType,
                                        self._dateGenRuleType)._generate()

###############################################################################

    def value(self,
              valueDate: FinDate,
              stockPrice: float,
              discountCurve: FinDiscountCurve,
              dividendYield: float,
              model:FinModel):
        ''' Value the cliquet option as a sequence of options using the Black-
        Scholes model. '''

        if valueDate > self._finalExpiryDate:
            raise FinError("Value date after final expiry date.")

        s0 = stockPrice
        q = dividendYield
        v_cliquet = 0.0

        self._v_options = []
        self._dfs = []
        self._actualDates = []

        if isinstance(model, FinModelBlackScholes):

            vol = model._volatility
            vol = max(vol, 1e-6)
            tprev = 0.0

            for dt in self._expiryDates:

                if dt > valueDate:
                    df = discountCurve.df(dt)
                    t = (dt - valueDate) / gDaysInYear
                    r = -np.log(df) / t
                    texp = t - tprev
                    dq = np.exp(-q * tprev)

                    if self._optionType == FinOptionTypes.EUROPEAN_CALL:
                        v = s0 * dq * bsValue(1.0, texp, 1.0, r, q, vol, FinOptionTypes.EUROPEAN_CALL.value)
                        v_cliquet += v
                    elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
                        v = s0 * dq * bsValue(1.0, texp, 1.0, r, q, vol, FinOptionTypes.EUROPEAN_PUT.value)
                        v_cliquet += v
                    else:
                        raise FinError("Unknown option type")

                    self._dfs.append(df)
                    self._v_options.append(v)
                    self._actualDates.append(dt)
                    tprev = t

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
        s += labelToString("START DATE", self._startDate)
        s += labelToString("FINAL EXPIRY DATE", self._finalExpiryDate)
        s += labelToString("OPTION TYPE", self._optionType)
        s += labelToString("FREQUENCY TYPE", self._freqType)
        s += labelToString("DAY COUNT TYPE", self._dayCountType)
        s += labelToString("CALENDAR TYPE", self._calendarType)
        s += labelToString("BUS DAY ADJUST TYPE", self._busDayAdjustType)
        s += labelToString("DATE GEN RULE TYPE", self._dateGenRuleType, "")
        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
