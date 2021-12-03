###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np
from enum import Enum

from ...utils.error import FinError
from ...utils.global_vars import gDaysInYear
from ...products.equity.equity_option import EquityOption
from ...models.process_simulator import FinProcessSimulator
from ...market.curves.discount_curve import DiscountCurve
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.date import Date


from ...utils.math import N

# TODO: SOME REDESIGN ON THE MONTE CARLO PROCESS IS PROBABLY NEEDED

###############################################################################


class EquityBarrierTypes(Enum):
    DOWN_AND_OUT_CALL = 1
    DOWN_AND_IN_CALL = 2
    UP_AND_OUT_CALL = 3
    UP_AND_IN_CALL = 4
    UP_AND_OUT_PUT = 5
    UP_AND_IN_PUT = 6
    DOWN_AND_OUT_PUT = 7
    DOWN_AND_IN_PUT = 8

###############################################################################


class EquityBarrierOption(EquityOption):
    """ Class to hold details of an Equity Barrier Option. It also
    calculates the option price using Black Scholes for 8 different
    variants on the Barrier structure in enum EquityBarrierTypes. """

    def __init__(self,
                 expiry_date: Date,
                 strike_price: float,
                 option_type: EquityBarrierTypes,
                 barrier_level: float,
                 num_observations_per_year: (int, float) = 252,
                 notional: float = 1.0):
        """ Create the EquityBarrierOption by specifying the expiry date,
        strike price, option type, barrier level, the number of observations
        per year and the notional. """

        check_argument_types(self.__init__, locals())

        self._expiry_date = expiry_date
        self._strike_price = float(strike_price)
        self._barrier_level = float(barrier_level)
        self._num_observations_per_year = int(num_observations_per_year)

        if option_type not in EquityBarrierTypes:
            raise FinError("Option Type " + str(option_type) + " unknown.")

        self._option_type = option_type
        self._notional = notional

###############################################################################

    def value(self,
              valuation_date: Date,
              stock_price: (float, np.ndarray),
              discount_curve: DiscountCurve,
              dividend_curve: DiscountCurve,
              model):
        """ This prices an Equity Barrier option using the formulae given in
        the paper by Clewlow, Llanos and Strickland December 1994 which can be
        found at

        https://warwick.ac.uk/fac/soc/wbs/subjects/finance/research/wpaperseries/1994/94-54.pdf
        """

        if isinstance(valuation_date, Date) == False:
            raise FinError("Valuation date is not a Date")

        if valuation_date > self._expiry_date:
            raise FinError("Valuation date after expiry date.")

        if discount_curve._valuation_date != valuation_date:
            raise FinError(
                "Discount Curve valuation date not same as option valuation date")

        if dividend_curve._valuation_date != valuation_date:
            raise FinError(
                "Dividend Curve valuation date not same as option valuation date")

        if isinstance(stock_price, int):
            stock_price = float(stock_price)

        if isinstance(stock_price, float):
            stock_prices = [stock_price]
        else:
            stock_prices = stock_price

        values = []
        for s in stock_prices:
            v = self._value_one(valuation_date, s, discount_curve,
                                dividend_curve, model)
            values.append(v)

        if isinstance(stock_price, float):
            return values[0]
        else:
            return np.array(values)

###############################################################################

    def _value_one(self,
                   valuation_date: Date,
                   stock_price: (float, np.ndarray),
                   discount_curve: DiscountCurve,
                   dividend_curve: DiscountCurve,
                   model):
        """ This values a single option. Because of its structure it cannot
        easily be vectorised which is why it has been wrapped. """

        texp = (self._expiry_date - valuation_date) / gDaysInYear

        if texp < 0:
            raise FinError("Option expires before value date.")

        texp = max(texp, 1e-6)

        lnS0k = np.log(stock_price / self._strike_price)
        sqrtT = np.sqrt(texp)

        r = discount_curve.cc_rate(self._expiry_date)
        q = dividend_curve.cc_rate(self._expiry_date)

        k = self._strike_price
        s = stock_price
        h = self._barrier_level

        volatility = model._volatility
        sigmaRootT = volatility * sqrtT
        v2 = volatility * volatility
        mu = r - q
        d1 = (lnS0k + (mu + v2 / 2.0) * texp) / sigmaRootT
        d2 = (lnS0k + (mu - v2 / 2.0) * texp) / sigmaRootT
        df = np.exp(-r * texp)
        dq = np.exp(-q * texp)

        c = s * dq * N(d1) - k * df * N(d2)
        p = k * df * N(-d2) - s * dq * N(-d1)
#        print("CALL:",c,"PUT:",p)

        if self._option_type == EquityBarrierTypes.DOWN_AND_OUT_CALL and s <= h:
            return 0.0
        elif self._option_type == EquityBarrierTypes.UP_AND_OUT_CALL and s >= h:
            return 0.0
        elif self._option_type == EquityBarrierTypes.UP_AND_OUT_PUT and s >= h:
            return 0.0
        elif self._option_type == EquityBarrierTypes.DOWN_AND_OUT_PUT and s <= h:
            return 0.0
        elif self._option_type == EquityBarrierTypes.DOWN_AND_IN_CALL and s <= h:
            return c
        elif self._option_type == EquityBarrierTypes.UP_AND_IN_CALL and s >= h:
            return c
        elif self._option_type == EquityBarrierTypes.UP_AND_IN_PUT and s >= h:
            return p
        elif self._option_type == EquityBarrierTypes.DOWN_AND_IN_PUT and s <= h:
            return p

        num_observations = 1 + texp * self._num_observations_per_year

        # Correction by Broadie, Glasserman and Kou, Mathematical Finance, 1997
        # Adjusts the barrier for discrete and not continuous observations
        h_adj = h
        t = texp / num_observations

        if self._option_type == EquityBarrierTypes.DOWN_AND_OUT_CALL:
            h_adj = h * np.exp(-0.5826 * volatility * np.sqrt(t))
        elif self._option_type == EquityBarrierTypes.DOWN_AND_IN_CALL:
            h_adj = h * np.exp(-0.5826 * volatility * np.sqrt(t))
        elif self._option_type == EquityBarrierTypes.UP_AND_IN_CALL:
            h_adj = h * np.exp(0.5826 * volatility * np.sqrt(t))
        elif self._option_type == EquityBarrierTypes.UP_AND_OUT_CALL:
            h_adj = h * np.exp(0.5826 * volatility * np.sqrt(t))
        elif self._option_type == EquityBarrierTypes.UP_AND_IN_PUT:
            h_adj = h * np.exp(0.5826 * volatility * np.sqrt(t))
        elif self._option_type == EquityBarrierTypes.UP_AND_OUT_PUT:
            h_adj = h * np.exp(0.5826 * volatility * np.sqrt(t))
        elif self._option_type == EquityBarrierTypes.DOWN_AND_OUT_PUT:
            h_adj = h * np.exp(-0.5826 * volatility * np.sqrt(t))
        elif self._option_type == EquityBarrierTypes.DOWN_AND_IN_PUT:
            h_adj = h * np.exp(-0.5826 * volatility * np.sqrt(t))
        else:
            raise FinError("Unknown barrier option type." +
                           str(self._option_type))

        h = h_adj

        if abs(volatility) < 1e-5:
            volatility = 1e-5

        l = (mu + v2 / 2.0) / v2
        y = np.log(h * h / (s * k)) / sigmaRootT + l * sigmaRootT
        x1 = np.log(s / h) / sigmaRootT + l * sigmaRootT
        y1 = np.log(h / s) / sigmaRootT + l * sigmaRootT
        hOverS = h / s

        if self._option_type == EquityBarrierTypes.DOWN_AND_OUT_CALL:
            if h >= k:
                c_do = s * dq * N(x1) - k * df * N(x1 - sigmaRootT) \
                    - s * dq * pow(hOverS, 2.0 * l) * N(y1) \
                    + k * df * pow(hOverS, 2.0 * l - 2.0) * N(y1 - sigmaRootT)
                price = c_do
            else:
                c_di = s * dq * pow(hOverS, 2.0 * l) * N(y) \
                    - k * df * pow(hOverS, 2.0 * l - 2.0) * N(y - sigmaRootT)
                price = c - c_di
        elif self._option_type == EquityBarrierTypes.DOWN_AND_IN_CALL:
            if h <= k:
                c_di = s * dq * pow(hOverS, 2.0 * l) * N(y) \
                    - k * df * pow(hOverS, 2.0 * l - 2.0) * N(y - sigmaRootT)
                price = c_di
            else:
                c_do = s * dq * N(x1) \
                    - k * df * N(x1 - sigmaRootT) \
                    - s * dq * pow(hOverS, 2.0 * l) * N(y1) \
                    + k * df * pow(hOverS, 2.0 * l - 2.0) * N(y1 - sigmaRootT)
                price = c - c_do
        elif self._option_type == EquityBarrierTypes.UP_AND_IN_CALL:
            if h >= k:
                c_ui = s * dq * N(x1) - k * df * N(x1 - sigmaRootT) \
                    - s * dq * pow(hOverS, 2.0 * l) * (N(-y) - N(-y1)) \
                    + k * df * pow(hOverS, 2.0 * l - 2.0) * \
                    (N(-y + sigmaRootT) - N(-y1 + sigmaRootT))
                price = c_ui
            else:
                price = c
        elif self._option_type == EquityBarrierTypes.UP_AND_OUT_CALL:
            if h > k:
                c_ui = s * dq * N(x1) - k * df * N(x1 - sigmaRootT) \
                    - s * dq * pow(hOverS, 2.0 * l) * (N(-y) - N(-y1)) \
                    + k * df * pow(hOverS, 2.0 * l - 2.0) * \
                    (N(-y + sigmaRootT) - N(-y1 + sigmaRootT))
                price = c - c_ui
            else:
                price = 0.0
        elif self._option_type == EquityBarrierTypes.UP_AND_IN_PUT:
            if h > k:
                p_ui = -s * dq * pow(hOverS, 2.0 * l) * N(-y) \
                    + k * df * pow(hOverS, 2.0 * l - 2.0) * N(-y + sigmaRootT)
                price = p_ui
            else:
                p_uo = -s * dq * N(-x1) \
                    + k * df * N(-x1 + sigmaRootT) \
                    + s * dq * pow(hOverS, 2.0 * l) * N(-y1) \
                       - k * df * pow(hOverS, 2.0 * l - 2.0) * \
                    N(-y1 + sigmaRootT)
                price = p - p_uo
        elif self._option_type == EquityBarrierTypes.UP_AND_OUT_PUT:
            if h >= k:
                p_ui = -s * dq * pow(hOverS, 2.0 * l) * N(-y) \
                    + k * df * pow(hOverS, 2.0 * l - 2.0) * N(-y + sigmaRootT)
                price = p - p_ui
            else:
                p_uo = -s * dq * N(-x1) \
                    + k * df * N(-x1 + sigmaRootT) \
                    + s * dq * pow(hOverS, 2.0 * l) * N(-y1) \
                       - k * df * pow(hOverS, 2.0 * l - 2.0) * \
                    N(-y1 + sigmaRootT)
                price = p_uo
        elif self._option_type == EquityBarrierTypes.DOWN_AND_OUT_PUT:
            if h >= k:
                price = 0.0
            else:
                p_di = -s * dq * N(-x1) \
                    + k * df * N(-x1 + sigmaRootT) \
                    + s * dq * pow(hOverS, 2.0 * l) * (N(y) - N(y1)) \
                       - k * df * pow(hOverS, 2.0 * l - 2.0) * \
                    (N(y - sigmaRootT) - N(y1 - sigmaRootT))
                price = p - p_di
        elif self._option_type == EquityBarrierTypes.DOWN_AND_IN_PUT:
            if h >= k:
                price = p
            else:
                p_di = -s * dq * N(-x1) \
                    + k * df * N(-x1 + sigmaRootT) \
                    + s * dq * pow(hOverS, 2.0 * l) * (N(y) - N(y1)) \
                       - k * df * pow(hOverS, 2.0 * l - 2.0) * \
                    (N(y - sigmaRootT) - N(y1 - sigmaRootT))
                price = p_di
        else:
            raise FinError("Unknown barrier option type." +
                           str(self._option_type))

        v = price * self._notional
        return v

###############################################################################

    def value_mc(self,
                 valuation_date: Date,
                 stock_price: float,
                 discount_curve: DiscountCurve,
                 dividend_curve: DiscountCurve,
                 process_type,
                 model_params,
                 numAnnObs: int = 252,
                 num_paths: int = 10000,
                 seed: int = 4242):
        """ A Monte-Carlo based valuation of the barrier option which simulates
        the evolution of the stock price of at a specified number of annual
        observation times until expiry to examine if the barrier has been
        crossed and the corresponding value of the final payoff, if any. It
        assumes a GBM model for the stock price. """

        texp = (self._expiry_date - valuation_date) / gDaysInYear
        num_time_steps = int(texp * numAnnObs)
        K = self._strike_price
        B = self._barrier_level
        option_type = self._option_type

        process = FinProcessSimulator()

        r = discount_curve.zero_rate(self._expiry_date)

        # TODO - NEED TO DECIDE IF THIS IS PART OF MODEL PARAMS OR NOT ??????????????

        r = discount_curve.cc_rate(self._expiry_date)
        q = dividend_curve.cc_rate(self._expiry_date)

        #######################################################################

        if option_type == EquityBarrierTypes.DOWN_AND_OUT_CALL and stock_price <= B:
            return 0.0
        elif option_type == EquityBarrierTypes.UP_AND_OUT_CALL and stock_price >= B:
            return 0.0
        elif option_type == EquityBarrierTypes.DOWN_AND_OUT_PUT and stock_price <= B:
            return 0.0
        elif option_type == EquityBarrierTypes.UP_AND_OUT_PUT and stock_price >= B:
            return 0.0

        #######################################################################

        simple_call = False
        simple_put = False

        if option_type == EquityBarrierTypes.DOWN_AND_IN_CALL and stock_price <= B:
            simple_call = True
        elif option_type == EquityBarrierTypes.UP_AND_IN_CALL and stock_price >= B:
            simple_call = True
        elif option_type == EquityBarrierTypes.UP_AND_IN_PUT and stock_price >= B:
            simple_put = True
        elif option_type == EquityBarrierTypes.DOWN_AND_IN_PUT and stock_price <= B:
            simple_put = True

        if simple_put or simple_call:
            Sall = process.get_process(
                process_type, texp, model_params, 1, num_paths, seed)

        if simple_call:
            c = (np.maximum(Sall[:, -1] - K, 0.0)).mean()
            c = c * np.exp(-r * texp)
            return c

        if simple_put:
            p = (np.maximum(K - Sall[:, -1], 0.0)).mean()
            p = p * np.exp(-r * texp)
            return p

        # Get full set of paths
        Sall = process.get_process(process_type, texp, model_params, num_time_steps,
                                   num_paths, seed)

        (num_paths, num_time_steps) = Sall.shape

        if option_type == EquityBarrierTypes.DOWN_AND_IN_CALL or \
           option_type == EquityBarrierTypes.DOWN_AND_OUT_CALL or \
           option_type == EquityBarrierTypes.DOWN_AND_IN_PUT or \
           option_type == EquityBarrierTypes.DOWN_AND_OUT_PUT:

            barrierCrossedFromAbove = [False] * num_paths

            for p in range(0, num_paths):
                barrierCrossedFromAbove[p] = np.any(Sall[p] <= B)

        if option_type == EquityBarrierTypes.UP_AND_IN_CALL or \
           option_type == EquityBarrierTypes.UP_AND_OUT_CALL or \
           option_type == EquityBarrierTypes.UP_AND_IN_PUT or \
           option_type == EquityBarrierTypes.UP_AND_OUT_PUT:

            barrierCrossedFromBelow = [False] * num_paths
            for p in range(0, num_paths):
                barrierCrossedFromBelow[p] = np.any(Sall[p] >= B)

        payoff = np.zeros(num_paths)
        ones = np.ones(num_paths)

        if option_type == EquityBarrierTypes.DOWN_AND_OUT_CALL:
            payoff = np.maximum(Sall[:, -1] - K, 0.0) * \
                (ones - barrierCrossedFromAbove)
        elif option_type == EquityBarrierTypes.DOWN_AND_IN_CALL:
            payoff = np.maximum(Sall[:, -1] - K, 0.0) * barrierCrossedFromAbove
        elif option_type == EquityBarrierTypes.UP_AND_IN_CALL:
            payoff = np.maximum(Sall[:, -1] - K, 0.0) * barrierCrossedFromBelow
        elif option_type == EquityBarrierTypes.UP_AND_OUT_CALL:
            payoff = np.maximum(Sall[:, -1] - K, 0.0) * \
                (ones - barrierCrossedFromBelow)
        elif option_type == EquityBarrierTypes.UP_AND_IN_PUT:
            payoff = np.maximum(K - Sall[:, -1], 0.0) * barrierCrossedFromBelow
        elif option_type == EquityBarrierTypes.UP_AND_OUT_PUT:
            payoff = np.maximum(K - Sall[:, -1], 0.0) * \
                (ones - barrierCrossedFromBelow)
        elif option_type == EquityBarrierTypes.DOWN_AND_OUT_PUT:
            payoff = np.maximum(K - Sall[:, -1], 0.0) * \
                (ones - barrierCrossedFromAbove)
        elif option_type == EquityBarrierTypes.DOWN_AND_IN_PUT:
            payoff = np.maximum(K - Sall[:, -1], 0.0) * barrierCrossedFromAbove
        else:
            raise FinError("Unknown barrier option type." +
                           str(self._option_type))

        v = payoff.mean() * np.exp(- r * texp)

        return v * self._notional

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self._expiry_date)
        s += label_to_string("STRIKE PRICE", self._strike_price)
        s += label_to_string("OPTION TYPE", self._option_type)
        s += label_to_string("BARRIER LEVEL", self._barrier_level)
        s += label_to_string("NUM OBSERVATIONS",
                             self._num_observations_per_year)
        s += label_to_string("NOTIONAL", self._notional, "")
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
