##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np


from ...utils.math import normcdf
from ...utils.global_vars import G_SMALL
from ...utils.error import FinError
from ...utils.date import Date

from ...models.gbm_process_simulator import get_paths_times
from ...products.equity.equity_option import EquityOption
from ...utils.helpers import label_to_string, check_argument_types
from ...market.curves.discount_curve import DiscountCurve
from ...utils.global_types import OptionTypes
from ...utils.helpers import option_years
from ...utils.check_values import check_curve_dt
from ...models.model import Model

##########################################################################
# TODO: Attempt control variate adjustment to monte carlo
# TODO: Sobol for Monte Carlo
# TODO: TIGHTEN UP LIMIT FOR W FROM 100
# TODO: Vectorise the analytical pricing formula
##########################################################################


##########################################################################
# FLOAT STRIKE LOOKBACK CALL PAYS MAX(S(T)-SMIN,0)
# FLOAT STRIKE LOOKBACK PUT PAYS MAX(SMAX-S(T),0)
##########################################################################


class EquityFloatLookbackOption(EquityOption):
    """This is an equity option in which the strike of the option is not fixed
    but is set at expiry to equal the minimum stock price in the case of a call
    or the maximum stock price in the case of a put. In other words the buyer
    of the call gets to buy the asset at the lowest price over the period
    before expiry while the buyer of the put gets to sell the asset at the
    highest price before expiry."""

    def __init__(self, expiry_dt: Date, opt_type: OptionTypes) -> None:
        """Create the FloatLookbackOption by specifying the expiry date and
        the OPTION_TYPE. The strike is determined internally as the maximum or
        minimum of the stock price depending on whether it is a put or a call
        option."""

        check_argument_types(self.__init__, locals())

        if opt_type not in [
            OptionTypes.EUROPEAN_CALL,
            OptionTypes.EUROPEAN_PUT,
        ]:
            raise FinError("OPTION_TYPE must be EUROPEAN_CALL or EUROPEAN_PUT")

        self.expiry_dt = expiry_dt
        self.opt_type = opt_type

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
        stock_min_max: float,
    ):
        """Valuation of the Floating Lookback option using Black-Scholes using
        the formulae derived by Goldman, Sosin and Gatto (1979)."""

        t_exp = option_years(value_dt, self.expiry_dt)

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        volatility = model.volatility

        if volatility < 0.0:
            raise FinError("Volatility must be non-negative.")

        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)

        v = volatility
        s0 = stock_price
        smin = 0.0
        smax = 0.0

        if self.opt_type == OptionTypes.EUROPEAN_CALL:
            smin = stock_min_max
            if smin > s0:
                raise FinError("Smin must be less than or equal to the stock price.")
        elif self.opt_type == OptionTypes.EUROPEAN_PUT:
            smax = stock_min_max
            if smax < s0:
                raise FinError("Smax must be greater than or equal to the stock price.")

        if abs(r - q) < G_SMALL:
            q = r + G_SMALL

        dq = np.exp(-q * t_exp)
        df = np.exp(-r * t_exp)
        b = r - q
        u = v * v / 2.0 / b
        w = 2.0 * b / v / v
        expbt = np.exp(b * t_exp)

        # Taken from Haug Page 142
        if self.opt_type == OptionTypes.EUROPEAN_CALL:

            a1 = (np.log(s0 / smin) + (b + (v**2) / 2.0) * t_exp) / v / np.sqrt(t_exp)
            a2 = a1 - v * np.sqrt(t_exp)

            if smin == s0:
                term = normcdf(-a1 + 2.0 * b * np.sqrt(t_exp) / v) - expbt * normcdf(-a1)
            elif s0 < smin and w < -100:
                term = -expbt * normcdf(-a1)
            else:
                term = ((s0 / smin) ** (-w)) * normcdf(-a1 + 2.0 * b * np.sqrt(t_exp) / v) - expbt * normcdf(-a1)

            v = s0 * dq * normcdf(a1) - smin * df * normcdf(a2) + s0 * df * u * term

        elif self.opt_type == OptionTypes.EUROPEAN_PUT:

            b1 = (np.log(s0 / smax) + (b + (v**2) / 2.0) * t_exp) / v / np.sqrt(t_exp)
            b2 = b1 - v * np.sqrt(t_exp)

            if smax == s0:
                term = -normcdf(b1 - 2.0 * b * np.sqrt(t_exp) / v) + expbt * normcdf(b1)
            elif s0 < smax and w > 100:
                term = expbt * normcdf(b1)
            else:
                term = (-((s0 / smax) ** (-w))) * normcdf(b1 - 2.0 * b * np.sqrt(t_exp) / v) + expbt * normcdf(b1)

            v = smax * df * normcdf(-b2) - s0 * dq * normcdf(-b1) + s0 * df * u * term

        else:
            raise FinError("Unknown lookback OPTION_TYPE:" + str(self.opt_type))

        if 1 == 0:
            print("\nLOOKBACK DEBUG")
            print("t_exp :", t_exp, type(t_exp), np.shape(t_exp))
            print("r     :", r, type(r), np.shape(r))
            print("q     :", q, type(q), np.shape(q))
            print("dq    :", dq, type(dq), np.shape(dq))
            print("df    :", df, type(df), np.shape(df))
            print("b     :", b, type(b), np.shape(b))
            print("u     :", u, type(u), np.shape(u))
            print("w     :", w, type(w), np.shape(w))
            print("expbt :", expbt, type(expbt), np.shape(expbt))
            print("a1    :", a1, type(a1), np.shape(a1))
            print("a2    :", a2, type(a2), np.shape(a2))
            print("term  :", term, type(term), np.shape(term))
            print("v     :", v, type(v), np.shape(v))

        return v

    ###########################################################################

    def value_mc(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
        stock_min_max: float,
        num_paths: int = 10000,
        num_steps_per_year: int = 252,
        seed: int = 4242,
    ):
        """Monte Carlo valuation of a floating strike lookback option using a
        Black-Scholes model that assumes the stock follows a GBM process."""

        t_exp = option_years(value_dt, self.expiry_dt)

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        if model.volatility < 0.0:
            raise FinError("Volatility must be non-negative.")

        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)
        df = discount_curve.df(self.expiry_dt)

        num_time_steps = max(
            1,
            int(np.ceil(t_exp * num_steps_per_year)),
        )

        mu = r - q

        _, s_all = get_paths_times(
            num_paths,
            num_time_steps,
            t_exp,
            mu,
            stock_price,
            model.volatility,
            seed,
        )

        opt_type = self.opt_type

        if opt_type == OptionTypes.EUROPEAN_CALL:

            s_min_vector = np.minimum(
                np.min(s_all, axis=1),
                stock_min_max,
            )

            payoff = np.maximum(
                s_all[:, -1] - s_min_vector,
                0.0,
            )

        elif opt_type == OptionTypes.EUROPEAN_PUT:

            s_max_vector = np.maximum(
                np.max(s_all, axis=1),
                stock_min_max,
            )

            payoff = np.maximum(
                s_max_vector - s_all[:, -1],
                0.0,
            )
        else:
            raise FinError("Unknown lookback OPTION_TYPE:" + str(opt_type))

        v = payoff.mean() * df
        return v

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("OPTION_TYPE", self.opt_type, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
