
##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import math

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.global_types import OptionTypes
from ...utils.helpers import label_to_string, option_years
from ...utils.check_values import (
    check_curve_dt,
    check_stock_price,
    check_strike_price,
)

from ...market.curves.discount_curve import DiscountCurve
from ...models.model import Model
from ...models.black_scholes import BlackScholes

from ...models.bs_european_analytic import (
    european_value,
    delta as bs_delta,
    gamma as bs_gamma,
    vega as bs_vega,
    theta as bs_theta,
)


class EquityQuantoVanillaOption:
    """European equity quanto call or put.

    The underlying equity is denominated in foreign currency.
    The payoff is converted into domestic currency at a fixed FX rate.

    FX convention: domestic currency per unit of foreign currency.

    Assumptions:
        - Black-Scholes equity dynamics.
        - Deterministic interest-rate and dividend curves.
        - Constant equity and FX volatilities.
        - Constant equity-FX correlation.
        - European exercise.
    """

    def __init__(
        self,
        expiry_dt: Date,
        strike_price: float,
        opt_type: OptionTypes,
        num_options: float = 1.0,
        quanto_fx_rate: float = 1.0,
    ) -> None:

        if not isinstance(expiry_dt, Date):
            raise FinError("Expiry date must be a Date.")

        if opt_type not in (
            OptionTypes.EUROPEAN_CALL,
            OptionTypes.EUROPEAN_PUT,
        ):
            raise FinError("Quanto option must be European.")

        check_strike_price(strike_price)

        if not math.isfinite(num_options) or num_options <= 0.0:
            raise FinError("Number of options must be positive.")

        if not math.isfinite(quanto_fx_rate) or quanto_fx_rate <= 0.0:
            raise FinError("Quanto FX rate must be positive.")

        self.expiry_dt = expiry_dt
        self.strike_price = float(strike_price)
        self.opt_type = opt_type
        self.opt_type_value = opt_type.value
        self.num_options = num_options
        self.quanto_fx_rate = quanto_fx_rate

    ###########################################################################

    def _quanto_parameters(
        self,
        value_dt: Date,
        stock_price: float,
        domestic_discount_curve: DiscountCurve,
        foreign_discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
        fx_volatility: float,
        correlation: float,
    ):
        """Return the Black-Scholes inputs with quanto-adjusted yield."""

        if not isinstance(value_dt, Date):
            raise FinError("Valuation date must be a Date.")

        if value_dt > self.expiry_dt:
            raise FinError("Valuation date after expiry.")

        check_stock_price(stock_price)

        check_curve_dt(value_dt, domestic_discount_curve)
        check_curve_dt(value_dt, foreign_discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        if not isinstance(model, BlackScholes):
            raise FinError("Quanto option requires BlackScholes.")

        sigma = model.volatility

        if not math.isfinite(sigma) or sigma <= 0.0:
            raise FinError("Equity volatility must be positive.")

        if not math.isfinite(fx_volatility) or fx_volatility < 0.0:
            raise FinError("FX volatility must be non-negative.")

        if not math.isfinite(correlation) or not -1.0 <= correlation <= 1.0:
            raise FinError("Correlation must be between -1 and 1.")

        t = option_years(value_dt, self.expiry_dt)

        if t <= 0.0:
            raise FinError("Option must have positive time to expiry.")

        rd = domestic_discount_curve.zero_rate_cc(self.expiry_dt)
        rf = foreign_discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)

        # Quanto-adjusted dividend yield
        q_quanto = q + rd - rf + correlation * sigma * fx_volatility

        return t, rd, q_quanto, sigma

    ###########################################################################

    def _evaluate(
        self,
        greek,
        value_dt,
        stock_price,
        domestic_discount_curve,
        foreign_discount_curve,
        dividend_curve,
        model,
        fx_volatility,
        correlation,
    ):
        """Common pricing and analytical Greek calculation."""

        t, rd, q_quanto, sigma = self._quanto_parameters(
            value_dt,
            stock_price,
            domestic_discount_curve,
            foreign_discount_curve,
            dividend_curve,
            model,
            fx_volatility,
            correlation,
        )

        s = stock_price
        k = self.strike_price
        opt_type = self.opt_type_value

        scale = self.num_options * self.quanto_fx_rate

        if greek == "value":
            result = european_value(
                s, t, k, rd, q_quanto, sigma, opt_type
            )

        elif greek == "delta":
            result = bs_delta(
                s, t, k, rd, q_quanto, sigma, opt_type
            )

        elif greek == "gamma":
            result = bs_gamma(
                s, t, k, rd, q_quanto, sigma, opt_type
            )

        elif greek == "vega":
            # Total equity vega includes the quanto drift correction.
            vanilla_vega = bs_vega(
                s, t, k, rd, q_quanto, sigma, opt_type
            )

            d = bs_delta(
                s, t, k, rd, q_quanto, sigma, opt_type
            )

            result = (
                vanilla_vega
                - t * s * d * correlation * fx_volatility
            )

        elif greek == "fx_vega":
            d = bs_delta(
                s, t, k, rd, q_quanto, sigma, opt_type
            )

            result = -t * s * d * correlation * sigma

        elif greek == "correlation_delta":
            d = bs_delta(
                s, t, k, rd, q_quanto, sigma, opt_type
            )

            result = -t * s * d * sigma * fx_volatility

        elif greek == "theta":
            # Constant continuously compounded rate convention.
            result = bs_theta(
                s, t, k, rd, q_quanto, sigma, opt_type
            )

        elif greek == "domestic_rho":
            # Domestic-rate changes affect discounting but not
            # the foreign-equity drift under the domestic measure.
            price = european_value(
                s, t, k, rd, q_quanto, sigma, opt_type
            )

            result = -t * price

        elif greek == "foreign_rho":
            d = bs_delta(
                s, t, k, rd, q_quanto, sigma, opt_type
            )

            result = t * s * d

        else:
            raise FinError("Unknown quanto Greek.")

        return result * scale

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        stock_price: float,
        domestic_discount_curve: DiscountCurve,
        foreign_discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
        fx_volatility: float,
        correlation: float,
    ):
        """Value European quanto call or put in domestic currency."""

        return self._evaluate(
            "value",
            value_dt,
            stock_price,
            domestic_discount_curve,
            foreign_discount_curve,
            dividend_curve,
            model,
            fx_volatility,
            correlation,
        )

    ###########################################################################

    # Greek methods accept the same market arguments as value().

    def delta(self, *args, **kwargs):
        """Domestic-currency delta with respect to foreign stock price."""
        return self._evaluate("delta", *args, **kwargs)

    def gamma(self, *args, **kwargs):
        """Second stock-price derivative of domestic option value."""
        return self._evaluate("gamma", *args, **kwargs)

    def vega(self, *args, **kwargs):
        """Total equity volatility sensitivity, including quanto drift."""
        return self._evaluate("vega", *args, **kwargs)

    def fx_vega(self, *args, **kwargs):
        """Sensitivity to FX volatility."""
        return self._evaluate("fx_vega", *args, **kwargs)

    def correlation_delta(self, *args, **kwargs):
        """Sensitivity to equity-FX correlation."""
        return self._evaluate("correlation_delta", *args, **kwargs)

    def theta(self, *args, **kwargs):
        """Calendar-time decay with flat model rates held constant."""
        return self._evaluate("theta", *args, **kwargs)

    def domestic_rho(self, *args, **kwargs):
        """Domestic continuously compounded zero-rate sensitivity."""
        return self._evaluate("domestic_rho", *args, **kwargs)

    def foreign_rho(self, *args, **kwargs):
        """Foreign continuously compounded zero-rate sensitivity."""
        return self._evaluate("foreign_rho", *args, **kwargs)

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("STRIKE PRICE", self.strike_price)
        s += label_to_string("OPTION TYPE", self.opt_type)
        s += label_to_string("NUMBER", self.num_options)
        s += label_to_string("QUANTO FX RATE", self.quanto_fx_rate, "")
        return s

    ###########################################################################

    def _print(self):
        """Print for backward compatibility."""
        print(self)
