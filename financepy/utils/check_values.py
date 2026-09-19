import numpy as np


from .date import Date
from .global_vars import G_DAYS_IN_YEAR
from .error import FinError

###########################################################################
# Validation for specific input types that handles arrays and scalars
###########################################################################


def check_curve_dt(anchor_dt, *curves):
    """Check that curves are valid for the requested valuation date."""

    for curve in curves:
        if curve.anchor_dt > anchor_dt:
            raise FinError(
                f"{type(curve).__name__} valuation date {curve.anchor_dt} " f"is after valuation date {anchor_dt}."
            )


###########################################################################


def check_t_exp(value_dt: Date, expiry_dt):
    """Calculate time to expiry in years."""

    if not isinstance(value_dt, Date):
        raise FinError("Valuation date must be a Date.")

    if isinstance(expiry_dt, Date):

        t_exp = (expiry_dt - value_dt) / G_DAYS_IN_YEAR

        if t_exp < 0.0:
            raise FinError("Expiry date falls before valuation date.")

        return t_exp

    elif isinstance(expiry_dt, list):

        t_exps = []

        for exp_dt in expiry_dt:
            t_exp = check_t_exp(value_dt, exp_dt)

            t_exps.append(t_exp)

        return np.array(t_exps)

    else:

        raise FinError("Expiry date must be a Date or list of Dates.")


########################################################################################


def check_stock_price(stock_price):

    s0 = np.asarray(stock_price, dtype=float)

    if s0.ndim > 1:
        raise FinError("Stock price must be a scalar or one-dimensional array.")

    if s0.size == 0:
        raise FinError("Stock price cannot be empty.")

    if not np.all(np.isfinite(s0)):
        raise FinError("Stock price must be finite.")

    if np.any(s0 <= 0.0):
        raise FinError("Stock price must be greater than zero.")


########################################################################################


def check_volatility(volatility):

    v = np.asarray(volatility, dtype=float)

    if v.ndim > 1:
        raise FinError("Volatility must be a scalar or one-dimensional array.")

    if v.size == 0:
        raise FinError("Volatility cannot be empty.")

    if not np.all(np.isfinite(v)):
        raise FinError("Volatility must be finite.")

    if np.any(v <= 0.0):
        raise FinError("Volatility must be greater than zero.")


########################################################################################


def check_strike_price(strike_price):

    k = np.asarray(strike_price, dtype=float)

    if k.ndim > 1:
        raise FinError("Strike price must be a scalar or one-dimensional array.")

    if k.size == 0:
        raise FinError("Strike price cannot be empty.")

    if not np.all(np.isfinite(k)):
        raise FinError("Strike price must be finite.")

    if np.any(k <= 0.0):
        raise FinError("Strike price must be greater than zero.")

    return strike_price


########################################################################################


def check_shapes(*args):
    """Allow scalars and equal-length vectors."""

    shape = None

    for arg in args:

        ndim = np.ndim(arg)

        if ndim > 1:
            raise FinError("Arguments must be scalars or one-dimensional arrays.")

        if ndim == 1:

            arg_shape = np.shape(arg)

            if shape is None:
                shape = arg_shape

            elif arg_shape != shape:
                raise FinError("Vector arguments must have the same length.")


########################################################################################
