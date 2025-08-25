# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np
from numba import njit

from ..utils.error import FinError
from ..utils.math import normcdf
from ..market.curves.interpolator import InterpTypes, _uinterpolate
from ..utils.helpers import label_to_string

INTERP = InterpTypes.FLAT_FWD_RATES.value

# dr = theta(t) dt + sigma * dW

########################################################################################


@njit(fastmath=True, cache=True)
def p_fast(t, t_mat, rt, delta, pt, ptd, p_mat, sigma):
    """Forward discount factor as seen at some time t which may be in the
    future for payment at time T where Rt is the delta-period short rate
    seen at time t and pt is the discount factor to time t, ptd is the one
    period discount factor to time t+dt and pT is the discount factor from
    now until the payment of the 1 dollar of the discount factor."""

    bt_t = t_mat - t
    bt_delta = delta
    term1 = np.log(p_mat / pt) - (bt_t / bt_delta) * np.log(ptd / pt)
    term2 = (sigma**2) * t * bt_t * (bt_t - bt_delta) / (2.0)

    log_ahat = term1 - term2
    bhatt_t = (bt_t / bt_delta) * delta
    p = np.exp(log_ahat - bhatt_t * rt)
    return p


########################################################################################


class ModelRatesHoLee:

    ####################################################################################

    def __init__(self, sigma):
        """Construct Ho-Lee model using single parameter of volatility. The
        dynamical equation is dr = theta(t) dt + sigma * dW. Any no-arbitrage
        fitting is done within functions below."""

        if sigma < 0.0:
            raise FinError("Negative volatility not allowed.")

        self._sigma = sigma

    ####################################################################################

    def zcb(self, rt1, t1, t2, discount_curve):

        delta = t2 - t1
        dt = 1e-10
        pt1 = discount_curve.df(t1)
        pt1p = discount_curve.df(t1 + dt)
        pt2 = discount_curve.df(t2)
        z = p_fast(t1, t2, rt1, delta, pt1, pt1p, pt2, self._sigma)
        return z

    ####################################################################################

    def option_on_zcb(
        self, t_exp, t_mat, strike_price, face_amount, df_times, df_values
    ):
        """Price an option on a zero coupon bond using analytical solution of
        Hull-White model. User provides bond face and option strike and expiry
        date and maturity date."""

        if t_exp > t_mat:
            raise FinError("Option expiry after bond matures.")

        if t_exp < 0.0:
            raise FinError("Option expiry time negative.")

        pt_exp = _uinterpolate(t_exp, df_times, df_values, INTERP)
        pt_mat = _uinterpolate(t_mat, df_times, df_values, INTERP)

        sigma = self._sigma

        sigmap = sigma * (t_mat - t_exp) * np.sqrt(t_exp)

        h = (
            np.log((face_amount * pt_mat) / (strike_price * pt_exp)) / sigmap
            + sigmap / 2.0
        )

        call_value = face_amount * pt_mat * normcdf(
            h
        ) - strike_price * pt_exp * normcdf(h - sigmap)
        put_value = strike_price * pt_exp * normcdf(
            -h + sigmap
        ) - face_amount * pt_mat * normcdf(-h)

        return {"call": call_value, "put": put_value}

    ####################################################################################

    def __repr__(self):
        """Return string with class details."""

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("Sigma", self._sigma)
        return s
