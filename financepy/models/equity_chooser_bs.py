import numpy as np
from scipy import optimize

from ..utils.math import M
from ..utils.global_vars import G_SMALL
from ..utils.global_types import OptionTypes
from .black_scholes_analytic import european_value

DEBUG_MODE = False

#########################################################################################


def _f(ss, *args):
    """Complex chooser option solve for critical stock price that makes the
    forward starting call and put options have the same price on the chooser
    date."""

    t_choose = args[0]
    tc = args[1]
    tp = args[2]
    rtc = args[3]
    rtp = args[4]
    kc = args[5]
    kp = args[6]
    v = args[7]
    qc = args[8]
    qp = args[9]

    call_int = OptionTypes.EUROPEAN_CALL.value
    put_int = OptionTypes.EUROPEAN_PUT.value

    if tc == t_choose:
        v_call = np.maximum(ss - kc, 0.0)
    else:

        v_call = european_value(ss, tc - t_choose, kc, rtc, qc, v, call_int)

    if tp == t_choose:
        v_put = np.maximum(kp - ss, 0.0)
    else:
        v_put = european_value(ss, tp - t_choose, kp, rtp, qp, v, put_int)

    v = v_call - v_put
    return v


#########################################################################################


def equity_chooser_value(
    t_choose, t_call, t_put, k_c, k_p, s0, df_t, df_c, df_p, dq_t, dq_c, dq_p, rfc, qfc, rfp, qfp, vol
):
    """Value the complex chooser option using an approach by Rubinstein
    (1991). See also Haug page 129 for complex chooser options."""

    scalar_input = np.isscalar(s0)

    vol = max(vol, G_SMALL)
    vol2 = vol * vol

    argtuple = (t_choose, t_call, t_put, rfc, rfp, k_c, k_p, vol, qfc, qfp)

    if DEBUG_MODE:
        print("args", argtuple)

    x_init = 0.5 * (k_c + k_p)
    istar = optimize.newton(_f, x0=x_init, args=argtuple, tol=1e-8, maxiter=50)

    if DEBUG_MODE:
        print("istar", istar)

    sqrt_tc = np.sqrt(t_call)
    sqrt_tp = np.sqrt(t_put)
    sqrt_t_choose = np.sqrt(t_choose)

    d1 = (np.log(s0 / istar) + np.log(dq_t / df_t) + 0.5 * vol2 * t_choose) / vol / sqrt_t_choose
    d2 = d1 - vol * sqrt_t_choose

    if DEBUG_MODE:
        print("d1", d1)
        print("d2", d2)

    y1 = (np.log(s0 / k_c) + np.log(dq_c / df_c) + 0.5 * vol2 * t_call) / vol / sqrt_tc
    y2 = (np.log(s0 / k_p) + np.log(dq_p / df_p) + 0.5 * vol2 * t_put) / vol / sqrt_tp

    if DEBUG_MODE:
        print("y1", y1)
        print("y2", y2)

    rho1 = sqrt_t_choose / sqrt_tc
    rho2 = sqrt_t_choose / sqrt_tp

    if DEBUG_MODE:
        print("rho1", rho1)
        print("rho2", rho2)

    if 1 == 1:
        w = s0 * dq_c * M(d1, y1, rho1)
        w = w - k_c * df_c * M(d2, y1 - vol * sqrt_tc, rho1)
        w = w - s0 * dq_p * M(-d1, -y2, rho2)
        w = w + k_p * df_p * M(-d2, -y2 + vol * sqrt_tp, rho2)
    else:
        m1 = np.array([M(a, b, rho1) for a, b in zip(d1, y1)])
        m2 = np.array([M(a, b, rho1) for a, b in zip(d2, y1 - vol * sqrt_tc)])
        m3 = np.array([M(a, b, rho2) for a, b in zip(-d1, -y2)])
        m4 = np.array([M(a, b, rho2) for a, b in zip(-d2, -y2 + vol * sqrt_tp)])

        w = s0 * dq_c * m1
        w -= k_c * df_c * m2
        w -= s0 * dq_p * m3
        w += k_p * df_p * m4

    if scalar_input:
        return w[0]

    return w
