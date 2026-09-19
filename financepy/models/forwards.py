import numpy as np
from numba import njit

# -------------------------
# FX forwards
# -------------------------


@njit(cache=True, fastmath=True)
def fx_forward_price(spot, r_dom, r_for, tau):
    return spot * np.exp((r_dom - r_for) * tau)


@njit(cache=True, fastmath=True)
def fx_forward_value(spot, strike, notional, r_dom, r_for, tau):
    fwd = spot * np.exp((r_dom - r_for) * tau)
    return notional * np.exp(-r_dom * tau) * (fwd - strike)


@njit(cache=True, fastmath=True)
def fx_forward_price_df(spot, df_dom, df_for):
    return spot * df_for / df_dom


@njit(cache=True, fastmath=True)
def fx_forward_value_df(spot, strike, notional, df_dom, df_for):
    return notional * (spot * df_for - strike * df_dom)


# -------------------------
# Equity forwards
# -------------------------


@njit(cache=True, fastmath=True)
def equity_forward_price(spot, r, q, tau):
    return spot * np.exp((r - q) * tau)


@njit(cache=True, fastmath=True)
def equity_forward_value(spot, strike, quantity, r, q, tau):
    fwd = spot * np.exp((r - q) * tau)
    return quantity * np.exp(-r * tau) * (fwd - strike)


@njit(cache=True, fastmath=True)
def equity_forward_price_df(spot, discount_df, dividend_df):
    return spot * dividend_df / discount_df


@njit(cache=True, fastmath=True)
def equity_forward_value_df(spot, strike, quantity, df, dividend_df):
    return quantity * (spot * dividend_df - strike * df)


# -------------------------
# Commodity forwards
# -------------------------

# Continuous storage cost u and convenience yield y:
#
# F = S * exp((r + u - y) * T)


@njit(cache=True, fastmath=True)
def commodity_forward_price(
    spot,
    r,
    storage_cost,
    convenience_yield,
    tau,
):
    return spot * np.exp((r + storage_cost - convenience_yield) * tau)


@njit(cache=True, fastmath=True)
def commodity_forward_value(
    spot,
    strike,
    quantity,
    r,
    storage_cost,
    convenience_yield,
    tau,
):
    fwd = spot * np.exp((r + storage_cost - convenience_yield) * tau)
    return quantity * np.exp(-r * tau) * (fwd - strike)


# Discount-factor version:
#
# df             = exp(-r*T)
# storage_df     = exp(-u*T)
# convenience_df = exp(-y*T)
#
# F = S * convenience_df / (df * storage_df)


@njit(cache=True, fastmath=True)
def commodity_forward_price_df(
    spot,
    df,
    storage_df,
    convenience_df,
):
    return spot * convenience_df / (df * storage_df)


@njit(cache=True, fastmath=True)
def commodity_forward_value_df(
    spot,
    strike,
    quantity,
    df,
    storage_df,
    convenience_df,
):
    return quantity * (spot * convenience_df / storage_df - strike * df)
