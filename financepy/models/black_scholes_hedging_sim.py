##############################################################################
# Black-Scholes delta-hedging simulation
##############################################################################

from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd

from numba import njit

from financepy.models.black_scholes_analytic import european_value
from financepy.models.black_scholes_analytic import delta

from financepy.utils.global_types import OptionTypes
from financepy.utils.error import FinError

##############################################################################
# Constants
##############################################################################

EUROPEAN_CALL = OptionTypes.EUROPEAN_CALL.value
EUROPEAN_PUT = OptionTypes.EUROPEAN_PUT.value


##############################################################################
# Results
##############################################################################


@dataclass
class HedgePathResult:
    """Results from one delta-hedged European option path."""

    history: Optional[pd.DataFrame]
    hedging_error: float
    terminal_stock_price: float
    option_payoff: float
    initial_option_value: float
    total_transaction_costs: float


@dataclass
class HedgeSimulationResult:
    """Results from a multi-path delta-hedging simulation."""

    hedging_errors: np.ndarray
    terminal_stock_prices: np.ndarray
    option_payoffs: np.ndarray
    transaction_costs: np.ndarray
    initial_option_value: float

    @property
    def num_paths(self) -> int:
        return len(self.hedging_errors)

    @property
    def mean_hedging_error(self) -> float:
        return float(np.mean(self.hedging_errors))

    @property
    def std_hedging_error(self) -> float:
        return float(
            np.std(
                self.hedging_errors,
                ddof=1,
            )
        )

    @property
    def rmse(self) -> float:
        return float(np.sqrt(np.mean(self.hedging_errors**2)))

    @property
    def mean_transaction_cost(self) -> float:
        return float(np.mean(self.transaction_costs))

    def summary(self) -> pd.Series:
        errors = self.hedging_errors

        return pd.Series(
            {
                "num_paths": len(errors),
                "initial_option_value": self.initial_option_value,
                "mean_hedging_error": np.mean(errors),
                "std_hedging_error": np.std(errors, ddof=1),
                "rmse": np.sqrt(np.mean(errors**2)),
                "min_hedging_error": np.min(errors),
                "p01_hedging_error": np.percentile(errors, 1.0),
                "p05_hedging_error": np.percentile(errors, 5.0),
                "median_hedging_error": np.median(errors),
                "p95_hedging_error": np.percentile(errors, 95.0),
                "p99_hedging_error": np.percentile(errors, 99.0),
                "max_hedging_error": np.max(errors),
                "mean_transaction_cost": np.mean(self.transaction_costs),
                "mean_terminal_stock_price": np.mean(self.terminal_stock_prices),
            }
        )

    ###########################################################################

    def plot_hedging_error_vs_terminal_stock_price(self):
        """Plot terminal stock price against terminal hedging error."""

        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 6))

        ax.scatter(
            self.terminal_stock_prices,
            self.hedging_errors,
            s=5,
            alpha=0.25,
        )

        ax.axhline(
            0.0,
            linewidth=1.0,
            linestyle="--",
        )

        ax.set_xlabel("Terminal Stock Price $S_T$")
        ax.set_ylabel("Hedging Error")
        ax.set_title("Black-Scholes Hedging Error vs Terminal Stock Price")
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        return fig, ax


##############################################################################
# Validation
##############################################################################


def _validate_inputs(
    num_options: int,
    option_type_int: int,
    stock_price: float,
    strike_price: float,
    implied_volatility: float,
    realized_volatility: float,
    time_to_expiry: float,
    num_steps: int,
    transaction_cost_rate: float,
):
    if num_options <= 0:
        raise FinError("Number of options must be positive.")

    if option_type_int not in (
        EUROPEAN_CALL,
        EUROPEAN_PUT,
    ):
        raise FinError("Only European calls and puts are supported.")

    if stock_price <= 0.0:
        raise FinError("Stock price must be positive.")

    if strike_price <= 0.0:
        raise FinError("Strike price must be positive.")

    if implied_volatility <= 0.0:
        raise FinError("Implied volatility must be positive.")

    if realized_volatility < 0.0:
        raise FinError("Realized volatility cannot be negative.")

    if time_to_expiry <= 0.0:
        raise FinError("Time to expiry must be positive.")

    if num_steps <= 0:
        raise FinError("Number of steps must be positive.")

    if transaction_cost_rate < 0.0:
        raise FinError("Transaction cost rate cannot be negative.")


##############################################################################
# Numba helpers
##############################################################################


@njit(cache=True)
def _option_payoff(
    stock_price,
    strike_price,
    option_type_int,
):
    if option_type_int == EUROPEAN_CALL:
        return max(
            stock_price - strike_price,
            0.0,
        )

    return max(
        strike_price - stock_price,
        0.0,
    )


##############################################################################
# Core hedge engine
##############################################################################


@njit(cache=True)
def _simulate_single_path_kernel(
    normal_draws,
    num_options,
    option_type_int,
    stock_price_0,
    strike_price,
    risk_free_rate,
    dividend_yield,
    implied_volatility,
    realized_volatility,
    time_to_expiry,
    stock_drift,
    transaction_cost_rate,
    store_history,
    history,
):
    """
    Canonical single-path hedge engine.

    ALL hedge accounting lives here.

    If store_history is True, history is filled with one row per hedge date.

    History columns
    ---------------
     0 step
     1 time
     2 remaining_time
     3 normal_draw
     4 stock_price
     5 option_value
     6 option_position_value
     7 target_delta
     8 shares_held
     9 shares_traded
    10 stock_trade_notional
    11 stock_trade_cash_flow
    12 transaction_cost
    13 cumulative_transaction_cost
    14 interest_accrued
    15 bank_after_trade
    16 hedge_portfolio_value
    17 option_payoff
    18 net_hedged_value
    """

    num_steps = len(normal_draws)

    dt = time_to_expiry / num_steps

    sqrt_dt = np.sqrt(dt)

    growth_factor = np.exp(risk_free_rate * dt)

    # ------------------------------------------------------------------
    # Initial option price and delta
    # ------------------------------------------------------------------
    option_value = european_value(
        stock_price_0,
        time_to_expiry,
        strike_price,
        risk_free_rate,
        dividend_yield,
        implied_volatility,
        option_type_int,
    )

    initial_option_value = option_value

    target_delta = delta(
        stock_price_0,
        time_to_expiry,
        strike_price,
        risk_free_rate,
        dividend_yield,
        implied_volatility,
        option_type_int,
    )

    stock_price = stock_price_0

    option_position_value = option_value * num_options

    # ------------------------------------------------------------------
    # Initial hedge
    # ------------------------------------------------------------------
    shares_before = 0.0

    shares = target_delta * num_options

    shares_traded = shares - shares_before

    stock_trade_notional = abs(shares_traded) * stock_price

    transaction_cost = transaction_cost_rate * stock_trade_notional

    total_transaction_costs = transaction_cost

    stock_trade_cash_flow = -shares_traded * stock_price

    # Writer receives option premium.
    bank_balance = option_position_value + stock_trade_cash_flow - transaction_cost

    hedge_portfolio_value = bank_balance + shares * stock_price

    net_hedged_value = hedge_portfolio_value - option_position_value

    # ------------------------------------------------------------------
    # Store initial state
    # ------------------------------------------------------------------
    if store_history:

        history[0, 0] = 0.0
        history[0, 1] = 0.0
        history[0, 2] = time_to_expiry
        history[0, 3] = np.nan
        history[0, 4] = stock_price
        history[0, 5] = option_value
        history[0, 6] = option_position_value
        history[0, 7] = target_delta
        history[0, 8] = shares
        history[0, 9] = shares_traded
        history[0, 10] = stock_trade_notional
        history[0, 11] = stock_trade_cash_flow
        history[0, 12] = 0.0
        history[0, 13] = transaction_cost
        history[0, 14] = total_transaction_costs
        history[0, 15] = 0.0
        history[0, 16] = bank_balance
        history[0, 17] = hedge_portfolio_value
        history[0, 18] = np.nan
        history[0, 19] = net_hedged_value

    # ------------------------------------------------------------------
    # Hedge path
    # ------------------------------------------------------------------
    for step in range(
        1,
        num_steps + 1,
    ):

        time = step * dt

        remaining_time = max(
            time_to_expiry - time,
            0.0,
        )

        # --------------------------------------------------------------
        # Cash accrues interest
        # --------------------------------------------------------------
        bank_before_interest = bank_balance

        bank_balance *= growth_factor

        interest_accrued = bank_balance - bank_before_interest

        # --------------------------------------------------------------
        # Dividend income on stock held over the interval
        # --------------------------------------------------------------

        dividend_income = shares * stock_price * (np.exp(dividend_yield * dt) - 1.0)

        bank_balance += dividend_income

        # --------------------------------------------------------------
        # Stock evolution
        # --------------------------------------------------------------
        normal_draw = normal_draws[step - 1]

        stock_price *= np.exp(
            (stock_drift - dividend_yield - 0.5 * realized_volatility**2) * dt
            + realized_volatility * sqrt_dt * normal_draw
        )

        shares_before = shares

        # --------------------------------------------------------------
        # Rehedge
        # --------------------------------------------------------------
        if step < num_steps:

            option_value = european_value(
                stock_price,
                remaining_time,
                strike_price,
                risk_free_rate,
                dividend_yield,
                implied_volatility,
                option_type_int,
            )

            option_position_value = option_value * num_options

            target_delta = delta(
                stock_price,
                remaining_time,
                strike_price,
                risk_free_rate,
                dividend_yield,
                implied_volatility,
                option_type_int,
            )

            target_shares = target_delta * num_options

            shares_traded = target_shares - shares_before

            stock_trade_notional = abs(shares_traded) * stock_price

            transaction_cost = transaction_cost_rate * stock_trade_notional

            stock_trade_cash_flow = -shares_traded * stock_price

            bank_balance += stock_trade_cash_flow - transaction_cost

            total_transaction_costs += transaction_cost

            shares = target_shares

            hedge_portfolio_value = bank_balance + shares * stock_price

            net_hedged_value = hedge_portfolio_value - option_position_value

            option_payoff = np.nan

        # --------------------------------------------------------------
        # Maturity
        # --------------------------------------------------------------
        else:

            payoff_per_option = _option_payoff(
                stock_price,
                strike_price,
                option_type_int,
            )

            option_value = payoff_per_option

            option_payoff = payoff_per_option * num_options

            option_position_value = option_payoff

            target_delta = np.nan

            # Liquidate hedge.
            shares_traded = -shares_before

            stock_trade_notional = abs(shares_traded) * stock_price

            transaction_cost = transaction_cost_rate * stock_trade_notional

            stock_trade_cash_flow = -shares_traded * stock_price

            bank_balance += stock_trade_cash_flow - transaction_cost

            total_transaction_costs += transaction_cost

            shares = 0.0

            hedge_portfolio_value = bank_balance

            net_hedged_value = hedge_portfolio_value - option_payoff

        # --------------------------------------------------------------
        # History
        # --------------------------------------------------------------
        if store_history:

            history[step, 0] = step
            history[step, 1] = time
            history[step, 2] = remaining_time
            history[step, 3] = normal_draw
            history[step, 4] = stock_price
            history[step, 5] = option_value
            history[step, 6] = option_position_value
            history[step, 7] = target_delta
            history[step, 8] = shares
            history[step, 9] = shares_traded
            history[step, 10] = stock_trade_notional
            history[step, 11] = stock_trade_cash_flow
            history[step, 12] = dividend_income
            history[step, 13] = transaction_cost
            history[step, 14] = total_transaction_costs
            history[step, 15] = interest_accrued
            history[step, 16] = bank_balance
            history[step, 17] = hedge_portfolio_value
            history[step, 18] = option_payoff
            history[step, 19] = net_hedged_value

    return (
        net_hedged_value,
        stock_price,
        option_payoff,
        total_transaction_costs,
        initial_option_value,
    )


##############################################################################
# Multi-path kernel
##############################################################################


@njit(cache=True)
def _simulate_hedge_paths_kernel(
    normal_draws,
    num_options,
    option_type_int,
    stock_price,
    strike_price,
    risk_free_rate,
    dividend_yield,
    implied_volatility,
    realized_volatility,
    time_to_expiry,
    stock_drift,
    transaction_cost_rate,
):
    """Run many simulations using the single canonical hedge kernel."""

    num_paths = normal_draws.shape[0]

    hedging_errors = np.empty(num_paths)

    terminal_stock_prices = np.empty(num_paths)

    option_payoffs = np.empty(num_paths)

    transaction_costs = np.empty(num_paths)

    # Dummy array because history is disabled.
    history = np.empty((1, 1))

    initial_option_value = 0.0

    for path in range(num_paths):

        (
            hedging_error,
            terminal_stock_price,
            option_payoff,
            total_transaction_cost,
            option_value_0,
        ) = _simulate_single_path_kernel(
            normal_draws[path],
            num_options,
            option_type_int,
            stock_price,
            strike_price,
            risk_free_rate,
            dividend_yield,
            implied_volatility,
            realized_volatility,
            time_to_expiry,
            stock_drift,
            transaction_cost_rate,
            False,
            history,
        )

        hedging_errors[path] = hedging_error

        terminal_stock_prices[path] = terminal_stock_price

        option_payoffs[path] = option_payoff

        transaction_costs[path] = total_transaction_cost

        initial_option_value = option_value_0

    return (
        hedging_errors,
        terminal_stock_prices,
        option_payoffs,
        transaction_costs,
        initial_option_value,
    )


##############################################################################
# Public multi-path API
##############################################################################


def simulate_hedge_paths(
    num_paths: int,
    num_options: int,
    option_type_int: int,
    stock_price: float,
    strike_price: float,
    risk_free_rate: float,
    dividend_yield: float,
    implied_volatility: float,
    realized_volatility: float,
    time_to_expiry: float,
    num_steps: int = 252,
    stock_drift: Optional[float] = None,
    transaction_cost_rate: float = 0.0,
    seed: Optional[int] = None,
) -> HedgeSimulationResult:

    if num_paths <= 0:
        raise FinError("Number of paths must be positive.")

    _validate_inputs(
        num_options,
        option_type_int,
        stock_price,
        strike_price,
        implied_volatility,
        realized_volatility,
        time_to_expiry,
        num_steps,
        transaction_cost_rate,
    )

    if stock_drift is None:
        stock_drift = risk_free_rate

    rng = np.random.default_rng(seed)

    normal_draws = rng.standard_normal(
        (
            num_paths,
            num_steps,
        )
    )

    (
        hedging_errors,
        terminal_stock_prices,
        option_payoffs,
        transaction_costs,
        initial_option_value,
    ) = _simulate_hedge_paths_kernel(
        normal_draws,
        num_options,
        option_type_int,
        stock_price,
        strike_price,
        risk_free_rate,
        dividend_yield,
        implied_volatility,
        realized_volatility,
        time_to_expiry,
        stock_drift,
        transaction_cost_rate,
    )

    return HedgeSimulationResult(
        hedging_errors=hedging_errors,
        terminal_stock_prices=terminal_stock_prices,
        option_payoffs=option_payoffs,
        transaction_costs=transaction_costs,
        initial_option_value=float(initial_option_value),
    )


##############################################################################
# Public single-path API
##############################################################################


def simulate_single_hedge_path(
    num_options: int,
    option_type_int: int,
    stock_price: float,
    strike_price: float,
    risk_free_rate: float,
    dividend_yield: float,
    implied_volatility: float,
    realized_volatility: float,
    time_to_expiry: float,
    num_steps: int = 252,
    stock_drift: Optional[float] = None,
    transaction_cost_rate: float = 0.0,
    rng: Optional[np.random.Generator] = None,
    store_history: bool = True,
) -> HedgePathResult:
    """
    Simulate one delta-hedged European option path.

    The hedge calculation itself is performed entirely by
    _simulate_single_path_kernel().
    """

    _validate_inputs(
        num_options,
        option_type_int,
        stock_price,
        strike_price,
        implied_volatility,
        realized_volatility,
        time_to_expiry,
        num_steps,
        transaction_cost_rate,
    )

    if stock_drift is None:
        stock_drift = risk_free_rate

    if rng is None:
        rng = np.random.default_rng()

    normal_draws = rng.standard_normal(num_steps)

    # Allocate history storage.
    if store_history:

        history_array = np.empty(
            (
                num_steps + 1,
                20,
            ),
            dtype=np.float64,
        )

    else:

        # Dummy array; kernel will not write into it.
        history_array = np.empty(
            (1, 1),
            dtype=np.float64,
        )

    (
        hedging_error,
        terminal_stock_price,
        option_payoff,
        total_transaction_costs,
        initial_option_value,
    ) = _simulate_single_path_kernel(
        normal_draws,
        num_options,
        option_type_int,
        stock_price,
        strike_price,
        risk_free_rate,
        dividend_yield,
        implied_volatility,
        realized_volatility,
        time_to_expiry,
        stock_drift,
        transaction_cost_rate,
        store_history,
        history_array,
    )

    if store_history:

        history = pd.DataFrame(
            history_array,
            columns=[
                "step",
                "time",
                "remaining_time",
                "normal_draw",
                "stock_price",
                "option_value",
                "option_position_value",
                "target_delta",
                "shares_held",
                "shares_traded",
                "stock_trade_notional",
                "stock_trade_cash_flow",
                "dividend_income",
                "transaction_cost",
                "cumulative_transaction_cost",
                "interest_accrued",
                "bank_after_trade",
                "hedge_portfolio_value",
                "option_payoff",
                "net_hedged_value",
            ],
        )

        # Step is conceptually an integer.
        history["step"] = history["step"].astype(int)

    else:

        history = None

    return HedgePathResult(
        history=history,
        hedging_error=float(hedging_error),
        terminal_stock_price=float(terminal_stock_price),
        option_payoff=float(option_payoff),
        initial_option_value=float(initial_option_value),
        total_transaction_costs=float(total_transaction_costs),
    )


#####################################################################################

if __name__ == "__main__":

    rng = np.random.default_rng(42)

    result = simulate_single_hedge_path(
        num_options=100,
        option_type_int=OptionTypes.EUROPEAN_CALL.value,
        stock_price=100.0,
        strike_price=100.0,
        risk_free_rate=0.03,
        dividend_yield=0.0,
        implied_volatility=0.20,
        realized_volatility=0.30,
        time_to_expiry=1.0,
        num_steps=252,
        transaction_cost_rate=0.001,
        rng=rng,
    )

    print(result.history)
    print(result.hedging_error)
    print(result.total_transaction_costs)

    result = simulate_hedge_paths(
        num_paths=10_000,
        num_options=1,
        option_type_int=OptionTypes.EUROPEAN_CALL.value,
        stock_price=100.0,
        strike_price=100.0,
        risk_free_rate=0.03,
        dividend_yield=0.02,
        implied_volatility=0.20,
        realized_volatility=0.20,
        time_to_expiry=1.0,
        num_steps=252,
        transaction_cost_rate=0.00,
        seed=42,
    )

    result.plot_hedging_error_vs_terminal_stock_price()
