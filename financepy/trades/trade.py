##############################################################################
# Copyright (C) 2018-2026 Dominic O'Kane
##############################################################################

# WORK IN PROGRESS - DO NOT USE

from typing import Any, Optional

from ..products.product import Product
from ..utils.date import Date
from ..utils.error import FinError


class Trade:
    """Represents a transaction in a financial product.

    A Trade holds transaction-specific information, while
    the Product is responsible for pricing and risk.

    Attributes
    ----------
    product : Product
        Financial product associated with the trade.

    initial_trade_dt : Date
        Date on which the trade was executed.

    trade_id : str, optional
        Unique trade identifier.

    counterparty : str, optional
        Counterparty to the transaction.

    quantity : float
        Number of units of the product.

    """

    def __init__(
        self,
        product: Product,
        initial_trade_dt: Date,
        trade_id: Optional[str] = None,
        counterparty: Optional[str] = None,
        quantity: float = 1.0,
    ) -> None:

        if not isinstance(product, Product):
            raise FinError("Product must inherit from Product.")

        if not isinstance(initial_trade_dt, Date):
            raise FinError("Initial trade date must be a Date.")

        if not isinstance(quantity, (int, float)):
            raise FinError("Quantity must be numeric.")

        if quantity == 0.0:
            raise FinError("Quantity cannot be zero.")

        self.product = product
        self.initial_trade_dt = initial_trade_dt
        self.trade_id = trade_id
        self.counterparty = counterparty
        self.quantity = float(quantity)

    ###########################################################################
    # VALUATION
    ###########################################################################

    def value(self, *args, **kwargs) -> Any:
        """Calculate the position value of the trade.

        All valuation arguments are forwarded to the
        underlying product.

        A positive quantity represents a long position.
        A negative quantity represents a short position.
        """

        product_value = self.product.value(*args, **kwargs)

        return self._scale_result(product_value)

    ###########################################################################
    # RISK
    ###########################################################################

    def delta(self, *args, **kwargs) -> Any:
        """Calculate the position delta."""

        return self._risk("delta", *args, **kwargs)

    def gamma(self, *args, **kwargs) -> Any:
        """Calculate the position gamma."""

        return self._risk("gamma", *args, **kwargs)

    def vega(self, *args, **kwargs) -> Any:
        """Calculate the position vega."""

        return self._risk("vega", *args, **kwargs)

    def theta(self, *args, **kwargs) -> Any:
        """Calculate the position theta."""

        return self._risk("theta", *args, **kwargs)

    def rho(self, *args, **kwargs) -> Any:
        """Calculate the position rho."""

        return self._risk("rho", *args, **kwargs)

    def _risk(self, method_name: str, *args, **kwargs) -> Any:
        """Delegate a risk calculation to the product."""

        method = getattr(self.product, method_name, None)

        if not callable(method):
            raise FinError(
                f"{type(self.product).__name__} "
                f"does not implement {method_name}."
            )

        result = method(*args, **kwargs)

        return self._scale_result(result)

    ###########################################################################
    # INTERNAL HELPERS
    ###########################################################################

    def _scale_result(self, result: Any) -> Any:
        """Scale product results by the trade quantity.

        Supports numeric results and dictionaries containing
        numeric valuation or risk results.

        Non-numeric metadata is preserved unchanged.
        """

        if isinstance(result, dict):
            return {
                key: self._scale_result(value)
                if isinstance(value, (int, float, dict))
                else value
                for key, value in result.items()
            }

        if isinstance(result, (int, float)):
            return result * self.quantity

        raise FinError(
            "Unsupported valuation or risk result type."
        )

    ###########################################################################
    # REPRESENTATION
    ###########################################################################

    def __repr__(self) -> str:
        """Return the trade representation."""

        return (
            f"Trade(\n"
            f"  trade_id={self.trade_id!r},\n"
            f"  initial_trade_dt={self.initial_trade_dt!r},\n"
            f"  counterparty={self.counterparty!r},\n"
            f"  quantity={self.quantity},\n"
            f"  product={self.product!r}\n"
            f")"
        )

    def _print(self) -> None:
        """Print trade details."""
        print(self)
