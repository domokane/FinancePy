##############################################################################
# Copyright (C) 2018-2026 Dominic O'Kane
##############################################################################

from abc import ABC, abstractmethod
from typing import Any

from financepy.utils.date import Date


class Product(ABC):
    """Abstract base class for financial products.

    Defines a common interface for valuation and
    representation of financial instruments.
    """

    @abstractmethod
    def value(self, value_dt: Date, *args, **kwargs) -> Any:
        """Calculate the value of the financial product."""
        pass

    @abstractmethod
    def __repr__(self) -> str:
        """Return a string representation of the product."""
        pass

    def _print(self) -> None:
        """Print the product representation."""
        print(self)
