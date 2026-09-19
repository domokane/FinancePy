##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from numbers import Real

from financepy.utils.error import FinError
from financepy.utils.currency import CurrencyTypes


class Amount:
    """A monetary amount together with its currency."""

    def __init__(
        self,
        amount: float,
        currency_type: CurrencyTypes,
    ):
        if not isinstance(amount, Real):
            raise FinError("Amount must be numeric.")

        if not isinstance(currency_type, CurrencyTypes):
            raise FinError("Invalid currency type.")

        self.amount = float(amount)
        self.currency_type = currency_type

    def __repr__(self):
        return (
            f"{self.currency_type.name} "
            f"{self.amount:,.2f}"
        )

    def __float__(self):
        return self.amount

    def __neg__(self):
        return Amount(
            -self.amount,
            self.currency_type,
        )

    def __abs__(self):
        return Amount(
            abs(self.amount),
            self.currency_type,
        )

    def __add__(self, other):
        if not isinstance(other, Amount):
            return NotImplemented

        self._check_currency(other)

        return Amount(
            self.amount + other.amount,
            self.currency_type,
        )

    def __sub__(self, other):
        if not isinstance(other, Amount):
            return NotImplemented

        self._check_currency(other)

        return Amount(
            self.amount - other.amount,
            self.currency_type,
        )

    def __mul__(self, scalar):
        if not isinstance(scalar, Real):
            return NotImplemented

        return Amount(
            self.amount * scalar,
            self.currency_type,
        )

    def __rmul__(self, scalar):
        return self.__mul__(scalar)

    def __truediv__(self, other):
        if isinstance(other, Real):
            return Amount(
                self.amount / other,
                self.currency_type,
            )

        if isinstance(other, Amount):
            self._check_currency(other)
            return self.amount / other.amount

        return NotImplemented

    def __eq__(self, other):
        if not isinstance(other, Amount):
            return False

        return (
            self.currency_type == other.currency_type
            and self.amount == other.amount
        )

    def _check_currency(self, other):
        if self.currency_type != other.currency_type:
            raise FinError(
                "Cannot combine amounts with different currencies."
            )