##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from numbers import Real

from financepy.utils.error import FinError
from financepy.utils.check_values import check_currency


class Amount:
    """A monetary amount together with its currency."""

    def __init__(
        self,
        amount: float,
        ccy_name: str,
    ) -> None:
        if not isinstance(amount, Real):
            raise FinError("Amount must be numeric.")

        check_currency(ccy_name)
        self.amount = float(amount)
        self.ccy_name = ccy_name

    def __repr__(self):
        return (
            f"{self.ccy_name} "
            f"{self.amount:,.2f}"
        )

    def __float__(self):
        return self.amount

    def __neg__(self):
        return Amount(
            -self.amount,
            self.ccy_name,
        )

    def __abs__(self):
        return Amount(
            abs(self.amount),
            self.ccy_name,
        )

    def __add__(self, other):
        if not isinstance(other, Amount):
            return NotImplemented

        self._check_currency(other)

        return Amount(
            self.amount + other.amount,
            self.ccy_name,
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
            self.ccy_name,
        )

    def __rmul__(self, scalar):
        return self.__mul__(scalar)

    def __truediv__(self, other):
        if isinstance(other, Real):
            return Amount(
                self.amount / other,
                self.ccy_name,
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
        if self.ccy_name != other.ccy_name:
            raise FinError(
                "Cannot combine amounts with different currencies."
            )
