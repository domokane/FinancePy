##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from .helpers import check_argument_types
from .currency import CurrencyTypes
from .math import ONE_MILLION

########################################################################################


class Amount:
    """An Amount is a holder for an amount in a specific currency."""

    def __init__(
        self,
        notional: float = ONE_MILLION,
        currency_type: CurrencyTypes = CurrencyTypes.NONE,
    ):
        """Create Amount object."""

        check_argument_types(self.__init__, locals())

        self.notional = notional
        self.currency_type = currency_type

    def __repr__(self):
        """Print out amount object."""

        s = ""
        if self.currency_type != CurrencyTypes.NONE:
            s += self.currency_type.name
            s += " "

        s += "{:,.2f}".format(self.notional)

        return s

    def amount(self):
        """Get the notional amount."""
        return self.notional

    def _print(self):
        """Print out the details of the schedule and the actual dates. This
        can be used for providing transparency on schedule calculations."""
        print(self)


########################################################################################
