from enum import Enum

from financepy.utils.error import FinError


class TenorUnit(Enum):
    NONE = 0
    DAYS = 1
    WEEKS = 7
    MONTHS = 7 * 4
    YEARS = 12 * 7 * 4


class Tenor:
    """
    A class to represent a Tenor such as '1D' or '10Y'. There is a unit
    (such as a day or a month) and a number representing how many units are in
    the tenor. For example 10Y means 10 years while 120m also means 10 years

    We support limited conversion from bigger units to smaller units, which
    also uses pretty basic rules such as
    1Y = 12M = 12*4 W = 12*4*7 D
    so use conversion at your own risk. Specifically Y->M is ok and W->D is
    ok but Y->D is clearly not that great

    We also support limited arithmetic operations that make Tenors a linear
    space over integers (ie you can add tenors
    and multiply them by integers )
    """

    def __init__(self, tenor_string: str = None):
        self._units = TenorUnit.NONE
        self._num_periods = 0

        if not tenor_string:
            # create an invalid one
            return

        tenor_string = tenor_string.upper()

        if (
            tenor_string == "ON"
        ):  # overnight - should be used only if spot days = 0
            self._units = TenorUnit.DAYS
            self._num_periods = 1
        elif (
            tenor_string == "TN"
        ):  # overnight - should be used when spot days > 0
            self._units = TenorUnit.DAYS
            self._num_periods = 1
        elif tenor_string[-1] == "D":
            self._units = TenorUnit.DAYS
            self._num_periods = int(tenor_string[:-1])
        elif tenor_string[-1] == "W":
            self._units = TenorUnit.WEEKS
            self._num_periods = int(tenor_string[:-1])
        elif tenor_string[-1] == "M":
            self._units = TenorUnit.MONTHS
            self._num_periods = int(tenor_string[:-1])
        elif tenor_string[-1] == "Y":
            self._units = TenorUnit.YEARS
            self._num_periods = int(tenor_string[:-1])
        else:
            raise FinError("Unknown tenor type in " + tenor_string)

    @classmethod
    def as_tenor(cls, str_or_tenor):
        if isinstance(str_or_tenor, Tenor):
            return str_or_tenor
        return Tenor(str_or_tenor)

    def is_valid(self):
        return self._units != TenorUnit.NONE

    def __mul__(self, other: int):
        """Multiply a tenor by an integer

        Args:
            other (int): Multiplier

        Returns:
            The result of multiplication
        """
        res = Tenor()
        res._units = self._units
        res._num_periods = self._num_periods * other
        return res

    # right multiply is the same as left multiply
    __rmul__ = __mul__

    def __add__(self, other):
        """Adds two tenors. Use the more fine-grained unit of the two for
        the result

        Args:
            other (Tenor): Tenor to add

        Returns:
            The sum of two Tenors as a Tenor
        """

        t1, t2 = Tenor.convert_to_same_units(self, other)

        common_unit = t1._units
        res = Tenor()
        res._units = common_unit
        res._num_periods = t1._num_periods + t2._num_periods
        return res

    def convert_to_unit(self, new_units: TenorUnit):
        """Convert a tenor to an equivalent tenor in new units. Only
        "downcasting" is allowed ie W->D ok, D->W not ok. This can obviously
        be improved. Very nominal conversations used, see the code for details

        Args:
            new_unit (TenorUnit): New unit for the tenor

        Raises:
            FinError: If units cannot be converted

        Returns:
            The same tenor in new units
        """

        if self._units == new_units:
            return self

        if self._units.value < new_units.value:
            raise FinError(
                f"Cannot convert units from {self._units.name}"
                + "to {new_units.name}"
            )

        res = Tenor()
        res._units = new_units
        res._num_periods = (
            self._num_periods * self._units.value / new_units.value
        )
        return res

    @classmethod
    def convert_to_same_units(cls, t1_in, t2_in):
        (t1, t2) = (
            (t1_in, t2_in)
            if t1_in._units.value <= t2_in._units.value
            else (t2_in, t1_in)
        )
        common_unit = t1._units

        t1 = t1.convert_to_unit(common_unit)
        t2 = t2.convert_to_unit(common_unit)
        return t1, t2

    def __eq__(self, other):
        """Equality in the strict sense, ie the units and the nummber of
        periods must be the same. In particular 12M != 1Y
        Perhaps this should be relaxed, but at the moment the idea is that a
        more relaxed equality could be explictly called
        if the user first calls Tenor.convert_to_same_units method

        Args:
            other (Tenor): rhs of equality

        Returns:
            True or False depending if teh two tenors are equal
        """

        return (
            self._units == other._units
            and self._num_periods == other._num_periods
        )

    def __repr__(self):
        if self.is_valid():
            return f"{self._num_periods}{self._units.name[0]}"
        else:
            return "*Invalid*"
