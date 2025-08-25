##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from enum import Enum

########################################################################################


class DateFormatTypes(Enum):
    """Enumeration of date format types."""

    BLOOMBERG = 1
    US_SHORT = 2
    US_MEDIUM = 3
    US_LONG = 4
    US_LONGEST = 5
    UK_SHORT = 6
    UK_MEDIUM = 7
    UK_LONG = 8
    UK_LONGEST = 9
    DATETIME = 10


########################################################################################

# Private module-level variable to hold current global date format
_date_format: DateFormatTypes = DateFormatTypes.UK_LONG

########################################################################################


def set_date_format(format_type: DateFormatTypes):
    """
    Set the current global date format type.

    Parameters
    ----------
    format_type : DateFormatTypes
        The new date format to set globally.
    """
    global _date_format  # pylint: disable=W0603
    _date_format = format_type


########################################################################################


def get_date_format() -> DateFormatTypes:
    """
    Get the current global date format type.

    Returns
    -------
    DateFormatTypes
        The currently set global date format.
    """
    return _date_format
