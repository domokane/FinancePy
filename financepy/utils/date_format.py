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
G_DATE_TYPE_FORMAT: DateFormatTypes = DateFormatTypes.UK_LONG

########################################################################################


def set_date_format(fmt: str) -> None:
    """Set the global date format."""
    global G_DATE_TYPE_FORMAT
    G_DATE_TYPE_FORMAT = fmt


########################################################################################


def get_date_format() -> str:
    """Return the current global date format."""
    return G_DATE_TYPE_FORMAT


########################################################################################


def print_date_format():
    global G_DATE_TYPE_FORMAT
    print("TEST TYPE", G_DATE_TYPE_FORMAT)


########################################################################################
