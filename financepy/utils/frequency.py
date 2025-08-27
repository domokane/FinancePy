##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from enum import Enum

from ..utils.error import FinError

########################################################################################


class FrequencyTypes(Enum):
    """Enumeration of frequency types."""

    ZERO = -1
    SIMPLE = 0
    ANNUAL = 1
    SEMI_ANNUAL = 2
    TRI_ANNUAL = 3
    QUARTERLY = 4
    MONTHLY = 12
    CONTINUOUS = 99


########################################################################################


def annual_frequency(freq_type: FrequencyTypes) -> float:
    """
    Return the number of payment periods per year for a given FrequencyType.
    CONTINUOUS returns -1. ZERO returns 1 (to avoid division by zero).
    """
    if not isinstance(freq_type, FrequencyTypes):
        print("FinFrequency:", freq_type)
        raise FinError("Unknown frequency type")

    if freq_type == FrequencyTypes.CONTINUOUS:
        return -1.0
    if freq_type == FrequencyTypes.SIMPLE:  # THE RETURN VALUE IS NEVER USED
        return -99
    if freq_type == FrequencyTypes.ZERO:
        return 1.0
    if freq_type == FrequencyTypes.ANNUAL:
        return 1.0
    if freq_type == FrequencyTypes.SEMI_ANNUAL:
        return 2.0
    if freq_type == FrequencyTypes.TRI_ANNUAL:
        return 3.0
    if freq_type == FrequencyTypes.QUARTERLY:
        return 4.0
    if freq_type == FrequencyTypes.MONTHLY:
        return 12.0

    raise FinError("Invalid frequency type")


########################################################################################
