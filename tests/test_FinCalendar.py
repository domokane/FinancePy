###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.calendar import Calendar, CalendarTypes
from financepy.utils.date import set_date_format, DateFormatTypes
from financepy.utils.date import Date
import sys

# Between 3rd of January 2020 and 3rd of January 2030
bus_days_in_decade = {
    "CalendarTypes.NONE": 2609,
    "CalendarTypes.WEEKEND": 2609,
    "CalendarTypes.AUSTRALIA": 2517,
    "CalendarTypes.CANADA": 2502,
    "CalendarTypes.FRANCE": 2507,
    "CalendarTypes.GERMANY": 2529,
    "CalendarTypes.ITALY": 2519,
    "CalendarTypes.JAPAN": 2466,
    "CalendarTypes.NEW_ZEALAND": 2520,
    "CalendarTypes.NORWAY": 2526,
    "CalendarTypes.SWEDEN": 2514,
    "CalendarTypes.SWITZERLAND": 2530,
    "CalendarTypes.TARGET": 2562,
    "CalendarTypes.UNITED_STATES": 2507,
    "CalendarTypes.UNITED_KINGDOM": 2527
}


def test_add_business_day():
    for calendar_type in CalendarTypes:
        num_days = bus_days_in_decade[str(calendar_type)]
        cal = Calendar(calendar_type)
        start = Date(3, 1, 2020)
        end = Date(3, 1, 2030)

        assert cal.add_business_days(start, num_days) == end, \
            f"Landed on incorrect business day using {calendar_type}"
