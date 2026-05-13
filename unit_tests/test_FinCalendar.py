# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.utils.calendar import Calendar, CalendarTypes
from financepy.utils.date import Date


bus_days_in_decade = {
    CalendarTypes.NONE: 3653,
    CalendarTypes.WEEKEND: 2609,
    CalendarTypes.AUSTRALIA: 2517,
    CalendarTypes.CANADA: 2502,
    CalendarTypes.FRANCE: 2507,
    CalendarTypes.GERMANY: 2529,
    CalendarTypes.ITALY: 2519,
    CalendarTypes.JAPAN: 2466,
    CalendarTypes.NEW_ZEALAND: 2520,
    CalendarTypes.NORWAY: 2526,
    CalendarTypes.SWEDEN: 2514,
    CalendarTypes.SWITZERLAND: 2530,
    CalendarTypes.TARGET: 2562,
    CalendarTypes.UNITED_STATES: 2507,
    CalendarTypes.UNITED_KINGDOM: 2527,
}

########################################################################################


def test_add_business_day():
    start = Date(3, 1, 2020)
    end = Date(3, 1, 2030)

    for cal_type in CalendarTypes:

        num_days = bus_days_in_decade[cal_type]
        cal = Calendar(cal_type)

        assert cal.add_business_days(start, num_days) == end

        assert cal.add_business_days(end, -num_days) == start


test_add_business_day()