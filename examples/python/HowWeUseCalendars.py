# -*- coding: utf-8 -*-
"""
Created on Sun Jan 13 21:52:16 2019

@author: Dominic O'Kane
"""

from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinCalendar import FinCalendar, FinCalendarTypes


cal_none = FinCalendar(FinCalendarTypes.NONE)
cal_weekend = FinCalendar(FinCalendarTypes.WEEKEND)
cal_uk = FinCalendar(FinCalendarTypes.UK)
cal_us = FinCalendar(FinCalendarTypes.US)
cal_target = FinCalendar(FinCalendarTypes.TARGET)

print("TEST WHETHER THESE DATES ARE BUSINESS DATES")

# Let us check some dates - we can create a weekend
print("Weekend")
dt = FinDate(2019, 1, 5)
print(dt, cal_none._type, cal_none.isBusinessDay(dt))
print(dt, cal_weekend._type, cal_weekend.isBusinessDay(dt))
print(dt, cal_uk._type, cal_uk.isBusinessDay(dt))
print(dt, cal_us._type, cal_us.isBusinessDay(dt))
print(dt, cal_target._type, cal_target.isBusinessDay(dt))

print("New Year's Day")
dt = FinDate(2019, 1, 1)
print(dt, cal_none._type, cal_none.isBusinessDay(dt))
print(dt, cal_weekend._type, cal_weekend.isBusinessDay(dt))
print(dt, cal_uk._type, cal_uk.isBusinessDay(dt))
print(dt, cal_us._type, cal_us.isBusinessDay(dt))
print(dt, cal_target._type, cal_target.isBusinessDay(dt))

print("Christmas Day")
dt = FinDate(2019, 12, 25)
print(dt, cal_none._type, cal_none.isBusinessDay(dt))
print(dt, cal_weekend._type, cal_weekend.isBusinessDay(dt))
print(dt, cal_uk._type, cal_uk.isBusinessDay(dt))
print(dt, cal_us._type, cal_us.isBusinessDay(dt))
print(dt, cal_target._type, cal_target.isBusinessDay(dt))

print("May Day")
dt = FinDate(2019, 5, 1)
print(dt, cal_none._type, cal_none.isBusinessDay(dt))
print(dt, cal_weekend._type, cal_weekend.isBusinessDay(dt))
print(dt, cal_uk._type, cal_uk.isBusinessDay(dt))
print(dt, cal_us._type, cal_us.isBusinessDay(dt))
print(dt, cal_target._type, cal_target.isBusinessDay(dt))

print("US Independence Day")
dt = FinDate(2019, 7, 4)
print(dt, cal_none._type, cal_none.isBusinessDay(dt))
print(dt, cal_weekend._type, cal_weekend.isBusinessDay(dt))
print(dt, cal_uk._type, cal_uk.isBusinessDay(dt))
print(dt, cal_us._type, cal_us.isBusinessDay(dt))
print(dt, cal_target._type, cal_target.isBusinessDay(dt))

print("Bastille Day")
dt = FinDate(2020, 7, 14)
print(dt, cal_none._type, cal_none.isBusinessDay(dt))
print(dt, cal_weekend._type, cal_weekend.isBusinessDay(dt))
print(dt, cal_uk._type, cal_uk.isBusinessDay(dt))
print(dt, cal_us._type, cal_us.isBusinessDay(dt))
print(dt, cal_target._type, cal_target.isBusinessDay(dt))

print("Veterans Day")
dt = FinDate(2020, 11, 11)
print(dt, cal_none._type, cal_none.isBusinessDay(dt))
print(dt, cal_weekend._type, cal_weekend.isBusinessDay(dt))
print(dt, cal_uk._type, cal_uk.isBusinessDay(dt))
print(dt, cal_us._type, cal_us.isBusinessDay(dt))
print(dt, cal_target._type, cal_target.isBusinessDay(dt))
