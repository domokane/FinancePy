# -*- coding: utf-8 -*-
"""
Created on Sun Jan 13 21:52:16 2019

@author: Dominic O'Kane
"""

from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinCalendar import FinCalendar, FinCalendarTypes
from financepy.finutils.FinCalendar import FinBusDayAdjustTypes

# Specify US holiday calendar
cal_us = FinCalendar(FinCalendarTypes.US)

print("Business Day is not changed")
dt = FinDate(2020, 7, 8)
print(dt, cal_us._type, cal_us.isBusinessDay(dt))
nextDate = cal_us.adjust(dt, FinBusDayAdjustTypes.FOLLOWING)
print(dt, nextDate)
print("")

print("US Independence Day is changed forwards")
dt = FinDate(2020, 7, 4)
adj = FinBusDayAdjustTypes.FOLLOWING
print(dt, cal_us._type, cal_us.isBusinessDay(dt))
nextDate = cal_us.adjust(dt, adj)
print(dt, adj, nextDate)
print("")

print("US Independence Day is changed back - to Thursday as 3rd is Friday!")
dt = FinDate(2020, 7, 4)
adj = FinBusDayAdjustTypes.PRECEDING
print(dt, cal_us._type, cal_us.isBusinessDay(dt))
nextDate = cal_us.adjust(dt, adj)
print(dt, adj, nextDate)
