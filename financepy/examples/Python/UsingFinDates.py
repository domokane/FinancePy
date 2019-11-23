# -*- coding: utf-8 -*-
"""
Created on Sun Jan 13 21:52:16 2019

@author: Dominic O'Kane
"""

from financepy.finutils.FinDate import FinDate

# We can creat a FinDate as follows

dt1 = FinDate(2019, 10, 10)

# To see what this is, print it

print("PRINT DATES")
print(dt1)

# To add two days we can do

print("ADD CALENDAR DAYS")
dt2 = dt1.addDays(2)
print(dt2)

# dt has not changed, we just created a new date

# To add business days we do the following - note this does not know 
# about regional or religious holidays - just weekends
print("ADD WORKDAYS")
print(dt1, dt1.addWorkDays(2))

# The weekend has now been skipped
# To add a month do

print("ADD MONTHS")
print(dt1, dt1.addMonths(2))

# An invalid date will throw an error
# dt5 = FinDate(2019, 1, 31)
# print(dt5)

# You can use tenors - a number and a 'd', 'm' or 'y' in upper or lower case
print("TENORS")
print(dt1.addTenor("1d"))
print(dt1.addTenor("2D"))
print(dt1.addTenor("3M"))
print(dt1.addTenor("4m"))
print(dt1.addTenor("5Y"))
print(dt1.addTenor("6y"))

# You can subtract dates
print("SUBTRACT DATES")
dt6 = dt1.addTenor("5Y")
dt7 = dt1.addTenor("10Y")
dd = dt7 - dt6
print(dt6, dt7, dd)

# check if a date is on a weekend
print("WEEKEND TEST")
print(dt1, dt1.isWeekend())
print(dt2, dt2.isWeekend())

# You can get CDS dates now or in the future
print("CDS Date")
print(dt1, dt1.nextCDSDate())

# You can get the next IMM dates now
print("IMM Date")
print(dt1, dt1.nextIMMDate())
