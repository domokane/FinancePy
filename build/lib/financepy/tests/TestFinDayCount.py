# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 10:10:19 2018

@author: Dominic O'Kane
"""

from financepy.finutils.FinDayCount import FinDayCount, FinDayCountTypes
from financepy.finutils.FinDate import FinDate

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__,globalTestCaseMode)

def test_FinDayCount():

    testCases.header("DAY_COUNT_METHOD","START","END","ALPHA")

    for dayCountMethod in FinDayCountTypes:

        startDate = FinDate(2019,1,1)
        nextDate = startDate
        numDays = 20
        dayCount = FinDayCount(dayCountMethod)
    
        for i in range(0,numDays):
            nextDate = nextDate.addDays(7)
            alpha = dayCount.yearFrac(startDate,nextDate)
            testCases.print(str(dayCountMethod),str(startDate),str(nextDate),alpha)

test_FinDayCount()
testCases.compareTestCases()




