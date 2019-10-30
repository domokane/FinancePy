# -*- coding: utf-8 -*-
"""
Created on Sat Feb 06 07:26:46 2016

@author: Dominic O'Kane
"""

import sys
sys.path.append("..//..")

from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode

testCases = FinTestCases(__file__,globalTestCaseMode)

def test_FinDate():

    startDate = FinDate(2018,1,1)

    testCases.header("DATE","MONTHS","CDS DATE") 

    for numMonths in range(0,120):
        nextCDSDate = startDate.nextCDSDate(numMonths)
        testCases.print(str(startDate),numMonths,str(nextCDSDate))
        
    startDate = FinDate(2018,1,1)

    testCases.header("STARTDATE","MONTHS","CDS DATE") 

    for numMonths in range(0,365):
        startDate = startDate.addDays(1)
        nextIMMDate = startDate.nextIMMDate()    
        testCases.print(numMonths,str(startDate),str(nextIMMDate))

test_FinDate()
testCases.compareTestCases()