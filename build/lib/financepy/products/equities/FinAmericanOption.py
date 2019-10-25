# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 10:04:57 2019

@author: Dominic
"""

import numpy as np
from numba import njit, float64, int64

from ...finutils.FinError import FinError
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinDate import FinDate
from ...products.equities.FinVanillaOption import FinOptionTypes

###############################################################################
# TODO
#   Add term structure of rates
###############################################################################

@njit(float64[:](float64,float64,float64,float64,int64,float64,int64,float64),
      fastmath=True,cache=True)
def valueOnce(stockPrice,
              riskFreeRate,
              dividendYield,
              volatility,
              numSteps,
              timeToExpiry,
              optionType,
              strikePrice):
        ''' Value an American option using a Binomial Treee '''
        if numSteps < 3:
            numSteps = 3

        # this is the size of the step
        dt = timeToExpiry/numSteps
        r = riskFreeRate
        q = dividendYield

        # the number of nodes on the tree
        numNodes = int(0.5 * (numSteps+1) * (numSteps+2))
        stockValues = np.zeros(numNodes)
        stockValues[0] = stockPrice

        optionValues = np.zeros(numNodes)
        u = np.exp(volatility * np.sqrt(dt))
        d = 1.0/u
        sLow = stockPrice

        probs = np.zeros(numSteps)
        periodDiscountFactors = np.zeros(numSteps)

        # store time independent information for later use in tree
        for iTime in range(0,numSteps):
            a = np.exp((r-q)*dt)
            probs[iTime] = (a - d)/(u - d)
            periodDiscountFactors[iTime] = np.exp(-r*dt)

        for iTime in range(1,numSteps+1):
            sLow *= d
            s = sLow
            for iNode in range(0,iTime+1):
                index = 0.5 * iTime * (iTime+1)
                stockValues[int(index+iNode)] = s
                s = s * (u*u)
        
        # work backwards by first setting values at expiry date
        index = int(0.5 * numSteps * (numSteps+1))

        for iNode in range(0, iTime+1):

            s = stockValues[index + iNode]

            if optionType == FinOptionTypes.EUROPEAN_CALL.value:
                optionValues[index + iNode] = max(s - strikePrice,0)
            elif optionType == FinOptionTypes.EUROPEAN_PUT.value:
                optionValues[index + iNode] = max(strikePrice -s ,0)
            elif optionType == FinOptionTypes.AMERICAN_CALL.value:
                optionValues[index + iNode] = max(s - strikePrice,0)
            elif optionType == FinOptionTypes.AMERICAN_PUT.value:
                optionValues[index + iNode] = max(strikePrice -s ,0)

        # begin backward steps from expiry to value date
        for iTime in range(numSteps-1,-1,-1):

            index = int(0.5 * iTime * (iTime+1))

            for iNode in range(0,iTime+1):

                s = stockValues[index + iNode]

                if optionType == FinOptionTypes.EUROPEAN_CALL.value:
                    exerciseValue = 0.0
                elif optionType == FinOptionTypes.EUROPEAN_PUT.value:
                    exerciseValue = 0.0
                elif optionType == FinOptionTypes.AMERICAN_CALL.value:
                    exerciseValue = max(s-strikePrice,0.0)
                elif optionType == FinOptionTypes.AMERICAN_PUT.value:
                    exerciseValue = max(strikePrice-s,0.0)

                nextIndex = int(0.5 * (iTime+1) * (iTime + 2))

                nextNodeDn = nextIndex + iNode
                nextNodeUp = nextIndex + iNode + 1

                vUp = optionValues[nextNodeUp]
                vDn = optionValues[nextNodeDn]
                futureExpectedValue = probs[iTime] * vUp
                futureExpectedValue += (1.0 - probs[iTime]) * vDn
                holdValue = periodDiscountFactors[iTime] * futureExpectedValue
                
                if optionType == FinOptionTypes.EUROPEAN_CALL.value:
                    optionValues[index + iNode] = holdValue
                elif optionType == FinOptionTypes.EUROPEAN_PUT.value:
                    optionValues[index + iNode] = holdValue
                elif optionType == FinOptionTypes.AMERICAN_CALL.value:
                    optionValues[index + iNode] = max(exerciseValue, holdValue)
                elif optionType == FinOptionTypes.AMERICAN_PUT.value:
                    optionValues[index + iNode] = max(exerciseValue, holdValue)
                    
        # We calculate all of the important Greeks in one go
        price = optionValues[0]
        delta = (optionValues[2]-optionValues[1])/(stockValues[2] - stockValues[1])
        deltaUp = (optionValues[5]-optionValues[4])/(stockValues[5]-stockValues[4])
        deltaDn = (optionValues[4]-optionValues[3])/(stockValues[4]-stockValues[3])
        gamma = (deltaUp - deltaDn)/(stockValues[2] - stockValues[1])
        theta = (optionValues[4] - optionValues[0])/(2.0 * dt)
        results = np.array([price,delta,gamma,theta])
        return results

###############################################################################

class FinAmericanOption():

    ''' Class that performs the valuation of an American style option on a 
    dividend paying stock. Can easily be extended to price American style FX 
    options. The dividend is assumed to be continuous. '''

    def __init__(self,
             expiryDate,
             strikePrice,
             optionType):

        if strikePrice < 0.0:
            raise FinError("Strike price must be positive.")

        self._expiryDate = expiryDate
        self._strikePrice = float(strikePrice)
        self._optionType = optionType

    def value(self,
              valueDate,
              stockPrice,
              riskFreeRate,
              dividendYield,
              volatility,
              numSteps = 100):

            if valueDate > self._expiryDate:
                raise FinError("Value date is after expiry date.")

            # do some validation
            timeToExpiry = FinDate.datediff(valueDate,self._expiryDate)/gDaysInYear

            # Convergence is improved if I take the average of an even and 
            # odd number of steps in the Binomial Tree

            value1 = valueOnce(stockPrice,
                               riskFreeRate,
                               dividendYield,
                               volatility,
                               numSteps,
                               timeToExpiry,
                               self._optionType.value,
                               self._strikePrice)

            value2 = valueOnce(stockPrice,
                                riskFreeRate,
                                dividendYield,
                                volatility,
                                numSteps+1,
                                timeToExpiry,
                                self._optionType.value,
                                self._strikePrice)

            v = (value1 + value2)/2.0
            return v

###############################################################################