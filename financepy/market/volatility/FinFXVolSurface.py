##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from scipy.optimize import fmin_powell
from scipy.optimize import newton
import matplotlib.pyplot as plt
from numba import njit, float64

from ...finutils.FinError import FinError
from ...finutils.FinDate import FinDate
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinGlobalTypes import FinOptionTypes
from ...finutils.FinDistribution import FinDistribution
from ...products.fx.FinFXVanillaOption import FinFXVanillaOption
from ...models.FinModelOptionImpliedDbn import optionImpliedDbn
from ...products.fx.FinFXMktConventions import FinFXATMMethod
from ...products.fx.FinFXMktConventions import FinFXDeltaMethod
from ...finutils.FinHelperFunctions import checkArgumentTypes, labelToString
from ...market.curves.FinDiscountCurve import FinDiscountCurve

from ...products.fx.FinFXModelTypes import FinFXModelBlackScholes
from ...market.volatility.FinOptionVolatilityFns import volFunctionClarke

from ...finutils.FinMath import N

###############################################################################
# TODO: Speed up search for strike by providing derivative function to go with
#       delta fit.
###############################################################################

def g2(K, *args):
    ''' This is the objective function used in the determination of the FX
    Option implied strike which is computed in the class below. '''

    self = args[0]
    valueDate = args[1]
    stockPrice = args[2]
    domDF = args[3]
    forDF = args[4]
    delta = args[5]
    deltaType = args[6]
    volatility = args[7]

    self._strikeFXRate = K

    deltaDict = self.fastDelta(valueDate,
                               stockPrice,
                               domDF,
                               forDF,
                               volatility)

    deltaOut = deltaDict[deltaType]
    objFn = delta - deltaOut
    return objFn

###############################################################################

def obj(cvec, *args):
    ''' Return a function that is minimised when the ATM, MS and RR vols have
    been best fitted using the parametric volatility curve represented by cvec
    '''

    self = args[0]
    i = self._tenorIndex
    self._parameters[i, :] = cvec

    call = args[1]
    put = args[2]

    # We first need to solve for the strikes at the 25 delta points using the
    # new volatility curve

    # Match the at-the-money option volatility
    atmCurveVol = self.volFunction(self._K_ATM[i], i)
    targetATMVol = self._atmVols[i]
    term1 = (targetATMVol - atmCurveVol)**2

    # Match the market strangle value but this has to be at the MS strikes
    call._strikeFXRate = self._K_25_D_C_MS[i]
    sigma_K_25_D_C_MS = self.volFunction(self._K_25_D_C_MS[i], i)
    model = FinFXModelBlackScholes(sigma_K_25_D_C_MS)
    V_C_25_MS = call.value(self._valueDate,
                           self._spotFXRate,
                           self._domDiscountCurve,
                           self._forDiscountCurve,
                           model)['v']

    put._strikeFXRate = self._K_25_D_P_MS[i]
    sigma_K_25_D_P_MS = self.volFunction(self._K_25_D_P_MS[i], i)
    model = FinFXModelBlackScholes(sigma_K_25_D_P_MS)
    V_P_25_MS = put.value(self._valueDate,
                          self._spotFXRate,
                          self._domDiscountCurve,
                          self._forDiscountCurve,
                          model)['v']

    V_25_D_MS = V_C_25_MS + V_P_25_MS
    target25StrangleValue = self._V_25_D_MS[i]
    term2 = (V_25_D_MS - target25StrangleValue)**2

    # Match the risk reversal volatility
    K_25_D_C = self.solveForSmileStrike(call, 0.2500, i, self._spotFXRate)
    sigma_K_25_D_C = self.volFunction(K_25_D_C, i)

    K_25_D_P = self.solveForSmileStrike(put, -0.2500, i, self._spotFXRate)
    sigma_K_25_D_P = self.volFunction(K_25_D_P, i)

    targetRRVol = self._riskReversal25DeltaVols[i]
    sigma_25_D_RR = (sigma_K_25_D_C - sigma_K_25_D_P)
    term3 = (sigma_25_D_RR - targetRRVol)**2

    tot = term1 + term2 + term3

#    print(cvec, tot)
    return tot

###############################################################################

# @njit(float64(float64, float64, float64, float64, float64, float64),
#            fastmath=True, cache=True)
# def volFunctionNumba(c0, c1, c2, F, K, texp):
#     x = np.log(F / K)
#     sigma0 = np.exp(c0)
#     arg = x / (sigma0 * np.sqrt(texp))
#     deltax = N(arg) - 0.50  # The -0.50 seems to be missing in book
#     f = c0 + c1 * deltax + c2 * (deltax * deltax)
#     return np.exp(f)

###############################################################################

def deltaFit(K, *args):
    ''' This is the objective function used in the determination of the FX
    Option implied strike which is computed in the class below. '''

    self = args[0]
    option = args[1]
    deltaTarget = args[2]
    tenorIndex = args[3]

    volatility = self.volFunction(K, tenorIndex)
    model = FinFXModelBlackScholes(volatility)

    option._strikeFXRate = K
    deltaDict = option.delta(self._valueDate,
                             self._spotFXRate,
                             self._domDiscountCurve,
                             self._forDiscountCurve,
                             model)

    deltaOut = deltaDict[self._deltaMethodString]
    objFn = deltaTarget - deltaOut

#    print("DeltaFit", K, objFn, deltaTarget, tenorIndex)

    return objFn

###############################################################################

def solveForStrike(valueDate,
                   vanillaOption,
                   spotFXRate,
                   domDiscountCurve,
                   forDiscountCurve,
                   delta,
                   deltaType,
                   volatility):
    ''' This function determines the implied strike of an FX option
    given a delta and the other option details. It uses a one-dimensional
    Newton root search algorith to determine the strike that matches an
    input volatility. '''

#    argtuple = (vanillaOption, valueDate, spotFXRate, domDiscountCurve,
#                forDiscountCurve, delta, deltaType, volatility)

#    sigma = optimize.newton(g, x0=spotFXRate, args=argtuple,
#                            tol=1e-5, maxiter=100, fprime2=None)

    # TODO: SHOULD I SHIFT THE VALUE DATE
    tdel = (vanillaOption._deliveryDate - valueDate) / gDaysInYear
    domDF = domDiscountCurve._df(tdel)
    forDF = forDiscountCurve._df(tdel)

    argtuple2 = (vanillaOption, valueDate, spotFXRate, domDF,
                forDF, delta, deltaType, volatility)

    sigma2 = newton(g2, x0=spotFXRate, args=argtuple2,
                    tol=1e-5, maxiter=100, fprime2=None)

#    print("==> Solve for strike", spotFXRate, sigma2)
    return sigma2

###############################################################################


class FinFXVolSurface():
    ''' Class to hold information characterising an FX volatility surface. '''

    def __init__(self,
                 valueDate: FinDate,
                 spotFXRate: float,
                 currencyPair: str,
                 notionalCurrency,
                 domDiscountCurve: FinDiscountCurve,
                 forDiscountCurve: FinDiscountCurve,
                 tenors,
                 atmVols,
                 mktStrangle25DeltaVols,
                 riskReversal25DeltaVols,
                 atmMethod=FinFXATMMethod.FWD_DELTA_NEUTRAL,
                 deltaMethod=FinFXDeltaMethod.SPOT_DELTA):

        checkArgumentTypes(self.__init__, locals())

        self._valueDate = valueDate
        self._spotFXRate = spotFXRate
        self._currencyPair = currencyPair

        if len(currencyPair) != 6:
            raise FinError("Currency pair must be 6 characters.")

        self._forName = self._currencyPair[0:3]
        self._domName = self._currencyPair[3:6]

        self._notionalCurrency = notionalCurrency
        self._domDiscountCurve = domDiscountCurve
        self._forDiscountCurve = forDiscountCurve
        self._numVolCurves = len(tenors)

        if len(atmVols) != self._numVolCurves:
            raise FinError("Number ATM vols must equal number of tenors")

        if len(atmVols) != self._numVolCurves:
            raise FinError("Number ATM vols must equal number of tenors")

        if len(mktStrangle25DeltaVols) != self._numVolCurves:
            raise FinError("Number MS25D vols must equal number of tenors")

        if len(riskReversal25DeltaVols) != self._numVolCurves:
            raise FinError("Number RR25D vols must equal number of tenors")

        self._tenors = tenors
        self._atmVols = np.array(atmVols)/100.0
        self._mktStrangle25DeltaVols = np.array(mktStrangle25DeltaVols)/100.0
        self._riskReversal25DeltaVols = np.array(riskReversal25DeltaVols)/100.0

        self._atmMethod = atmMethod
        self._deltaMethod = deltaMethod

        if self._deltaMethod == FinFXDeltaMethod.SPOT_DELTA:
            self._deltaMethodString = "pips_spot_delta"
        elif self._deltaMethod == FinFXDeltaMethod.FORWARD_DELTA:
            self._deltaMethodString = "pips_fwd_delta"
        elif self._deltaMethod == FinFXDeltaMethod.SPOT_DELTA_PREM_ADJ:
            self._deltaMethodString = "pct_spot_delta_prem_adj"
        elif self._deltaMethod == FinFXDeltaMethod.FORWARD_DELTA_PREM_ADJ:
            self._deltaMethodString = "pct_fwd_delta_prem_adj"
        else:
            raise FinError("Unknown Delta Type")

        self._tenorIndex = 0

        self._expiryDates = []
        for i in range(0, self._numVolCurves):
            expiryDate = valueDate.addTenor(tenors[i])
            self._expiryDates.append(expiryDate)

        self.buildVolSurface()

###############################################################################

    def volFunction(self, K, tenorIndex):
        ''' Return the volatility for a strike using a given polynomial
        interpolation following Section 3.9 of Iain Clark book. '''

        params = self._parameters[tenorIndex]
        F = self._F0T[tenorIndex]
        texp = self._texp[tenorIndex]

        vol = np.float64(volFunctionClarke(params, F, K, texp))

        if vol.any() < 0.0:
            raise ValueError("Negative volatility. Not permitted.")

        return vol

###############################################################################

    def impliedDbns(self, lowFX, highFX, numIntervals):
        ''' Calculate the pdf for each tenor horizon. Returns a list of 
        FinDistribution objects, one for each tenor horizon. '''

        dbns = []

        for iTenor in range(0, len(self._tenors)):
            
            F = self._F0T[iTenor]
            texp = self._texp[iTenor]

            dFX = (highFX - lowFX)/ numIntervals

            domDF = self._domDiscountCurve._df(texp)
            forDF = self._forDiscountCurve._df(texp)

            rd = -np.log(domDF) / texp
            rf = -np.log(forDF) / texp

            params = self._parameters[iTenor]

            strikes = []
            vols = []

            for iK in range(0, numIntervals):
                strike = lowFX + iK*dFX                
                vol = volFunctionClarke(params, F, strike, texp)
                strikes.append(strike) 
                vols.append(vol)
            
            strikes = np.array(strikes)
            vols = np.array(vols)

            dbn = optionImpliedDbn(self._spotFXRate, texp, rd, rf, strikes, vols)            

            dbns.append(dbn)

        return dbns

###############################################################################

    def buildVolSurface(self):

        S0 = self._spotFXRate
        numVolCurves = self._numVolCurves
        numParameters = 3

        self._F0T = np.zeros(numVolCurves)
        self._K_ATM = np.zeros(numVolCurves)
        self._deltaATM = np.zeros(numVolCurves)

        self._K_25_D_C = np.zeros(numVolCurves)
        self._K_25_D_P = np.zeros(numVolCurves)
        self._K_25_D_C_MS = np.zeros(numVolCurves)
        self._K_25_D_P_MS = np.zeros(numVolCurves)
        self._V_25_D_MS = np.zeros(numVolCurves)

        self._parameters = np.zeros([numVolCurves, numParameters])
        self._texp = np.zeros(numVolCurves)

        atmVol = self._atmVols[0]
        xopt = [np.log(atmVol), 0.050, 0.800]
        c0 = xopt

        for i in range(0, self._numVolCurves):

            self._tenorIndex = i

            expiryDate = self._expiryDates[i]
            texp = (expiryDate - self._valueDate) / gDaysInYear
            self._texp[i] = texp

            domDF = self._domDiscountCurve._df(texp)
            forDF = self._forDiscountCurve._df(texp)

            F0T = S0 * forDF/domDF
            self._F0T[i] = F0T
            atmVol = self._atmVols[i]

            c0 = xopt

            # This follows exposition in Clarke Page 52
            if self._atmMethod == FinFXATMMethod.SPOT:
                self._K_ATM[i] = S0
            elif self._atmMethod == FinFXATMMethod.FWD:
                self._K_ATM[i] = F0T
            elif self._atmMethod == FinFXATMMethod.FWD_DELTA_NEUTRAL:
                self._K_ATM[i] = F0T * np.exp(atmVol*atmVol*texp/2.0)
            elif self._atmMethod == FinFXATMMethod.FWD_DELTA_NEUTRAL_PREM_ADJ:
                self._K_ATM[i] = F0T * np.exp(-atmVol*atmVol*texp/2.0)
            else:
                raise FinError("Unknown Delta Type")

            K_dummy = 999.0

            call = FinFXVanillaOption(expiryDate,
                                      K_dummy,
                                      self._currencyPair,
                                      FinOptionTypes.EUROPEAN_CALL,
                                      1.0,
                                      self._notionalCurrency)

            put = FinFXVanillaOption(expiryDate,
                                     K_dummy,
                                     self._currencyPair,
                                     FinOptionTypes.EUROPEAN_PUT,
                                     1.0,
                                     self._notionalCurrency)

            # Determine the price of a market strangle from market strangle
            # Need to price a call and put that agree with market strangle

            vol25DeltaMS = self._atmVols[i] + self._mktStrangle25DeltaVols[i]

            K_25_D_C_MS = solveForStrike(self._valueDate,
                                         call,
                                         self._spotFXRate,
                                         self._domDiscountCurve,
                                         self._forDiscountCurve,
                                         0.2500,
                                         self._deltaMethodString,
                                         vol25DeltaMS)

            K_25_D_P_MS = solveForStrike(self._valueDate,
                                         put,
                                         self._spotFXRate,
                                         self._domDiscountCurve,
                                         self._forDiscountCurve,
                                         -0.2500,
                                         self._deltaMethodString,
                                         vol25DeltaMS)

            call._strikeFXRate = K_25_D_C_MS
            put._strikeFXRate = K_25_D_P_MS

            # Store the set of strikes in the class
            self._K_25_D_C_MS[i] = K_25_D_C_MS
            self._K_25_D_P_MS[i] = K_25_D_P_MS

            # USE MARKET STRANGLE VOL TO DETERMINE PRICE OF A MARKET STRANGLE
            model = FinFXModelBlackScholes(vol25DeltaMS)

            V_25_C_MS = call.value(self._valueDate,
                                   self._spotFXRate,
                                   self._domDiscountCurve,
                                   self._forDiscountCurve,
                                   model)['v']

            V_25_P_MS = put.value(self._valueDate,
                                  self._spotFXRate,
                                  self._domDiscountCurve,
                                  self._forDiscountCurve,
                                  model)['v']

            # Market price of strangle in the domestic currency
            V_25_D_MS = V_25_C_MS + V_25_P_MS
            self._V_25_D_MS[i] = V_25_D_MS

            if 1 == 0:
                print("Value Date:", self._valueDate)
                print("Expiry Date:", self._expiryDates[i])
                print("T:", self._texp[i])
                print("S0:", self._spotFXRate)
                print("F0T:", self._F0T[i])
                print("ATM Method:", self._atmMethod)
                print("DeltaATM:", self._deltaATM[i])
                print("DeltaMethod:", self._deltaMethod)

            # Determine parameters of vol surface using Powell minimisation
            args = (self, call, put)
            tol = 1e-7
            xopt = fmin_powell(obj, c0, args, disp=False, ftol=tol)
            v = obj(xopt, *args)

            print(i, texp, xopt, v)

            # The function is quadratic. So the quadratic value should be small
            if abs(v) > tol:
                print("Using an xtol:", tol, " maybe make this smaller")
                raise FinError("Failed to fit volatility smile curve.")

            # Calculate the 25 Delta call and put strikes from new vol surface
            if i > 0:
                x0 = self._K_25_D_C[i-1]
            else:
                x0 = self._spotFXRate

            K_25_D_C = self.solveForSmileStrike(call, 0.25, i, x0)

            if i > 0:
                x0 = self._K_25_D_P[i-1]
            else:
                x0 = self._spotFXRate

            K_25_D_P = self.solveForSmileStrike(put, -0.25, i, x0)

            # Store the set of strikes in the class
            self._K_25_D_C[i] = K_25_D_C
            self._K_25_D_P[i] = K_25_D_P

            self._parameters[i, :] = np.array(xopt)

###############################################################################

    def solveForSmileStrike(self,
                            vanillaOption,
                            deltaTarget,
                            tenorIndex, 
                            initialValue):
        ''' Solve for the strike that sets the delta of the option equal to the
        target value of delta allowing the volatility to be a function of the
        strike. '''

        argtuple = (self, vanillaOption, deltaTarget, tenorIndex)

        sigma = newton(deltaFit, x0=initialValue, args=argtuple,
                       rtol=1e-5, maxiter=50, fprime2=None)

        return sigma

###############################################################################

    def checkCalibration(self, checkCalibrationFlag: bool):

        if checkCalibrationFlag:
            print("==========================================================")
            print("====== CHECK CALIBRATION =================================")
            print("==========================================================")
        else:
            return

        K_dummy = 0.0

        for i in range(0, self._numVolCurves):

            expiryDate = self._expiryDates[i]

            call = FinFXVanillaOption(expiryDate,
                                      K_dummy,
                                      self._currencyPair,
                                      FinOptionTypes.EUROPEAN_CALL,
                                      1.0,
                                      self._notionalCurrency, )

            put = FinFXVanillaOption(expiryDate,
                                     K_dummy,
                                     self._currencyPair,
                                     FinOptionTypes.EUROPEAN_PUT,
                                     1.0,
                                     self._notionalCurrency)

            sigma_K_25_D_C = self.volFunction(self._K_25_D_C[i], i)

            model = FinFXModelBlackScholes(sigma_K_25_D_C)
            call._strikeFXRate = self._K_25_D_C[i]
            delta_call = call.delta(self._valueDate,
                                    self._spotFXRate,
                                    self._domDiscountCurve,
                                    self._forDiscountCurve,
                                    model)[self._deltaMethodString]

            sigma_K_25_D_P = self.volFunction(self._K_25_D_P[i], i)
            model = FinFXModelBlackScholes(sigma_K_25_D_P)
            put._strikeFXRate = self._K_25_D_P[i]
            delta_put = put.delta(self._valueDate,
                                  self._spotFXRate,
                                  self._domDiscountCurve,
                                  
                                  self._forDiscountCurve,
                                  model)[self._deltaMethodString]

            print("Value Date:", self._valueDate)
            print("Expiry Date:", self._expiryDates[i])
            print("T:", self._texp[i])
            print("S0:", self._spotFXRate)
            print("F0T:", self._F0T[i])
            print("ATM Method:", self._atmMethod)
            print("DeltaMethod:", self._deltaMethod)
            print("DeltaATM:", self._deltaATM[i])

            # Print key strikes and their volatilities
            print("Parameters:", self._parameters[i])

            sigma_ATM = self.volFunction(self._K_ATM[i], i)
            print("K_ATM:       %9.7f    Vol: %9.4f"
                  % (self._K_ATM[i], 100.0*sigma_ATM))

            sigma_K_25_D_P = self.volFunction(self._K_25_D_P[i], i)
            print("K_25_D_P:    %9.7f    Vol: %9.4f    Delta: %9.8f"
                  % (self._K_25_D_P[i], 100.0*sigma_K_25_D_P, delta_put))

            sigma_K_25_D_C = self.volFunction(self._K_25_D_C[i], i)
            print("K_25_D_C:    %9.7f    Vol: %9.4f    Delta: %9.8f"
                  % (self._K_25_D_C[i], 100.0*sigma_K_25_D_C, delta_call))

            sigma_RR = sigma_K_25_D_C - sigma_K_25_D_P            
            print("RR: %9.5f" % (100.0*sigma_RR))

            sigma_K_25_D_P_MS = self.volFunction(self._K_25_D_P_MS[i], i)
            print("K_25_D_P_MS: %9.7f    Vol: %9.4f"
                  % (self._K_25_D_P_MS[i], 100.0*sigma_K_25_D_P_MS))

            sigma_K_25_D_C_MS = self.volFunction(self._K_25_D_C_MS[i], i)            
            print("K_25_D_C_MS: %9.7f    Vol: %9.4f"
                  % (self._K_25_D_C_MS[i], 100.0*sigma_K_25_D_C_MS))

            sigma_SS = 0.5*(sigma_K_25_D_C + sigma_K_25_D_P) - sigma_ATM
            print("SS: %9.5f" % (100.0*sigma_SS))

            print("V_25_D_MS: %9.6f" % self._V_25_D_MS[i])

###############################################################################

    def plotVolCurves(self):

        plt.figure()

        for tenorIndex in range(0, self._numVolCurves):

            atmVol = self._atmVols[tenorIndex]*100
            msVol = self._mktStrangle25DeltaVols[tenorIndex]*100
            rrVol = self._riskReversal25DeltaVols[tenorIndex]*100

            lowK = self._K_25_D_P[tenorIndex] * 0.95
            highK = self._K_25_D_C[tenorIndex] * 1.05

            strikes = []
            vols = []

            numIntervals = 30
            K = lowK
            dK = (highK - lowK)/numIntervals
            for i in range(0, numIntervals):
                sigma = self.volFunction(K, tenorIndex) * 100.0
                strikes.append(K)
                vols.append(sigma)
                K = K + dK

            labelStr = self._tenors[tenorIndex]
            labelStr += " ATM: " + str(atmVol)[0:6]
            labelStr += " MS: " + str(msVol)[0:6]
            labelStr += " RR: " + str(rrVol)[0:6]

            plt.plot(strikes, vols, label=labelStr)
            plt.xlabel("Strike")
            plt.ylabel("Volatility")

            title = self._currencyPair

            keyStrikes = []
            keyStrikes.append(self._K_25_D_P[tenorIndex])
            keyStrikes.append(self._K_25_D_P_MS[tenorIndex])
            keyStrikes.append(self._K_ATM[tenorIndex])
            keyStrikes.append(self._K_25_D_C[tenorIndex])
            keyStrikes.append(self._K_25_D_C_MS[tenorIndex])

            keyVols = []

            for K in keyStrikes:
                sigma = self.volFunction(K, tenorIndex) * 100.0
                keyVols.append(sigma)

            plt.plot(keyStrikes, keyVols, 'ro')

        plt.title(title)
        plt.legend(loc='upper right')

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("VALUE DATE", self._valueDate)
        s += labelToString("FX RATE", self._spotFXRate)
        s += labelToString("CCY PAIR", self._currencyPair)
        s += labelToString("NOTIONAL CCY", self._notionalCurrency)
        s += labelToString("NUM TENORS", self._numVolCurves)
        s += labelToString("ATM METHOD", self._atmMethod)
        s += labelToString("DELTA METHOD", self._deltaMethod)

        for i in range(0, self._numVolCurves):

            s += "\n"

            s += labelToString("TENOR", self._tenors[i])
            s += labelToString("EXPIRY DATE", self._expiryDates[i])
            s += labelToString("TIME (YRS)", self._texp[i])
            s += labelToString("FWD FX", self._F0T[i])

            s += labelToString("ATM VOLS", self._atmVols[i]*100.0)
            s += labelToString("MS VOLS", self._mktStrangle25DeltaVols[i]*100.0)
            s += labelToString("RR VOLS", self._riskReversal25DeltaVols[i]*100.0)

            s += labelToString("ATM Strike", self._K_ATM[i])
            s += labelToString("ATM Delta", self._deltaATM[i])
 
            s += labelToString("K_ATM", self._K_ATM[i])
            s += labelToString("MS 25D Call Strike", self._K_25_D_C_MS)
            s += labelToString("MS 25D Put Strike", self._K_25_D_P_MS)
            s += labelToString("NEW 25D CALL STRIKE", self._K_25_D_C)
            s += labelToString("NEW 25D PUT STRIKE", self._K_25_D_P)
            s += labelToString("PARAMS", self._parameters[i])

        return s

###############################################################################

    def _print(self):
        ''' Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. '''
        print(self)

###############################################################################


