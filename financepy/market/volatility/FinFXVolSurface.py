##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from scipy.optimize import newton
from scipy.optimize import minimize

import matplotlib.pyplot as plt
from numba import njit, float64, int64

from ...finutils.FinError import FinError
from ...finutils.FinDate import FinDate
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinGlobalTypes import FinOptionTypes
from ...products.fx.FinFXVanillaOption import FinFXVanillaOption
from ...models.FinModelOptionImpliedDbn import optionImpliedDbn
from ...products.fx.FinFXMktConventions import FinFXATMMethod
from ...products.fx.FinFXMktConventions import FinFXDeltaMethod
from ...finutils.FinHelperFunctions import checkArgumentTypes, labelToString
from ...market.curves.FinDiscountCurve import FinDiscountCurve

from ...products.fx.FinFXModelTypes import FinFXModelBlackScholes
from ...market.volatility.FinOptionVolatilityFns import volFunctionClarke
from ...market.volatility.FinOptionVolatilityFns import volFunctionSABR
from ...market.volatility.FinOptionVolatilityFns import FinVolFunctionTypes

from ...finutils.FinMath import norminvcdf

from ...models.FinModelBlackScholesAnalytical import bsValue
from ...products.fx.FinFXVanillaOption import fastDelta
    
###############################################################################
# TODO: Speed up search for strike by providing derivative function to go with
#       delta fit.
###############################################################################

def g(K, *args):
    ''' This is the objective function used in the determination of the FX
    option implied strike which is computed in the class below. '''

    s = args[0]
    t = args[1]
    rd = args[2]
    rf = args[3]
    volatility = args[4]
    deltaMethodValue = args[5]
    optionTypeValue = args[6]
    deltaTarget = args[7]

    deltaOut = fastDelta(s, t, K, rd, rf,
                         volatility,
                         deltaMethodValue,
                         optionTypeValue)
    
    objFn = deltaTarget - deltaOut    
    return objFn

###############################################################################

# def obj(cvec, *args):
#     ''' Return a function that is minimised when the ATM, MS and RR vols have
#     been best fitted using the parametric volatility curve represented by cvec
#     '''
#     self = args[0]
#     i = self._tenorIndex
#     self._parameters[i, :] = cvec

#     # We first need to solve for the strikes at the 25 delta points using the
#     # new volatility curve

#     # Match the at-the-money option volatility
#     atmCurveVol = self.volFunction(self._K_ATM[i], i)
#     targetATMVol = self._atmVols[i]
#     term1 = (targetATMVol - atmCurveVol)**2

#     # Match the market strangle value but this has to be at the MS strikes
#     # call._strikeFXRate = self._K_25D_C_MS[i]
#     sigma_K_25D_C_MS = self.volFunction(self._K_25D_C_MS[i], i)
#     # model = FinFXModelBlackScholes(sigma_K_25D_C_MS)

#     # V_C_25_MS = call.value(self._valueDate,
#     #                        self._spotFXRate,
#     #                        self._domDiscountCurve,
#     #                        self._forDiscountCurve,
#     #                        model)['v']

#     V_25D_C_MS = bsValue(self._spotFXRate,
#                           self._texp[i],
#                           self._K_25D_C_MS[i],
#                           self._rd[i],
#                           self._rf[i],
#                           sigma_K_25D_C_MS,
#                           FinOptionTypes.EUROPEAN_CALL.value)

#     # put._strikeFXRate = self._K_25D_P_MS[i]
#     sigma_K_25D_P_MS = self.volFunction(self._K_25D_P_MS[i], i)
#     # model = FinFXModelBlackScholes(sigma_K_25D_P_MS)

#     # V_25D_P_MS = put.value(self._valueDate,
#     #                       self._spotFXRate,
#     #                       self._domDiscountCurve,
#     #                       self._forDiscountCurve,
#     #                       model)['v']

#     V_25D_P_MS = bsValue(self._spotFXRate,
#                           self._texp[i],
#                           self._K_25D_P_MS[i],
#                           self._rd[i],
#                           self._rf[i],
#                           sigma_K_25D_P_MS,
#                           FinOptionTypes.EUROPEAN_PUT.value)

#     V_25D_MS = V_25D_C_MS + V_25D_P_MS
#     target25StrangleValue = self._V_25D_MS[i]
#     term2 = (V_25D_MS - target25StrangleValue)**2

#     # Match the risk reversal volatility
    
#     K0 = self._K_25D_C_MS[i]
#     K_25D_C = self.solveForSmileStrike(FinOptionTypes.EUROPEAN_CALL.value, +0.2500, i, K0)
#     sigma_K_25D_C = self.volFunction(K_25D_C, i)

#     K0 = self._K_25D_P_MS[i]
#     K_25D_P = self.solveForSmileStrike(FinOptionTypes.EUROPEAN_PUT.value, -0.2500, i, K0)
#     sigma_K_25D_P = self.volFunction(K_25D_P, i)

#     targetRRVol = self._riskReversal25DeltaVols[i]
#     sigma_25D_RR = (sigma_K_25D_C - sigma_K_25D_P)
#     term3 = (sigma_25D_RR - targetRRVol)**2

#     self._K_25D_C[i] = K_25D_C
#     self._K_25D_P[i] = K_25D_P

#     tot = term1 + term2 + term3

#     print("OBJSLOW", cvec, tot)

# #    print("TOT:", tot)
#     return tot

###############################################################################

def objFAST(params, *args):
    ''' Return a function that is minimised when the ATM, MS and RR vols have
    been best fitted using the parametric volatility curve represented by cvec
    '''

    s = args[0]
    t = args[1]
    rd = args[2]
    rf = args[3]
    K_ATM = args[4]
    atmVol = args[5]
    K_25D_C_MS = args[6]
    K_25D_P_MS = args[7]
    V_25D_MS_target = args[8]
    deltaMethodValue = args[9]
    targetRRVol = args[10]
    volTypeValue = args[11]

    f = s * np.exp((rd-rf)*t)
    # We first need to solve for the strikes at the 25 delta points using the
    # new volatility curve

    # Match the at-the-money option volatility
    atmCurveVol = volFunctionFAST(volTypeValue, params, f, K_ATM, t)
    term1 = (atmVol - atmCurveVol)**2

    ###########################################################################
    # Match the market strangle value but this has to be at the MS strikes
    ###########################################################################

    sigma_K_25D_C_MS = volFunctionFAST(volTypeValue, params, f, K_25D_C_MS, t)

    V_25D_C_MS = bsValue(s, t, K_25D_C_MS, rd, rf, sigma_K_25D_C_MS,
                         FinOptionTypes.EUROPEAN_CALL.value)

    sigma_K_25D_P_MS = volFunctionFAST(volTypeValue, params, f, K_25D_P_MS, t)

    V_25D_P_MS = bsValue(s, t, K_25D_P_MS, rd, rf, sigma_K_25D_P_MS,
                         FinOptionTypes.EUROPEAN_PUT.value)

    V_25D_MS = V_25D_C_MS + V_25D_P_MS
    term2 = (V_25D_MS - V_25D_MS_target)**2

    ###########################################################################
    # Match the risk reversal volatility    
    ###########################################################################

    K_25D_C = solveForSmileStrikeFAST(s, t, rd, rf, 
                                      FinOptionTypes.EUROPEAN_CALL.value, 
                                      volTypeValue, +0.2500, 
                                      deltaMethodValue, K_25D_C_MS, 
                                      params)
 
    sigma_K_25D_C = volFunctionFAST(volTypeValue, params, f, K_25D_C, t)

    K_25D_P = solveForSmileStrikeFAST(s, t, rd, rf, 
                                      FinOptionTypes.EUROPEAN_PUT.value,
                                      volTypeValue, -0.2500, 
                                      deltaMethodValue, K_25D_P_MS,
                                      params)

    sigma_K_25D_P = volFunctionFAST(volTypeValue, params, f, K_25D_P, t)

    sigma_25D_RR = (sigma_K_25D_C - sigma_K_25D_P)
    term3 = (sigma_25D_RR - targetRRVol)**2

    ###########################################################################

    tot = term1 + term2 + term3

    return tot

###############################################################################

def solveToHorizonFAST(s, t, rd, rf, K_ATM,
                       atmVol, ms25DVol, rr25DVol, 
                       deltaMethodValue, volTypeValue, 
                       xopt):

    c0 = xopt

    # Determine the price of a market strangle from market strangle
    # Need to price a call and put that agree with market strangle

    vol_25D_MS = atmVol + ms25DVol

    K_25D_C_MS = solveForStrike(s, t, rd, rf,
                                FinOptionTypes.EUROPEAN_CALL.value,
                                +0.2500,
                                deltaMethodValue,
                                vol_25D_MS)

    K_25D_P_MS = solveForStrike(s, t, rd, rf,
                                FinOptionTypes.EUROPEAN_PUT.value,
                                -0.2500,
                                deltaMethodValue,
                                vol_25D_MS)

    # USE MARKET STRANGLE VOL TO DETERMINE PRICE OF A MARKET STRANGLE
    V_25D_C_MS = bsValue(s, t, K_25D_C_MS, rd, rf, vol_25D_MS,
                         FinOptionTypes.EUROPEAN_CALL.value)

    V_25D_P_MS = bsValue(s, t, K_25D_P_MS, rd, rf, vol_25D_MS,
                         FinOptionTypes.EUROPEAN_PUT.value)

    # Market price of strangle in the domestic currency
    V_25D_MS = V_25D_C_MS + V_25D_P_MS

    # Determine parameters of vol surface using minimisation
    tol = 1e-8

    args = (s, t, rd, rf, K_ATM, atmVol, K_25D_C_MS, K_25D_P_MS, V_25D_MS,
            deltaMethodValue, rr25DVol, volTypeValue)

    opt = minimize(objFAST, c0, args, method="CG", tol=tol)
    xopt = opt.x

    v = objFAST(xopt, *args)

    # The function is quadratic. So the quadratic value should be small
    if abs(v) > tol:
        print("Using an xtol:", tol, " maybe make this smaller")
        print("Reported value of v:", v)
        raise FinError("Failed to fit volatility smile curve.")

    params = np.array(xopt)

    K_25D_C = solveForSmileStrikeFAST(s, t, rd, rf, 
                                      FinOptionTypes.EUROPEAN_CALL.value, 
                                      volTypeValue, +0.2500, 
                                      deltaMethodValue, K_25D_C_MS, 
                                      params)
 
    K_25D_P = solveForSmileStrikeFAST(s, t, rd, rf, 
                                      FinOptionTypes.EUROPEAN_PUT.value,
                                      volTypeValue, -0.2500, 
                                      deltaMethodValue, K_25D_P_MS,
                                      params)

    return (params, K_25D_C_MS, K_25D_P_MS, K_25D_C, K_25D_P)

###############################################################################

@njit(float64(int64, float64[:], float64, float64, float64), 
      cache=True, fastmath=True)
def volFunctionFAST(volFunctionTypeValue, params, f, k, t):
    ''' Return the volatility for a strike using a given polynomial
    interpolation following Section 3.9 of Iain Clark book. '''

    if volFunctionTypeValue == FinVolFunctionTypes.CLARKE.value:
        vol = volFunctionClarke(params, f, k, t)
        return vol
    elif volFunctionTypeValue == FinVolFunctionTypes.SABR.value:
        vol = volFunctionSABR(params, f, k, t)
        return vol
    else:
        raise FinError("Unknown Model Type")

###############################################################################

def deltaFit(K, *args):
    ''' This is the objective function used in the determination of the FX
    Option implied strike which is computed in the class below. I map it into
    inverse normcdf space to avoid the flat slope of this function at low vol 
    and high K. It speeds up the code as it allows initial values close to 
    the solution to be used. '''

    volTypeValue = args[0]
    s = args[1]
    t = args[2]
    rd = args[3]
    rf = args[4]
    optionTypeValue = args[5]
    deltaTypeValue = args[6]
    inverseDeltaTarget = args[7]
    params = args[8]

    f = s*np.exp((rd-rf)*t)
    v = volFunctionFAST(volTypeValue, params, f, K, t)
    deltaOut = fastDelta(s, t, K, rd, rf, v, deltaTypeValue, optionTypeValue)
    inverseDeltaOut = norminvcdf(np.abs(deltaOut))   
    invObjFn = inverseDeltaTarget - inverseDeltaOut

    return invObjFn

###############################################################################

def solveForSmileStrikeFAST(s, t, rd, rf,
                            optionTypeValue,
                            volatilityTypeValue,
                            deltaTarget,
                            deltaMethodValue,
                            initialValue,
                            parameters):
    ''' Solve for the strike that sets the delta of the option equal to the
    target value of delta allowing the volatility to be a function of the
    strike. '''

    inverseDeltaTarget = norminvcdf(np.abs(deltaTarget))
    
    argtuple = (volatilityTypeValue, s, t, rd, rf, optionTypeValue, 
                deltaMethodValue, inverseDeltaTarget, parameters)

    K = newton(deltaFit, x0=initialValue, args=argtuple,
               tol=1e-8, maxiter=50, fprime2=None)

    return K
    
###############################################################################

def solveForStrike(spotFXRate, tdel, rd, rf,
                   optionTypeValue,
                   deltaTarget,
                   deltaMethodValue, 
                   volatility):
    ''' This function determines the implied strike of an FX option
    given a delta and the other option details. It uses a one-dimensional
    Newton root search algorith to determine the strike that matches an
    input volatility. '''
        
    # =========================================================================
    # IMPORTANT NOTE:
    # =========================================================================
    # For some delta quotation conventions I can solve for K explicitly.
    # Note that as I am using the function norminvdelta to calculate the 
    # inverse value of delta, this may not, on a round trip using N(x), give
    # back the value x as it is calculated to a different number of decimal 
    # places. It should however agree to 6-7 decimal places. Which is OK.
    # =========================================================================

    if deltaMethodValue == FinFXDeltaMethod.SPOT_DELTA.value:

        domDF = np.exp(-rd*tdel) 
        forDF = np.exp(-rf*tdel) 

        if optionTypeValue == FinOptionTypes.EUROPEAN_CALL.value:
            phi = +1.0
        else:
            phi = -1.0

        F0T = spotFXRate * forDF / domDF
        vsqrtt = volatility * np.sqrt(tdel)
        arg = np.abs(deltaTarget*phi/forDF)  # CHECK THIS !!!
        norminvdelta = norminvcdf(arg)
        K = F0T * np.exp(-vsqrtt *(phi*norminvdelta - vsqrtt/2.0))
        return K

    elif deltaMethodValue == FinFXDeltaMethod.FORWARD_DELTA.value:

        domDF = np.exp(-rd*tdel) 
        forDF = np.exp(-rf*tdel) 

        if optionTypeValue == FinOptionTypes.EUROPEAN_CALL.value:
            phi = +1.0
        else:
            phi = -1.0

        F0T = spotFXRate * forDF / domDF
        vsqrtt = volatility * np.sqrt(tdel)
        arg = np.abs(deltaTarget*phi)   # CHECK THIS!!!!!!!!      
        norminvdelta = norminvcdf(arg) 
        K = F0T * np.exp(-vsqrtt *(phi*norminvdelta - vsqrtt/2.0))
        return K
 
    elif deltaMethodValue == FinFXDeltaMethod.SPOT_DELTA_PREM_ADJ.value:

        argtuple = (spotFXRate, tdel, rd, rf, volatility,
                    deltaMethodValue, optionTypeValue, deltaTarget)
    
        K = newton(g, x0=spotFXRate, args=argtuple,
                   tol=1e-7, rtol=1e-7, maxiter=50, fprime2=None)

        return K

    elif deltaMethodValue == FinFXDeltaMethod.FORWARD_DELTA_PREM_ADJ.value:

        argtuple = (spotFXRate, tdel, rd, rf, volatility,
                    deltaMethodValue, optionTypeValue, deltaTarget)
    
        K = newton(g, x0=spotFXRate, args=argtuple,
                   tol=1e-7, rtol=1e-7, maxiter=50, fprime2=None)

        return K

    else:

        raise FinError("Unknown FinFXDeltaMethod")

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
                 deltaMethod=FinFXDeltaMethod.SPOT_DELTA, 
                 volatilityFunctionType=FinVolFunctionTypes.CLARKE):

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

        self._volatilityFunctionType = volatilityFunctionType
        self._tenorIndex = 0

        self._expiryDates = []
        for i in range(0, self._numVolCurves):
            expiryDate = valueDate.addTenor(tenors[i])
            self._expiryDates.append(expiryDate)

        self.buildVolSurface()

###############################################################################

    # def volFunction(self, K, tenorIndex):
    #     ''' Return the volatility for a strike using a given polynomial
    #     interpolation following Section 3.9 of Iain Clark book. '''

    #     params = self._parameters[tenorIndex]
    #     F = self._F0T[tenorIndex]
    #     texp = self._texp[tenorIndex]

    #     if self._volatilityFunctionType == FinVolFunctionTypes.CLARKE:
    #         vol = volFunctionClarke(params, F, K, texp)
    #     elif self._volatilityFunctionType == FinVolFunctionTypes.SABR:
    #         vol = volFunctionSABR(params, F, K, texp)
    #     else:
    #         raise FinError("Unknown Model Type")

    #     return vol

###############################################################################

    def buildVolSurface(self):

        s = self._spotFXRate
        numVolCurves = self._numVolCurves

        if self._volatilityFunctionType == FinVolFunctionTypes.CLARKE:
            numParameters = 3
        elif self._volatilityFunctionType == FinVolFunctionTypes.SABR:
            numParameters = 3
        else:
            raise FinError("Unknown Model Type")

        self._parameters = np.zeros([numVolCurves, numParameters])

        self._F0T = np.zeros(numVolCurves)
        self._rd = np.zeros(numVolCurves)
        self._rf = np.zeros(numVolCurves)
        self._K_ATM = np.zeros(numVolCurves)
        self._deltaATM = np.zeros(numVolCurves)

        self._K_25D_C = np.zeros(numVolCurves)
        self._K_25D_P = np.zeros(numVolCurves)
        self._K_25D_C_MS = np.zeros(numVolCurves)
        self._K_25D_P_MS = np.zeros(numVolCurves)
        self._V_25D_MS = np.zeros(numVolCurves)       
        self._texp = np.zeros(numVolCurves)
                        
        #######################################################################
        # TODO: ADD SPOT DAYS
        #######################################################################
        spotDate = self._valueDate

        for i in range(0, numVolCurves):

            expiryDate = self._expiryDates[i]
            t = (expiryDate - spotDate) / gDaysInYear

            domDF = self._domDiscountCurve._df(t)
            forDF = self._forDiscountCurve._df(t)
            f = s * forDF/domDF

            self._texp[i] = t
            self._rd[i] = -np.log(domDF) / t
            self._rf[i] = -np.log(forDF) / t
            self._F0T[i] = f

            atmVol = self._atmVols[i]

            # This follows exposition in Clarke Page 52
            if self._atmMethod == FinFXATMMethod.SPOT:
                self._K_ATM[i] = s
            elif self._atmMethod == FinFXATMMethod.FWD:
                self._K_ATM[i] = f
            elif self._atmMethod == FinFXATMMethod.FWD_DELTA_NEUTRAL:
                self._K_ATM[i] = f * np.exp(atmVol*atmVol*t/2.0)
            elif self._atmMethod == FinFXATMMethod.FWD_DELTA_NEUTRAL_PREM_ADJ:
                self._K_ATM[i] = f * np.exp(-atmVol*atmVol*t/2.0)
            else:
                raise FinError("Unknown Delta Type")

        #######################################################################
        # THE ACTUAL COMPUTATION LOOP STARTS HERE
        #######################################################################

        xinits = []
        for i in range(0, numVolCurves):

            atmVol = self._atmVols[i]

            if self._volatilityFunctionType == FinVolFunctionTypes.CLARKE:
                xinit = [np.log(atmVol), 0.050, 0.800]
            elif self._volatilityFunctionType == FinVolFunctionTypes.SABR:
                # SABR parameters are alpha, nu, rho
                xinit = [0.18, 0.82, -0.12]
            else:
                raise FinError("Unknown Model Type")

            xinits.append(xinit)


        deltaMethodValue = self._deltaMethod.value
        volTypeValue = self._volatilityFunctionType.value

        for i in range(0, numVolCurves):
            
            t = self._texp[i]
            rd = self._rd[i]
            rf = self._rf[i]
            K_ATM = self._K_ATM[i]
            atmVol = self._atmVols[i]
            ms25DVol = self._mktStrangle25DeltaVols[i]
            rr25DVol = self._riskReversal25DeltaVols[i]

#            print(t, rd, rf, K_ATM, atmVol, ms25DVol, rr25DVol)

            res = solveToHorizonFAST(s, t, rd, rf, K_ATM, 
                                   atmVol, ms25DVol, rr25DVol, 
                                   deltaMethodValue, volTypeValue,
                                   xinits[i])

            (self._parameters[i,:], 
             self._K_25D_C_MS[i], self._K_25D_P_MS[i], 
             self._K_25D_C[i], self._K_25D_P[i]) = res

###############################################################################

    # def solveToHorizon(self, i, xopt):

    #     print("Starting solve To Horizon", i)
    #     self._tenorIndex = i

    #     c0 = xopt

    #     # Determine the price of a market strangle from market strangle
    #     # Need to price a call and put that agree with market strangle

    #     vol_25D_MS = self._atmVols[i] + self._mktStrangle25DeltaVols[i]

    #     K_25D_C_MS = solveForStrike(self._spotFXRate,
    #                                 self._texp[i],                                         
    #                                 self._rd[i],
    #                                 self._rf[i],
    #                                 FinOptionTypes.EUROPEAN_CALL.value,
    #                                 +0.2500,
    #                                 self._deltaMethod.value,
    #                                 vol_25D_MS)

    #     K_25D_P_MS = solveForStrike(self._spotFXRate,
    #                                 self._texp[i],
    #                                 self._rd[i],
    #                                 self._rf[i],
    #                                 FinOptionTypes.EUROPEAN_PUT.value,
    #                                 -0.2500,
    #                                 self._deltaMethod.value,
    #                                 vol_25D_MS)

    #     # Store the set of strikes in the class
    #     self._K_25D_C_MS[i] = K_25D_C_MS
    #     self._K_25D_P_MS[i] = K_25D_P_MS

    #     # USE MARKET STRANGLE VOL TO DETERMINE PRICE OF A MARKET STRANGLE

    #     V_25D_C_MS = bsValue(self._spotFXRate,
    #                          self._texp[i],
    #                          K_25D_C_MS,
    #                          self._rd[i],
    #                          self._rf[i],
    #                          vol_25D_MS,
    #                          FinOptionTypes.EUROPEAN_CALL.value)

    #     V_25D_P_MS = bsValue(self._spotFXRate,
    #                          self._texp[i],
    #                          K_25D_P_MS,
    #                          self._rd[i],
    #                          self._rf[i],
    #                          vol_25D_MS,
    #                          FinOptionTypes.EUROPEAN_PUT.value)

    #     # Market price of strangle in the domestic currency
    #     V_25D_MS = V_25D_C_MS + V_25D_P_MS
    #     self._V_25D_MS[i] = V_25D_MS

    #     if 1 == 0:
    #         print("Value Date:", self._valueDate)
    #         print("Expiry Date:", self._expiryDates[i])
    #         print("T:", self._texp[i])
    #         print("S0:", self._spotFXRate)
    #         print("F0T:", self._F0T[i])
    #         print("ATM Method:", self._atmMethod)
    #         print("DeltaATM:", self._deltaATM[i])
    #         print("DeltaMethod:", self._deltaMethod)

    #     # Determine parameters of vol surface using minimisation
    #     args = (self) #, call, put)
    #     tol = 1e-8

    #     opt = minimize(obj, c0, args, method="CG", tol=tol)
    #     xopt = opt.x

    #     v = obj(xopt, args)

    #     # The function is quadratic. So the quadratic value should be small
    #     if abs(v) > tol:
    #         print("Using an xtol:", tol, " maybe make this smaller")
    #         raise FinError("Failed to fit volatility smile curve.")

    #     self._parameters[i, :] = np.array(xopt)

###############################################################################

    def solveForSmileStrike(self,
                            optionTypeValue,
                            deltaTarget,
                            tenorIndex, 
                            initialValue):
        ''' Solve for the strike that sets the delta of the option equal to the
        target value of delta allowing the volatility to be a function of the
        strike. '''

        s0 = self._spotFXRate
        tdel = self._texp[tenorIndex]
        rd = self._rd[tenorIndex]
        rf = self._rf[tenorIndex]

        inverseDeltaTarget = norminvcdf(np.abs(deltaTarget))
        argtuple = (self, s0, tdel, rd, rf, optionTypeValue, 
                    inverseDeltaTarget, tenorIndex)
        
        volTypeValue = self._volatilityFunctionType.value

        argtuple = (volTypeValue, s0, tdel, rd, rf, optionTypeValue, 
                    self._deltaMethod.value, 
                    inverseDeltaTarget, self._parameters[tenorIndex])

        K = newton(deltaFit, x0=initialValue, args=argtuple,
                       tol=1e-5, maxiter=50, fprime2=None)

        return K

###############################################################################

    def checkCalibration(self, verbose: bool, tol: float = 1e-6):

        if verbose:

            print("==========================================================")
            print("VALUE DATE:", self._valueDate)
            print("SPOT FX RATE:", self._spotFXRate)
            print("ATM METHOD:", self._atmMethod)
            print("DELTA METHOD:", self._deltaMethod)
            print("==========================================================")
        
        K_dummy = 999

        for i in range(0, self._numVolCurves):

            expiryDate = self._expiryDates[i]

            if verbose:
                print("TENOR:", self._tenors[i])
                print("EXPIRY DATE:", expiryDate)
                print("IN ATM VOL: %9.6f %%"% (100.0*self._atmVols[i]))
                print("IN MKT STRANGLE 25D VOL: %9.6f %%"% (100.0*self._mktStrangle25DeltaVols[i]))
                print("IN RSK REVERSAL 25D VOL: %9.6f %%"% (100.0*self._riskReversal25DeltaVols[i]))

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


            ###################################################################
            # AT THE MONEY
            ###################################################################

            if verbose:
                print("==========================================================")
                print("T_(YEARS): ", self._texp[i])
                print("CNT_CPD_RD:%9.6f %%"% (self._rd[i]*100))
                print("CNT_CPD_RF:%9.6f %%"% (self._rf[i]*100))
                print("FWD_RATE:  %9.6f"% (self._F0T[i]))
 
            sigma_ATM_out = volFunctionFAST(self._volatilityFunctionType.value,
                                        self._parameters[i],
                                        self._F0T[i],
                                        self._K_ATM[i],
                                        self._texp[i])

            if verbose:
                print("==========================================================")
                print("VOL FUNCTION", self._volatilityFunctionType)
                print("VOL_PARAMETERS:", self._parameters[i])
                print("==========================================================")
                print("OUT_K_ATM:  %9.6f" % (self._K_ATM[i]))
                print("OUT_ATM_VOL: %9.6f %%" 
                      % (100.0*sigma_ATM_out))

            diff = sigma_ATM_out - self._atmVols[i]

            if np.abs(diff) > tol:
                print("FAILED FIT TO ATM VOL IN: %9.6f  OUT: %9.6f  DIFF: %9.6f"% 
                      (self._atmVols[i]*100.0, sigma_ATM_out*100.0, 
                       diff * 100.0))

            call._strikeFXRate = self._K_ATM[i]
            put._strikeFXRate = self._K_ATM[i]

            model = FinFXModelBlackScholes(sigma_ATM_out)

            delta_call = call.delta(self._valueDate,
                                    self._spotFXRate,
                                    self._domDiscountCurve,
                                    self._forDiscountCurve,
                                    model)[self._deltaMethodString]

            delta_put = put.delta(self._valueDate,
                                  self._spotFXRate,
                                  self._domDiscountCurve,
                                  self._forDiscountCurve,
                                  model)[self._deltaMethodString]

            if verbose:
                print("CALL_DELTA: % 9.6f  PUT_DELTA: % 9.6f  NET_DELTA: % 9.6f" 
                      % (delta_call, delta_put, delta_call + delta_put))

            ###################################################################
            # NOW WE ASSIGN THE SAME VOLATILITY TO THE MS STRIKES
            # THESE STRIKES ARE DETERMINED BY SETTING DELTA TO 0.25/-0.25
            ###################################################################


            msVol = self._atmVols[i] + self._mktStrangle25DeltaVols[i]

            if verbose:

                print("==========================================================")
                print("MKT STRANGLE VOL IN: %9.6f %%"
                      % (100.0*self._mktStrangle25DeltaVols[i]))

            call._strikeFXRate = self._K_25D_C_MS[i]
            put._strikeFXRate = self._K_25D_P_MS[i]

            model = FinFXModelBlackScholes(msVol)

            delta_call = call.delta(self._valueDate,
                                    self._spotFXRate,
                                    self._domDiscountCurve,
                                    self._forDiscountCurve,
                                    model)[self._deltaMethodString]

            delta_put = put.delta(self._valueDate,
                                  self._spotFXRate,
                                  self._domDiscountCurve,
                                  self._forDiscountCurve,
                                  model)[self._deltaMethodString]

            if verbose:
                print("K_25D_C_MS: %9.6f  ATM + MSVOL: %9.6f %%   DELTA: %9.6f"
                      % (self._K_25D_C_MS[i], 100.0*msVol, delta_call))
    
                print("K_25D_P_MS: %9.6f  ATM + MSVOL: %9.6f %%   DELTA: %9.6f"
                      % (self._K_25D_P_MS[i], 100.0*msVol, delta_put))

            call_value = call.value(self._valueDate,
                                    self._spotFXRate,
                                    self._domDiscountCurve,
                                    self._forDiscountCurve,
                                    model)['v']

            put_value = put.value(self._valueDate,
                                  self._spotFXRate,
                                  self._domDiscountCurve,
                                  self._forDiscountCurve,
                                  model)['v']

            mktStrangleValue = call_value + put_value

            if verbose:
                print("CALL_VALUE: %9.6f  PUT_VALUE: %9.6f  MS_VALUE: % 9.6f" 
                      % (call_value, put_value, mktStrangleValue))

            ###################################################################
            # NOW WE ASSIGN A DIFFERENT VOLATILITY TO THE MS STRIKES
            # THE DELTAS WILL NO LONGER EQUAL 0.25, -0.25
            ###################################################################

            # CALL
            sigma_K_25D_C_MS = volFunctionFAST(self._volatilityFunctionType.value,
                                            self._parameters[i],
                                            self._F0T[i],
                                            self._K_25D_C_MS[i],
                                            self._texp[i])
 
            model = FinFXModelBlackScholes(sigma_K_25D_C_MS)
            call_value = call.value(self._valueDate,
                                    self._spotFXRate,
                                    self._domDiscountCurve,
                                    self._forDiscountCurve,
                                    model)['v']

            # THIS IS NOT GOING TO BE 0.25 AS WE HAVE USED A DIFFERENT SKEW VOL
            delta_call = call.delta(self._valueDate,
                                    self._spotFXRate,
                                    self._domDiscountCurve,
                                    self._forDiscountCurve,
                                    model)[self._deltaMethodString]

            # PUT
            sigma_K_25D_P_MS = volFunctionFAST(self._volatilityFunctionType.value,
                                            self._parameters[i],
                                            self._F0T[i],
                                            self._K_25D_P_MS[i],
                                            self._texp[i])

        
            model = FinFXModelBlackScholes(sigma_K_25D_P_MS)
            put_value = put.value(self._valueDate,
                                  self._spotFXRate,
                                  self._domDiscountCurve,
                                  self._forDiscountCurve,
                                  model)['v']

            # THIS IS NOT GOING TO BE -0.25 AS WE HAVE USED A DIFFERENT SKEW VOL
            delta_put = put.delta(self._valueDate,
                                  self._spotFXRate,
                                  self._domDiscountCurve,
                                  self._forDiscountCurve,
                                  model)[self._deltaMethodString]

            mktStrangleValueSkew = call_value + put_value

            if verbose:
                print("K_25D_C_MS: %9.6f  SURFACE_VOL: %9.6f %%   DELTA: %9.6f"
                      % (self._K_25D_C_MS[i], 100.0*sigma_K_25D_C_MS, delta_call))
    
                print("K_25D_P_MS: %9.6f  SURFACE_VOL: %9.6f %%   DELTA: %9.6f"
                      % (self._K_25D_P_MS[i], 100.0*sigma_K_25D_P_MS, delta_put))
    
                print("CALL_VALUE: %9.6f  PUT_VALUE: %9.6f  MS_SKEW_VALUE: % 9.6f" 
                      % (call_value, put_value, mktStrangleValueSkew))

            diff = mktStrangleValue - mktStrangleValueSkew
            if np.abs(diff) > tol:
                print("FAILED FIT TO 25D MS VAL: %9.6f  OUT: %9.6f  DIFF: % 9.6f"%
                      (mktStrangleValue, mktStrangleValueSkew, diff))

            ###################################################################
            # NOW WE SHIFT STRIKES SO THAT DELTAS NOW EQUAL 0.25, -0.25
            ###################################################################
        
            call._strikeFXRate = self._K_25D_C[i]
            put._strikeFXRate = self._K_25D_P[i]

            sigma_K_25D_C = volFunctionFAST(self._volatilityFunctionType.value,
                                            self._parameters[i],
                                            self._F0T[i],
                                            self._K_25D_C[i],
                                            self._texp[i])
 
            model = FinFXModelBlackScholes(sigma_K_25D_C)

            # THIS DELTA SHOULD BE +0.25
            delta_call = call.delta(self._valueDate,
                                    self._spotFXRate,
                                    self._domDiscountCurve,
                                    self._forDiscountCurve,
                                    model)[self._deltaMethodString]

            sigma_K_25D_P = volFunctionFAST(self._volatilityFunctionType.value,
                                            self._parameters[i],
                                            self._F0T[i],
                                            self._K_25D_P[i],
                                            self._texp[i])

            model = FinFXModelBlackScholes(sigma_K_25D_P)

            # THIS DELTA SHOULD BE -0.25
            delta_put = put.delta(self._valueDate,
                                  self._spotFXRate,
                                  self._domDiscountCurve,
                                  self._forDiscountCurve,
                                  model)[self._deltaMethodString]
            
            if verbose:
                print("K_25D_C: %9.7f  VOL: %9.6f  DELTA: % 9.6f"
                      % (self._K_25D_C[i], 100.0*sigma_K_25D_C, delta_call))
                
                print("K_25D_P: %9.7f  VOL: %9.6f  DELTA: % 9.6f"
                      % (self._K_25D_P[i], 100.0*sigma_K_25D_P, delta_put))


            sigma_RR = sigma_K_25D_C - sigma_K_25D_P            

            if verbose:
                print("==========================================================")
                print("RR = VOL_K_25_C - VOL_K_25_P => RR_IN: %9.6f %% RR_OUT: %9.6f %%" 
                  % (100.0 * self._riskReversal25DeltaVols[i], 100.0*sigma_RR))
                print("==========================================================")

            diff = sigma_RR - self._riskReversal25DeltaVols[i]

            if np.abs(diff) > tol:
                print("FAILED FIT TO 25D RRV IN: % 9.6f  OUT: % 9.6f  DIFF: % 9.6f"% 
                      (self._riskReversal25DeltaVols[i]*100.0,
                       sigma_RR*100.0, 
                       diff*100.0))

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

    def plotVolCurves(self):

        plt.figure()

        volTypeVal = self._volatilityFunctionType.value

        for tenorIndex in range(0, self._numVolCurves):

            atmVol = self._atmVols[tenorIndex]*100
            msVol = self._mktStrangle25DeltaVols[tenorIndex]*100
            rrVol = self._riskReversal25DeltaVols[tenorIndex]*100

            lowK = self._K_25D_P[tenorIndex] * 0.95
            highK = self._K_25D_C[tenorIndex] * 1.05

            strikes = []
            vols = []
            numIntervals = 30
            K = lowK
            dK = (highK - lowK)/numIntervals
            params = self._parameters[tenorIndex]
            t = self._texp[tenorIndex]
            f = self._F0T[tenorIndex]

            for i in range(0, numIntervals):
                sigma = volFunctionFAST(volTypeVal, params, f, K, t) * 100.0
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
            keyStrikes.append(self._K_25D_P[tenorIndex])
            keyStrikes.append(self._K_25D_P_MS[tenorIndex])
            keyStrikes.append(self._K_ATM[tenorIndex])
            keyStrikes.append(self._K_25D_C[tenorIndex])
            keyStrikes.append(self._K_25D_C_MS[tenorIndex])

            keyVols = []

            for K in keyStrikes:
                sigma = volFunctionFAST(volTypeVal, params, f, K, t) * 100.0
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
        s += labelToString("VOL FUNCTION", self._volatilityFunctionType)

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
            s += labelToString("MS 25D Call Strike", self._K_25D_C_MS[i])
            s += labelToString("MS 25D Put Strike", self._K_25D_P_MS[i])
            s += labelToString("NEW 25D CALL STRIKE", self._K_25D_C[i])
            s += labelToString("NEW 25D PUT STRIKE", self._K_25D_P[i])
            s += labelToString("PARAMS", self._parameters[i])

        return s

###############################################################################

    def _print(self):
        ''' Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. '''
        print(self)

###############################################################################

if __name__ == '__main__':
    pass
