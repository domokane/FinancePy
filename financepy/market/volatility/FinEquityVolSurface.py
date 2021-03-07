##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from scipy.optimize import minimize

import matplotlib.pyplot as plt
from numba import njit, float64, int64

from ...utils.FinError import FinError
from ...utils.date import Date
from ...utils.global_vars import gDaysInYear
from ...utils.global_types import FinOptionTypes
from ...models.FinModelOptionImpliedDbn import optionImpliedDbn
from ...utils.helpers import check_argument_types, labelToString
from ...market.discount.curve import DiscountCurve

from ...models.volatility_fns import FinVolFunctionTypes
from ...models.volatility_fns import volFunctionClark
from ...models.volatility_fns import volFunctionBloomberg
from ...models.volatility_fns import volFunctionSVI
from ...models.volatility_fns import volFunctionSSVI
from ...models.sabr import volFunctionSABR
from ...models.sabr import volFunctionSABR_BETA_ONE
from ...models.sabr import volFunctionSABR_BETA_HALF

from ...utils.math import norminvcdf

from ...models.black_scholes_analytic import bsDelta

from ...utils.distribution import FinDistribution

from ...utils.FinSolvers1D import newton_secant
from ...utils.FinSolversNM import nelder_mead
from ...utils.global_types import FinSolverTypes

###############################################################################
# ISSUES
###############################################################################
 
###############################################################################
# TODO: Speed up search for strike by providing derivative function to go with
#       delta fit.
###############################################################################


###############################################################################
# Do not cache this function - WRONG. IT WORKS BUT WHY WHEN IN FX IT FAILS ??
###############################################################################

@njit(fastmath=True, cache=True)
def _obj(params, *args):
    """ Return a value that is minimised when the ATM, MS and RR vols have
    been best fitted using the parametric volatility curve represented by
    params and specified by the volTypeValue. We fit at one time slice only.
    """

    s = args[0]
    t = args[1]
    r = args[2]
    q = args[3]
    strikes = args[4]
    index = args[5]
    volatilityGrid = args[6]
    volTypeValue = args[7]

    f = s * np.exp((r-q)*t)

    numStrikes = len(volatilityGrid[0])
 
    tot = 0.0

    for i in range(0, numStrikes):
            
        fittedVol = volFunction(volTypeValue, params, f, strikes[i], t)

        mktVol = volatilityGrid[index][i]
            
        diff = fittedVol - mktVol
            
        tot += diff**2

    return tot

###############################################################################
# Do not cache this function as it leads to complaints
###############################################################################

def _solveToHorizon(s, t, r, q,
                    strikes,
                    timeIndex,
                    volatilityGrid,
                    volTypeValue,
                    xinits,
                    finSolverType):

    ###########################################################################
    # Determine parameters of vol surface using minimisation
    ###########################################################################

    tol = 1e-6

    args = (s, t, r, q, strikes, timeIndex, volatilityGrid, volTypeValue)

    # Nelmer-Mead (both SciPy & Numba) is quicker, but occasionally fails 
    # to converge, so for those cases try again with CG
    # Numba version is quicker, but can be slightly away from CG output
    try:
        if finSolverType == FinSolverTypes.NELDER_MEAD_NUMBA:
            xopt = nelder_mead(_obj, np.array(xinits), 
                               bounds=np.array([[], []]).T, 
                               args=args, tol_f=tol,
                               tol_x=tol, max_iter=1000)
        elif finSolverType == FinSolverTypes.NELDER_MEAD:
            opt = minimize(_obj, xinits, args, method="Nelder-Mead", tol=tol)
            xopt = opt.x
        elif finSolverType == FinSolverTypes.CONJUGATE_GRADIENT:
            opt = minimize(_obj, xinits, args, method="CG", tol=tol)
            xopt = opt.x
    except:
        # If convergence fails try again with CG if necessary
        if finSolverType != FinSolverTypes.CONJUGATE_GRADIENT:
            print('Failed to converge, will try CG')
            opt = minimize(_obj, xinits, args, method="CG", tol=tol)
            xopt = opt.x

    params = np.array(xopt)    
    return params

###############################################################################


@njit(float64(int64, float64[:], float64, float64, float64), 
      cache=True, fastmath=True)
def volFunction(volFunctionTypeValue, params, f, k, t):
    """ Return the volatility for a strike using a given polynomial
    interpolation following Section 3.9 of Iain Clark book. """

    if volFunctionTypeValue == FinVolFunctionTypes.CLARK.value:
        vol = volFunctionClark(params, f, k, t)
        return vol
    elif volFunctionTypeValue == FinVolFunctionTypes.SABR_BETA_ONE.value:
        vol = volFunctionSABR_BETA_ONE(params, f, k, t)
        return vol
    elif volFunctionTypeValue == FinVolFunctionTypes.SABR_BETA_HALF.value:
        vol = volFunctionSABR_BETA_HALF(params, f, k, t)
        return vol
    elif volFunctionTypeValue == FinVolFunctionTypes.BBG.value:
        vol = volFunctionBloomberg(params, f, k, t)
        return vol
    elif volFunctionTypeValue == FinVolFunctionTypes.SABR.value:
        vol = volFunctionSABR(params, f, k, t)
        return vol
    elif volFunctionTypeValue == FinVolFunctionTypes.CLARK5.value:
        vol = volFunctionClark(params, f, k, t)
        return vol
    elif volFunctionTypeValue == FinVolFunctionTypes.SVI.value:
        vol = volFunctionSVI(params, f, k, t)
        return vol
    elif volFunctionTypeValue == FinVolFunctionTypes.SSVI.value:
        vol = volFunctionSSVI(params, f, k, t)
        return vol
    else:
        raise FinError("Unknown Model Type")

###############################################################################


@njit(cache=True, fastmath=True)
def _deltaFit(k, *args):
    """ This is the objective function used in the determination of the 
    option implied strike which is computed in the class below. I map it into
    inverse normcdf space to avoid the flat slope of this function at low vol
    and high K. It speeds up the code as it allows initial values close to
    the solution to be used. """

    volTypeValue = args[0]
    s = args[1]
    t = args[2]
    r = args[3]
    q = args[4]
    optionTypeValue = args[5]
    inverseDeltaTarget = args[6]
    params = args[7]

    f = s * np.exp((r-q)*t)
    v = volFunction(volTypeValue, params, f, k, t)
    deltaOut = bsDelta(s, t, k, r, q, v, optionTypeValue)
    inverseDeltaOut = norminvcdf(np.abs(deltaOut))
    invObjFn = inverseDeltaTarget - inverseDeltaOut

    return invObjFn

###############################################################################
# Unable to cache this function due to dynamic globals warning. Revisit.
###############################################################################


@njit(float64(float64, float64, float64, float64, int64, int64, float64,
              float64, float64[:]), fastmath=True)
def _solveForSmileStrike(s, t, r, q,
                         optionTypeValue,
                         volatilityTypeValue,
                         deltaTarget,
                         initialGuess,
                         parameters):
     """ Solve for the strike that sets the delta of the option equal to the
     target value of delta allowing the volatility to be a function of the
     strike. """

     inverseDeltaTarget = norminvcdf(np.abs(deltaTarget))

     argtuple = (volatilityTypeValue, s, t, r, q, 
                 optionTypeValue, 
                 inverseDeltaTarget, 
                 parameters)

     K = newton_secant(_deltaFit, x0=initialGuess, args=argtuple,
                       tol=1e-8, maxiter=50)

     return K

###############################################################################
# Unable to cache function and if I remove njit it complains about pickle
###############################################################################


class FinEquityVolSurface():
    """ Class to perform a calibration of a chosen parametrised surface to the
    prices of equity options at different strikes and expiry tenors. There is a 
    choice of volatility function from cubic in delta to full SABR and SSVI. 
    Check out FinVolFunctionTypes. Visualising the volatility curve is useful. 
    Also, there is no guarantee that the implied pdf will be positive."""

    def __init__(self,
                 valuation_date: Date,
                 stock_price: float,
                 discount_curve: DiscountCurve,
                 dividendCurve: DiscountCurve,
                 expiry_dates: (list),
                 strikes: (list, np.ndarray),
                 volatilityGrid: (list, np.ndarray),
                 volatilityFunctionType:FinVolFunctionTypes=FinVolFunctionTypes.CLARK,
                 finSolverType:FinSolverTypes=FinSolverTypes.NELDER_MEAD):
        """ Create the FinEquitySurface object by passing in market vol data
        for a list of strikes and expiry dates. """

        check_argument_types(self.__init__, locals())

        self._valuation_date = valuation_date
        self._stock_price = stock_price

        self._discount_curve = discount_curve
        self._dividendCurve = dividendCurve

        nExpiryDates = len(expiry_dates)
        nStrikes = len(strikes)
        n = len(volatilityGrid)
        m = len(volatilityGrid[0])
        
        if n != nExpiryDates:
            raise FinError("1st dimension of vol grid is not nExpiryDates")

        if m != nStrikes:
            raise FinError("2nd dimension of the vol matrix is not nStrikes")

        self._strikes = strikes
        self._numStrikes = len(strikes)

        self._expiry_dates = expiry_dates
        self._numExpiryDates = len(expiry_dates)

        self._volatilityGrid = volatilityGrid    
        self._volatilityFunctionType = volatilityFunctionType

        self._buildVolSurface(finSolverType=finSolverType)

###############################################################################

    def volatilityFromStrikeDate(self, K, expiry_date):
        """ Interpolates the Black-Scholes volatility from the volatility
        surface given call option strike and expiry date. Linear interpolation
        is done in variance space. The smile strikes at bracketed dates are 
        determined by determining the strike that reproduces the provided delta
        value. This uses the calibration delta convention, but it can be 
        overriden by a provided delta convention. The resulting volatilities 
        are then determined for each bracketing expiry time and linear 
        interpolation is done in variance space and then converted back to a 
        lognormal volatility."""

        texp = (expiry_date - self._valuation_date) / gDaysInYear

        volTypeValue = self._volatilityFunctionType.value

        index0 = 0 # lower index in bracket
        index1 = 0 # upper index in bracket

        numCurves = self._numExpiryDates

        if numCurves == 1:

            index0 = 0
            index1 = 0
            
        # If the time is below first time then assume a flat vol
        elif texp <= self._texp[0]:

            index0 = 0
            index1 = 0

        # If the time is beyond the last time then extrapolate with a flat vol
        elif texp >= self._texp[-1]:

            index0 = len(self._texp) - 1
            index1 = len(self._texp) - 1

        else: # Otherwise we look for bracketing times and interpolate

            for i in range(1, numCurves):

                if texp <= self._texp[i] and texp > self._texp[i-1]:
                    index0 = i-1
                    index1 = i
                    break

        fwd0 = self._F0T[index0]
        fwd1 = self._F0T[index1]
                
        t0 = self._texp[index0]
        t1 = self._texp[index1]

        vol0 = volFunction(volTypeValue, self._parameters[index0],
                               fwd0, K, t0)

        if index1 != index0:

            vol1 = volFunction(volTypeValue, self._parameters[index1],
                               fwd1, K, t1)

        else:
            
            vol1 = vol0

        # In the expiry time dimension, both volatilities are interpolated 
        # at the same strikes but different deltas.
        vart0 = vol0*vol0*t0
        vart1 = vol1*vol1*t1

        if np.abs(t1-t0) > 1e-6:
            vart = ((texp-t0) * vart1 + (t1-texp) * vart0) / (t1 - t0)

            if vart < 0.0:
                raise FinError("Negative variance.")

            volt = np.sqrt(vart/texp)

        else:
            volt = vol1

        return volt

###############################################################################

    # def deltaToStrike(self, callDelta, expiry_date, deltaMethod):
    #     """ Interpolates the strike at a delta and expiry date. Linear 
    #     interpolation is used in strike."""

    #     texp = (expiry_date - self._valuation_date) / gDaysInYear

    #     volTypeValue = self._volatilityFunctionType.value

    #     s = self._spotFXRate

    #     if deltaMethod is None:
    #         deltaMethodValue = self._deltaMethod.value
    #     else:
    #         deltaMethodValue = deltaMethod.value

    #     index0 = 0 # lower index in bracket
    #     index1 = 0 # upper index in bracket

    #     numCurves = self._numVolCurves

    #     # If there is only one time horizon then assume flat vol to this time
    #     if numCurves == 1:

    #         index0 = 0
    #         index1 = 0

    #     # If the time is below first time then assume a flat vol
    #     elif texp <= self._texp[0]:

    #         index0 = 0
    #         index1 = 0

    #     # If the time is beyond the last time then extrapolate with a flat vol
    #     elif texp > self._texp[-1]:
 
    #         index0 = len(self._texp) - 1
    #         index1 = len(self._texp) - 1

    #     else: # Otherwise we look for bracketing times and interpolate

    #         for i in range(1, numCurves):

    #             if texp <= self._texp[i] and texp > self._texp[i-1]:
    #                 index0 = i-1
    #                 index1 = i
    #                 break

    #     #######################################################################
                
    #     t0 = self._texp[index0]
    #     t1 = self._texp[index1]

    #     initialGuess = self._K_ATM[index0]

    #     K0 = _solveForSmileStrike(s, texp, self._rd[index0], self._rf[index0],
    #                               FinOptionTypes.EUROPEAN_CALL.value,
    #                               volTypeValue, callDelta,
    #                               deltaMethodValue,
    #                               initialGuess,
    #                               self._parameters[index0], 
    #                               self._strikes[index0], 
    #                               self._gaps[index0])

    #     if index1 != index0:

    #         K1 = _solveForSmileStrike(s, texp, 
    #                                   self._rd[index1], 
    #                                   self._rf[index1],
    #                                   FinOptionTypes.EUROPEAN_CALL.value,
    #                                   volTypeValue, callDelta,
    #                                   deltaMethodValue,
    #                                   initialGuess,
    #                                   self._parameters[index1], 
    #                                   self._strikes[index1], 
    #                                   self._gaps[index1])
    #     else:

    #         K1 = K0
            
    #     # In the expiry time dimension, both volatilities are interpolated 
    #     # at the same strikes but different deltas.
 
    #     if np.abs(t1-t0) > 1e-6:

    #         K = ((texp-t0) * K1 + (t1-texp) * K1) / (K1 - K0)

    #     else:

    #         K = K1

    #     return K

###############################################################################
        
    def volatilityFromDeltaDate(self, callDelta, expiry_date, 
                                deltaMethod = None):
        """ Interpolates the Black-Scholes volatility from the volatility
        surface given a call option delta and expiry date. Linear interpolation
        is done in variance space. The smile strikes at bracketed dates are 
        determined by determining the strike that reproduces the provided delta
        value. This uses the calibration delta convention, but it can be 
        overriden by a provided delta convention. The resulting volatilities 
        are then determined for each bracketing expiry time and linear 
        interpolation is done in variance space and then converted back to a 
        lognormal volatility."""

        texp = (expiry_date - self._valuation_date) / gDaysInYear

        volTypeValue = self._volatilityFunctionType.value

        s = self._stock_price

        index0 = 0 # lower index in bracket
        index1 = 0 # upper index in bracket

        numCurves = self._numExpiryDates

        # If there is only one time horizon then assume flat vol to this time
        if numCurves == 1:

            index0 = 0
            index1 = 0

        # If the time is below first time then assume a flat vol
        elif texp <= self._texp[0]:

            index0 = 0
            index1 = 0

        # If the time is beyond the last time then extrapolate with a flat vol
        elif texp > self._texp[-1]:
 
            index0 = len(self._texp) - 1
            index1 = len(self._texp) - 1

        else: # Otherwise we look for bracketing times and interpolate

            for i in range(1, numCurves):

                if texp <= self._texp[i] and texp > self._texp[i-1]:
                    index0 = i-1
                    index1 = i
                    break
        
        fwd0 = self._F0T[index0]
        fwd1 = self._F0T[index1]
                
        t0 = self._texp[index0]
        t1 = self._texp[index1]

        initialGuess = self._stock_price

        K0 = _solveForSmileStrike(s,
                                  texp,
                                  self._r[index0],
                                  self._q[index0],
                                  FinOptionTypes.EUROPEAN_CALL.value,
                                  volTypeValue, callDelta,
                                  initialGuess,
                                  self._parameters[index0])

        vol0 = volFunction(volTypeValue, self._parameters[index0],
                           fwd0, K0, t0)

        if index1 != index0:

            K1 = _solveForSmileStrike(s, texp, 
                                      self._r[index1], 
                                      self._q[index1],
                                      FinOptionTypes.EUROPEAN_CALL.value,
                                      volTypeValue, callDelta,
                                      initialGuess,
                                      self._parameters[index1])

            vol1 = volFunction(volTypeValue, self._parameters[index1], 
                               fwd1, K1, t1)
        else:
            vol1 = vol0
            
        # In the expiry time dimension, both volatilities are interpolated 
        # at the same strikes but different deltas.
        vart0 = vol0*vol0*t0
        vart1 = vol1*vol1*t1

        if np.abs(t1-t0) > 1e-6:

            vart = ((texp-t0) * vart1 + (t1-texp) * vart0) / (t1 - t0)
            kt = ((texp-t0) * K1 + (t1-texp) * K0) / (t1 - t0)

            if vart < 0.0:
                raise FinError("Failed interpolation due to negative variance.")

            volt = np.sqrt(vart/texp)

        else:

            volt = vol0
            kt = K0

        return volt, kt

###############################################################################

    def _buildVolSurface(self, finSolverType=FinSolverTypes.NELDER_MEAD):
        """ Main function to construct the vol surface. """

        s = self._stock_price

        if self._volatilityFunctionType == FinVolFunctionTypes.CLARK:
            numParameters = 3
        elif self._volatilityFunctionType == FinVolFunctionTypes.SABR_BETA_ONE:
            numParameters = 3
        elif self._volatilityFunctionType == FinVolFunctionTypes.SABR_BETA_HALF:
            numParameters = 3
        elif self._volatilityFunctionType == FinVolFunctionTypes.BBG:
            numParameters = 3
        elif self._volatilityFunctionType == FinVolFunctionTypes.SABR:
            numParameters = 4
        elif self._volatilityFunctionType == FinVolFunctionTypes.CLARK5:
            numParameters = 5
        elif self._volatilityFunctionType == FinVolFunctionTypes.SVI:
            numParameters = 5
        elif self._volatilityFunctionType == FinVolFunctionTypes.SSVI:
            numParameters = 5
        else:
            print(self._volatilityFunctionType)
            raise FinError("Unknown Model Type")

        numExpiryDates = self._numExpiryDates

        self._parameters = np.zeros([numExpiryDates, numParameters])        
        self._texp = np.zeros(numExpiryDates)

        self._F0T = np.zeros(numExpiryDates)
        self._r = np.zeros(numExpiryDates)
        self._q = np.zeros(numExpiryDates)

        #######################################################################
        # TODO: ADD SPOT DAYS
        #######################################################################

        spotDate = self._valuation_date

        for i in range(0, numExpiryDates):

            expiry_date = self._expiry_dates[i]
            texp = (expiry_date - spotDate) / gDaysInYear

            disDF = self._discount_curve._df(texp)
            divDF = self._dividendCurve._df(texp)
            f = s * divDF/disDF

            self._texp[i] = texp
            self._r[i] = -np.log(disDF) / texp
            self._q[i] = -np.log(divDF) / texp
            self._F0T[i] = f

        #######################################################################
        # THE ACTUAL COMPUTATION LOOP STARTS HERE
        #######################################################################

        volTypeValue = self._volatilityFunctionType.value

        xinits = []
        xinit = np.zeros(numParameters)
        xinits.append(xinit)

        for i in range(0, numExpiryDates):

            t = self._texp[i]
            r = self._r[i]
            q = self._q[i]

            res = _solveToHorizon(s, t, r, q,
                                  self._strikes,
                                  i,
                                  self._volatilityGrid, 
                                  volTypeValue,
                                  xinits[i],
                                  finSolverType)

            self._parameters[i,:] = res

            xinit = res
            xinits.append(xinit)

###############################################################################

    def checkCalibration(self, verbose: bool, tol: float = 1e-6):
        """ Compare calibrated vol surface with market and output a report
        which sets out the quality of fit to the ATM and 10 and 25 delta market
        strangles and risk reversals. """

        if verbose:

            print("==========================================================")
            print("VALUE DATE:", self._valuation_date)
            print("STOCK PRICE:", self._stock_price)
            print("==========================================================")

        K_dummy = 999

        for i in range(0, self._numExpiryDates):

            expiry_date = self._expiry_dates[i]
            print("==========================================================")

            for j in range(0, self._numStrikes):
                
                strike = self._strikes[j]

                fittedVol = self.volatilityFromStrikeDate(strike, expiry_date)

                mktVol = self._volatilityGrid[i][j]
                
                diff = fittedVol - mktVol
                
                print("%s %12.3f %7.4f %7.4f %7.5f"% 
                      (expiry_date, strike, 
                       fittedVol*100.0, mktVol*100, diff*100))

        print("==========================================================")
        
###############################################################################

    def impliedDbns(self, lowS, highS, numIntervals):
        """ Calculate the pdf for each tenor horizon. Returns a list of
        FinDistribution objects, one for each tenor horizon. """

        dbns = []

        for iTenor in range(0, self._numExpiryDates):

            f = self._F0T[iTenor]
            t = self._texp[iTenor]

            dS = (highS - lowS)/ numIntervals

            disDF = self._discount_curve._df(t)
            divDF = self._dividendCurve._df(t)

            r = -np.log(disDF) / t
            q = -np.log(divDF) / t

            Ks = []
            vols = []

            for iK in range(0, numIntervals):

                k = lowS + iK*dS

                vol = volFunction(self._volatilityFunctionType.value, 
                                  self._parameters[iTenor], 
                                  f, k, t)

                Ks.append(k)
                vols.append(vol)

            Ks = np.array(Ks)
            vols = np.array(vols)

            density = optionImpliedDbn(self._stock_price, t, r, q, Ks, vols)

            dbn = FinDistribution(Ks, density)
            dbns.append(dbn)

        return dbns

###############################################################################

    def plotVolCurves(self):
        """ Generates a plot of each of the vol discount implied by the market
        and fitted. """
        
        lowK = self._strikes[0] * 0.9
        highK = self._strikes[-1] * 1.1

        for tenorIndex in range(0, self._numExpiryDates):

            expiry_date = self._expiry_dates[tenorIndex]
            plt.figure()

            ks = []
            fittedVols = []

            numIntervals = 30
            K = lowK
            dK = (highK - lowK)/numIntervals
            
            for i in range(0, numIntervals):

                ks.append(K)
                fittedVol = self.volatilityFromStrikeDate(K, expiry_date) * 100.
                fittedVols.append(fittedVol)
                K = K + dK

            labelStr = "FITTED AT " + str(self._expiry_dates[tenorIndex])
            plt.plot(ks, fittedVols, label=labelStr)

            labelStr = "MARKET AT " + str(self._expiry_dates[tenorIndex])
            mktVols = self._volatilityGrid[tenorIndex] * 100.0
            plt.plot(self._strikes, mktVols, 'o', label=labelStr)

            plt.xlabel("Strike")
            plt.ylabel("Volatility")

            title = str(self._volatilityFunctionType)
            plt.title(title)
            plt.legend()

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("VALUE DATE", self._valuation_date)
        s += labelToString("STOCK PRICE", self._stock_price)
        s += labelToString("VOL FUNCTION", self._volatilityFunctionType)

        for i in range(0, self._numExpiryDates):
            s += labelToString("EXPIRY DATE", self._expiry_dates[i])

        for i in range(0, self._numStrikes):
            s += labelToString("STRIKE", self._strikes[i])

        s += labelToString("EQUITY VOL GRID", self._volatilityGrid)
 
        return s

###############################################################################

    def _print(self):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """
        print(self)

###############################################################################
