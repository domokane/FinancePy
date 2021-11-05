##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from scipy.optimize import minimize

import matplotlib.pyplot as plt
from numba import njit, float64, int64

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.global_vars import gDaysInYear
from ...utils.global_types import OptionTypes
from ...models.option_implied_dbn import option_implied_dbn
from ...utils.helpers import check_argument_types, label_to_string
from ...market.curves.discount_curve import DiscountCurve

from ...models.volatility_fns import VolFunctionTypes
from ...models.volatility_fns import vol_function_clark
from ...models.volatility_fns import vol_function_bloomberg
from ...models.volatility_fns import vol_function_svi
from ...models.volatility_fns import vol_function_ssvi
from ...models.sabr import vol_function_sabr
from ...models.sabr import vol_function_sabr_beta_one
from ...models.sabr import vol_function_sabr_beta_half

from ...utils.math import norminvcdf

from ...models.black_scholes_analytic import bs_delta

from ...utils.distribution import FinDistribution

from ...utils.solver_1d import newton_secant
from ...utils.solver_nm import nelder_mead
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
    params and specified by the vol_type_value. We fit at one time slice only.
    """

    s = args[0]
    t = args[1]
    r = args[2]
    q = args[3]
    strikes = args[4]
    index = args[5]
    volatility_grid = args[6]
    vol_type_value = args[7]

    f = s * np.exp((r-q)*t)

    num_strikes = len(volatility_grid[0])

    tot = 0.0

    for i in range(0, num_strikes):
        fittedVol = vol_function(vol_type_value, params, f, strikes[i], t)
        mkt_vol = volatility_grid[index][i]
        diff = fittedVol - mkt_vol
        tot += diff**2

    return tot

###############################################################################
# Do not cache this function as it leads to complaints
###############################################################################


def _solve_to_horizon(s, t, r, q,
                      strikes,
                      timeIndex,
                      volatility_grid,
                      vol_type_value,
                      x_inits,
                      finSolverType):

    ###########################################################################
    # Determine parameters of vol surface using minimisation
    ###########################################################################

    tol = 1e-6

    args = (s, t, r, q, strikes, timeIndex, volatility_grid, vol_type_value)

    # Nelmer-Mead (both SciPy & Numba) is quicker, but occasionally fails
    # to converge, so for those cases try again with CG
    # Numba version is quicker, but can be slightly away from CG output
    try:
        if finSolverType == FinSolverTypes.NELDER_MEAD_NUMBA:
            xopt = nelder_mead(_obj, np.array(x_inits),
                               bounds=np.array([[], []]).T,
                               args=args, tol_f=tol,
                               tol_x=tol, max_iter=1000)
        elif finSolverType == FinSolverTypes.NELDER_MEAD:
            opt = minimize(_obj, x_inits, args, method="Nelder-Mead", tol=tol)
            xopt = opt.x
        elif finSolverType == FinSolverTypes.CONJUGATE_GRADIENT:
            opt = minimize(_obj, x_inits, args, method="CG", tol=tol)
            xopt = opt.x
    except:
        # If convergence fails try again with CG if necessary
        if finSolverType != FinSolverTypes.CONJUGATE_GRADIENT:
            print('Failed to converge, will try CG')
            opt = minimize(_obj, x_inits, args, method="CG", tol=tol)
            xopt = opt.x

    params = np.array(xopt)
    return params

###############################################################################


@njit(float64(int64, float64[:], float64, float64, float64),
      cache=True, fastmath=True)
def vol_function(vol_function_type_value, params, f, k, t):
    """ Return the volatility for a strike using a given polynomial
    interpolation following Section 3.9 of Iain Clark book. """

    if vol_function_type_value == VolFunctionTypes.CLARK.value:
        vol = vol_function_clark(params, f, k, t)
        return vol
    elif vol_function_type_value == VolFunctionTypes.SABR_BETA_ONE.value:
        vol = vol_function_sabr_beta_one(params, f, k, t)
        return vol
    elif vol_function_type_value == VolFunctionTypes.SABR_BETA_HALF.value:
        vol = vol_function_sabr_beta_half(params, f, k, t)
        return vol
    elif vol_function_type_value == VolFunctionTypes.BBG.value:
        vol = vol_function_bloomberg(params, f, k, t)
        return vol
    elif vol_function_type_value == VolFunctionTypes.SABR.value:
        vol = vol_function_sabr(params, f, k, t)
        return vol
    elif vol_function_type_value == VolFunctionTypes.CLARK5.value:
        vol = vol_function_clark(params, f, k, t)
        return vol
    elif vol_function_type_value == VolFunctionTypes.SVI.value:
        vol = vol_function_svi(params, f, k, t)
        return vol
    elif vol_function_type_value == VolFunctionTypes.SSVI.value:
        vol = vol_function_ssvi(params, f, k, t)
        return vol
    else:
        raise FinError("Unknown Model Type")

###############################################################################


@njit(cache=True, fastmath=True)
def _delta_fit(k, *args):
    """ This is the objective function used in the determination of the 
    option implied strike which is computed in the class below. I map it into
    inverse normcdf space to avoid the flat slope of this function at low vol
    and high K. It speeds up the code as it allows initial values close to
    the solution to be used. """

    vol_type_value = args[0]
    s = args[1]
    t = args[2]
    r = args[3]
    q = args[4]
    option_type_value = args[5]
    inverseDeltaTarget = args[6]
    params = args[7]

    f = s * np.exp((r-q)*t)
    v = vol_function(vol_type_value, params, f, k, t)
    delta_out = bs_delta(s, t, k, r, q, v, option_type_value)
    inverseDeltaOut = norminvcdf(np.abs(delta_out))
    invObjFn = inverseDeltaTarget - inverseDeltaOut

    return invObjFn

###############################################################################
# Unable to cache this function due to dynamic globals warning. Revisit.
###############################################################################


@njit(float64(float64, float64, float64, float64, int64, int64, float64,
              float64, float64[:]), fastmath=True)
def _solver_for_smile_strike(s, t, r, q,
                             option_type_value,
                             volatilityTypeValue,
                             delta_target,
                             initialGuess,
                             parameters):
    """ Solve for the strike that sets the delta of the option equal to the
    target value of delta allowing the volatility to be a function of the
    strike. """

    inverseDeltaTarget = norminvcdf(np.abs(delta_target))

    argtuple = (volatilityTypeValue, s, t, r, q,
                option_type_value,
                inverseDeltaTarget,
                parameters)

    K = newton_secant(_delta_fit, x0=initialGuess, args=argtuple,
                      tol=1e-8, maxiter=50)

    return K

###############################################################################
# Unable to cache function and if I remove njit it complains about pickle
###############################################################################


class EquityVolSurface:
    """ Class to perform a calibration of a chosen parametrised surface to the
    prices of equity options at different strikes and expiry tenors. There is a 
    choice of volatility function from cubic in delta to full SABR and SSVI. 
    Check out VolFunctionTypes. Visualising the volatility curve is useful. 
    Also, there is no guarantee that the implied pdf will be positive."""

    def __init__(self,
                 valuation_date: Date,
                 stock_price: float,
                 discount_curve: DiscountCurve,
                 dividend_curve: DiscountCurve,
                 expiry_dates: (list),
                 strikes: (list, np.ndarray),
                 volatility_grid: (list, np.ndarray),
                 volatility_function_type: VolFunctionTypes = VolFunctionTypes.CLARK,
                 finSolverType: FinSolverTypes = FinSolverTypes.NELDER_MEAD):
        """ Create the EquitySurface object by passing in market vol data
        for a list of strikes and expiry dates. """

        check_argument_types(self.__init__, locals())

        self._valuation_date = valuation_date
        self._stock_price = stock_price

        self._discount_curve = discount_curve
        self._dividend_curve = dividend_curve

        nExpiryDates = len(expiry_dates)
        nStrikes = len(strikes)
        n = len(volatility_grid)
        m = len(volatility_grid[0])

        if n != nExpiryDates:
            raise FinError("1st dimension of vol grid is not nExpiryDates")

        if m != nStrikes:
            raise FinError("2nd dimension of the vol matrix is not nStrikes")

        self._strikes = strikes
        self._num_strikes = len(strikes)

        self._expiry_dates = expiry_dates
        self._numExpiryDates = len(expiry_dates)

        self._volatility_grid = volatility_grid
        self._volatility_function_type = volatility_function_type

        self._build_vol_surface(finSolverType=finSolverType)

###############################################################################

    def volatility_from_strike_date(self, K, expiry_date):
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

        vol_type_value = self._volatility_function_type.value

        index0 = 0  # lower index in bracket
        index1 = 0  # upper index in bracket

        num_curves = self._numExpiryDates

        if num_curves == 1:

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

        else:  # Otherwise we look for bracketing times and interpolate

            for i in range(1, num_curves):

                if texp <= self._texp[i] and texp > self._texp[i-1]:
                    index0 = i-1
                    index1 = i
                    break

        fwd0 = self._F0T[index0]
        fwd1 = self._F0T[index1]

        t0 = self._texp[index0]
        t1 = self._texp[index1]

        vol0 = vol_function(vol_type_value, self._parameters[index0],
                            fwd0, K, t0)

        if index1 != index0:

            vol1 = vol_function(vol_type_value, self._parameters[index1],
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

    # def delta_to_strike(self, callDelta, expiry_date, deltaMethod):
    #     """ Interpolates the strike at a delta and expiry date. Linear
    #     interpolation is used in strike."""

    #     texp = (expiry_date - self._valuation_date) / gDaysInYear

    #     vol_type_value = self._volatility_function_type.value

    #     s = self._spot_fx_rate

    #     if deltaMethod is None:
    #         delta_method_value = self._deltaMethod.value
    #     else:
    #         delta_method_value = deltaMethod.value

    #     index0 = 0 # lower index in bracket
    #     index1 = 0 # upper index in bracket

    #     num_curves = self._num_vol_curves

    #     # If there is only one time horizon then assume flat vol to this time
    #     if num_curves == 1:

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

    #         for i in range(1, num_curves):

    #             if texp <= self._texp[i] and texp > self._texp[i-1]:
    #                 index0 = i-1
    #                 index1 = i
    #                 break

    #     #######################################################################

    #     t0 = self._texp[index0]
    #     t1 = self._texp[index1]

    #     initialGuess = self._K_ATM[index0]

    #     K0 = _solver_for_smile_strike(s, texp, self._rd[index0], self._rf[index0],
    #                               OptionTypes.EUROPEAN_CALL.value,
    #                               vol_type_value, callDelta,
    #                               delta_method_value,
    #                               initialGuess,
    #                               self._parameters[index0],
    #                               self._strikes[index0],
    #                               self._gaps[index0])

    #     if index1 != index0:

    #         K1 = _solver_for_smile_strike(s, texp,
    #                                   self._rd[index1],
    #                                   self._rf[index1],
    #                                   OptionTypes.EUROPEAN_CALL.value,
    #                                   vol_type_value, callDelta,
    #                                   delta_method_value,
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

    def volatility_from_delta_date(self, callDelta, expiry_date,
                                   deltaMethod=None):
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

        vol_type_value = self._volatility_function_type.value

        s = self._stock_price

        index0 = 0  # lower index in bracket
        index1 = 0  # upper index in bracket

        num_curves = self._numExpiryDates

        # If there is only one time horizon then assume flat vol to this time
        if num_curves == 1:

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

        else:  # Otherwise we look for bracketing times and interpolate

            for i in range(1, num_curves):

                if texp <= self._texp[i] and texp > self._texp[i-1]:
                    index0 = i-1
                    index1 = i
                    break

        fwd0 = self._F0T[index0]
        fwd1 = self._F0T[index1]

        t0 = self._texp[index0]
        t1 = self._texp[index1]

        initialGuess = self._stock_price

        K0 = _solver_for_smile_strike(s,
                                      texp,
                                      self._r[index0],
                                      self._q[index0],
                                      OptionTypes.EUROPEAN_CALL.value,
                                      vol_type_value, callDelta,
                                      initialGuess,
                                      self._parameters[index0])

        vol0 = vol_function(vol_type_value, self._parameters[index0],
                            fwd0, K0, t0)

        if index1 != index0:

            K1 = _solver_for_smile_strike(s, texp,
                                          self._r[index1],
                                          self._q[index1],
                                          OptionTypes.EUROPEAN_CALL.value,
                                          vol_type_value, callDelta,
                                          initialGuess,
                                          self._parameters[index1])

            vol1 = vol_function(vol_type_value, self._parameters[index1],
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
                raise FinError(
                    "Failed interpolation due to negative variance.")

            volt = np.sqrt(vart/texp)

        else:

            volt = vol0
            kt = K0

        return volt, kt

###############################################################################

    def _build_vol_surface(self, finSolverType=FinSolverTypes.NELDER_MEAD):
        """ Main function to construct the vol surface. """

        s = self._stock_price

        numExpiryDates = self._numExpiryDates

        if self._volatility_function_type == VolFunctionTypes.CLARK:
            num_parameters = 3
            self._parameters = np.zeros([numExpiryDates, num_parameters])
        elif self._volatility_function_type == VolFunctionTypes.SABR_BETA_ONE:
            num_parameters = 3
            self._parameters = np.zeros([numExpiryDates, num_parameters])
        elif self._volatility_function_type == VolFunctionTypes.SABR_BETA_HALF:
            num_parameters = 3
            self._parameters = np.zeros([numExpiryDates, num_parameters])
        elif self._volatility_function_type == VolFunctionTypes.BBG:
            num_parameters = 3
            self._parameters = np.zeros([numExpiryDates, num_parameters])
        elif self._volatility_function_type == VolFunctionTypes.SABR:
            num_parameters = 4
            self._parameters = np.zeros([numExpiryDates, num_parameters])
        elif self._volatility_function_type == VolFunctionTypes.CLARK5:
            num_parameters = 5
            self._parameters = np.zeros([numExpiryDates, num_parameters])
        elif self._volatility_function_type == VolFunctionTypes.SVI:
            num_parameters = 5
            self._parameters = np.zeros([numExpiryDates, num_parameters])
        elif self._volatility_function_type == VolFunctionTypes.SSVI:
            num_parameters = 5
            self._parameters = np.zeros([numExpiryDates, num_parameters])
            self._parameters[:, 0] = 0.2  # sigma
            self._parameters[:, 1] = 0.8  # gamma
            self._parameters[:, 2] = -0.7  # rho
            self._parameters[:, 3] = 0.3
            self._parameters[:, 4] = 0.048
        else:
            print(self._volatilityFunctionType)
            raise FinError("Unknown Model Type")

        self._texp = np.zeros(numExpiryDates)

        self._F0T = np.zeros(numExpiryDates)
        self._r = np.zeros(numExpiryDates)
        self._q = np.zeros(numExpiryDates)

        #######################################################################
        # TODO: ADD SPOT DAYS
        #######################################################################

        spot_date = self._valuation_date

        for i in range(0, numExpiryDates):

            expiry_date = self._expiry_dates[i]
            texp = (expiry_date - spot_date) / gDaysInYear

            disDF = self._discount_curve._df(texp)
            divDF = self._dividend_curve._df(texp)
            f = s * divDF/disDF

            self._texp[i] = texp
            self._r[i] = -np.log(disDF) / texp
            self._q[i] = -np.log(divDF) / texp
            self._F0T[i] = f

        #######################################################################
        # THE ACTUAL COMPUTATION LOOP STARTS HERE
        #######################################################################

        vol_type_value = self._volatility_function_type.value

        x_inits = []
        x_init = np.zeros(num_parameters)
        x_inits.append(x_init)

        for i in range(0, numExpiryDates):

            t = self._texp[i]
            r = self._r[i]
            q = self._q[i]

            res = _solve_to_horizon(s, t, r, q,
                                    self._strikes,
                                    i,
                                    self._volatility_grid,
                                    vol_type_value,
                                    x_inits[i],
                                    finSolverType)

            self._parameters[i, :] = res

            x_init = res
            x_inits.append(x_init)

###############################################################################

    def check_calibration(self, verbose: bool, tol: float = 1e-6):
        """ Compare calibrated vol surface with market and output a report
        which sets out the quality of fit to the ATM and 10 and 25 delta market
        strangles and risk reversals. """

        if verbose:

            print("==========================================================")
            print("VALUE DATE:", self._valuation_date)
            print("STOCK PRICE:", self._stock_price)
            print("==========================================================")

        for i in range(0, self._numExpiryDates):

            expiry_date = self._expiry_dates[i]
            print("==========================================================")

            for j in range(0, self._num_strikes):

                strike = self._strikes[j]

                fittedVol = self.volatility_from_strike_date(
                    strike, expiry_date)

                mkt_vol = self._volatility_grid[i][j]

                diff = fittedVol - mkt_vol

                print("%s %12.3f %7.4f %7.4f %7.5f" %
                      (expiry_date, strike,
                       fittedVol*100.0, mkt_vol*100, diff*100))

        print("==========================================================")

###############################################################################

    def implied_dbns(self, lowS, highS, numIntervals):
        """ Calculate the pdf for each tenor horizon. Returns a list of
        FinDistribution objects, one for each tenor horizon. """

        dbns = []

        for iTenor in range(0, self._numExpiryDates):

            f = self._F0T[iTenor]
            t = self._texp[iTenor]

            dS = (highS - lowS) / numIntervals

            disDF = self._discount_curve._df(t)
            divDF = self._dividend_curve._df(t)

            r = -np.log(disDF) / t
            q = -np.log(divDF) / t

            Ks = []
            vols = []

            for iK in range(0, numIntervals):

                k = lowS + iK*dS

                vol = vol_function(self._volatility_function_type.value,
                                   self._parameters[iTenor],
                                   f, k, t)

                Ks.append(k)
                vols.append(vol)

            Ks = np.array(Ks)
            vols = np.array(vols)

            density = option_implied_dbn(self._stock_price, t, r, q, Ks, vols)

            dbn = FinDistribution(Ks, density)
            dbns.append(dbn)

        return dbns

###############################################################################

    def plot_vol_curves(self):
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
                fittedVol = self.volatility_from_strike_date(
                    K, expiry_date) * 100.
                fittedVols.append(fittedVol)
                K = K + dK

            labelStr = "FITTED AT " + str(self._expiry_dates[tenorIndex])
            plt.plot(ks, fittedVols, label=labelStr)

            labelStr = "MARKET AT " + str(self._expiry_dates[tenorIndex])
            mkt_vols = self._volatility_grid[tenorIndex] * 100.0
            plt.plot(self._strikes, mkt_vols, 'o', label=labelStr)

            plt.xlabel("Strike")
            plt.ylabel("Volatility")

            title = str(self._volatility_function_type)
            plt.title(title)
            plt.legend()

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("VALUE DATE", self._valuation_date)
        s += label_to_string("STOCK PRICE", self._stock_price)
        s += label_to_string("VOL FUNCTION", self._volatility_function_type)

        for i in range(0, self._numExpiryDates):
            s += label_to_string("EXPIRY DATE", self._expiry_dates[i])

        for i in range(0, self._num_strikes):
            s += label_to_string("STRIKE", self._strikes[i])

        s += label_to_string("EQUITY VOL GRID", self._volatility_grid)

        return s

###############################################################################

    def _print(self):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """
        print(self)

###############################################################################
