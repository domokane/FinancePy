###############################################################################
# Copyright (C) 2020 Saeed Amen, Dominic O'Kane
###############################################################################

from numba import njit, boolean, int64, float64, vectorize
import numpy as np

###############################################################################
## from https://quanteconpy.readthedocs.io/en/latest/_modules/quantecon/optimize/root_finding.html #####################

_ECONVERGED = 0
_ECONVERR = -1

_iter = 100
_xtol = 2e-12
_rtol = 4*np.finfo(float).eps

from collections import namedtuple

results = namedtuple('results', 'root function_calls iterations converged')

###############################################################################

@njit(cache=True, fastmath=True)
def _results(r):
    r"""Select from a tuple of(root, funccalls, iterations, flag)"""
    x, funcalls, iterations, flag = r
    return x# results(x, funcalls, iterations, flag == 0)

###############################################################################

@njit(fastmath=True, cache=True)
def newton_secant(func, x0, args=(), tol=1.48e-8, maxiter=50,
                  disp=True):
    """
    Find a zero from the secant method using the jitted version of
    Scipy's secant method.

    Note that `func` must be jitted via Numba.

    Parameters
    ----------
    func : callable and jitted
        The function whose zero is wanted. It must be a function of a
        single variable of the form f(x,a,b,c...), where a,b,c... are extra
        arguments that can be passed in the `args` parameter.
    x0 : float
        An initial estimate of the zero that should be somewhere near the
        actual zero.
    args : tuple, optional(default=())
        Extra arguments to be used in the function call.
    tol : float, optional(default=1.48e-8)
        The allowable error of the zero value.
    maxiter : int, optional(default=50)
        Maximum number of iterations.
    disp : bool, optional(default=True)
        If True, raise a RuntimeError if the algorithm didn't converge.

    Returns
    -------
    results : namedtuple
        A namedtuple containing the following items:
        ::

            root - Estimated location where function is zero.
            function_calls - Number of times the function was called.
            iterations - Number of iterations needed to find the root.
            converged - True if the routine converged
    """

    if tol <= 0:
        raise ValueError("tol is too small <= 0")
    if maxiter < 1:
        raise ValueError("maxiter must be greater than 0")

    # Convert to float (don't use float(x0); this works also for complex x0)
    p0 = 1.0 * x0
    funcalls = 0
    status = _ECONVERR

    # Secant method
    if x0 >= 0:
        p1 = x0 * (1 + 1e-4) + 1e-4
    else:
        p1 = x0 * (1 + 1e-4) - 1e-4
        q0 = func(p0, *args)
    funcalls += 1
    q1 = func(p1, *args)
    funcalls += 1
    for itr in range(maxiter):
        if q1 == q0:
            p = (p1 + p0) / 2.0
            status = _ECONVERGED
            break
        else:
            p = p1 - q1 * (p1 - p0) / (q1 - q0)
        if np.abs(p - p1) < tol:
            status = _ECONVERGED
            break
        p0 = p1
        q0 = q1
        p1 = p
        q1 = func(p1, *args)
        funcalls += 1

    if disp and status == _ECONVERR:
        msg = "Failed to converge"
        raise RuntimeError(msg)

    return p#_results((p, funcalls, itr + 1, status))

###############################################################################

@njit(fastmath=True, cache=True)
def brent_max(func, a, b, args=(), xtol=1e-5, maxiter=500):
    """
    Uses a jitted version of the maximization routine from SciPy's fminbound.
    The algorithm is identical except that it's been switched to maximization
    rather than minimization, and the tests for convergence have been stripped
    out to allow for jit compilation.

    Note that the input function `func` must be jitted or the call will fail.

    Parameters
    ----------
    func : jitted function
    a : scalar
        Lower bound for search
    b : scalar
        Upper bound for search
    args : tuple, optional
        Extra arguments passed to the objective function.
    maxiter : int, optional
        Maximum number of iterations to perform.
    xtol : float, optional
        Absolute error in solution `xopt` acceptable for convergence.

    Returns
    -------
    xf : float
        The maximizer
    fval : float
        The maximum value attained
    info : tuple
        A tuple of the form (status_flag, num_iter).  Here status_flag
        indicates whether or not the maximum number of function calls was
        attained.  A value of 0 implies that the maximum was not hit.
        The value `num_iter` is the number of function calls.

    Examples
    --------
    >>> @njit
    ... def f(x):
    ...     return -(x + 2.0)**2 + 1.0
    ...
    >>> xf, fval, info = brent_max(f, -2, 2)

    """
    if not np.isfinite(a):
        raise ValueError("a must be finite.")

    if not np.isfinite(b):
        raise ValueError("b must be finite.")

    if not a < b:
        raise ValueError("a must be less than b.")

    maxfun = maxiter
    status_flag = 0

    sqrt_eps = np.sqrt(2.2e-16)
    golden_mean = 0.5 * (3.0 - np.sqrt(5.0))

    fulc = a + golden_mean * (b - a)
    nfc, xf = fulc, fulc
    rat = e = 0.0
    x = xf
    fx = -func(x, *args)
    num = 1

    ffulc = fnfc = fx
    xm = 0.5 * (a + b)
    tol1 = sqrt_eps * np.abs(xf) + xtol / 3.0
    tol2 = 2.0 * tol1

    while (np.abs(xf - xm) > (tol2 - 0.5 * (b - a))):
        golden = 1
        # Check for parabolic fit
        if np.abs(e) > tol1:
            golden = 0
            r = (xf - nfc) * (fx - ffulc)
            q = (xf - fulc) * (fx - fnfc)
            p = (xf - fulc) * q - (xf - nfc) * r
            q = 2.0 * (q - r)
            if q > 0.0:
                p = -p
            q = np.abs(q)
            r = e
            e = rat

            # Check for acceptability of parabola
            if ((np.abs(p) < np.abs(0.5*q*r)) and (p > q*(a - xf)) and
                    (p < q * (b - xf))):
                rat = (p + 0.0) / q
                x = xf + rat

                if ((x - a) < tol2) or ((b - x) < tol2):
                    si = np.sign(xm - xf) + ((xm - xf) == 0)
                    rat = tol1 * si
            else:      # do a golden section step
                golden = 1

        if golden:  # Do a golden-section step
            if xf >= xm:
                e = a - xf
            else:
                e = b - xf
            rat = golden_mean*e

        if rat == 0:
            si = np.sign(rat) + 1
        else:
            si = np.sign(rat)

        x = xf + si * np.maximum(np.abs(rat), tol1)
        fu = -func(x, *args)
        num += 1

        if fu <= fx:
            if x >= xf:
                a = xf
            else:
                b = xf
            fulc, ffulc = nfc, fnfc
            nfc, fnfc = xf, fx
            xf, fx = x, fu
        else:
            if x < xf:
                a = x
            else:
                b = x
            if (fu <= fnfc) or (nfc == xf):
                fulc, ffulc = nfc, fnfc
                nfc, fnfc = x, fu
            elif (fu <= ffulc) or (fulc == xf) or (fulc == nfc):
                fulc, ffulc = x, fu

        xm = 0.5 * (a + b)
        tol1 = sqrt_eps * np.abs(xf) + xtol / 3.0
        tol2 = 2.0 * tol1

        if num >= maxfun:
            status_flag = 1
            break

    fval = -fx
    info = status_flag, num

    return xf, fval, info

###############################################################################
## https://github.com/linesd/minimize/blob/master/optimizer/minimize.py

# The function uses conjugate gradients and approximate linesearches based 
# on polynomial interpolation with Wolfe-Powel conditions

@njit(cache=True, fastmath=True)
def minimize_wolfe_powel(f, X, length, fargs=(), reduction=None, verbose=False, concise=False):
    """
    Minimize a differentiable multivariate function.
    Parameters
    ----------
    f : function to minimize. The function must return the value
        of the function (float) and a numpy array of partial
        derivatives of shape (D,) with respect to X, where D is
        the dimensionality of the function.
    X : numpy array - Shape : (D, 1)
        initial guess.
    length : int
        The length of the run. If positive, length gives the maximum
        number of line searches, if negative its absolute value gives
        the maximum number of function evaluations.
    args : tuple
        Tuple of parameters to be passed to the function f.
    reduction : float
        The expected reduction in the function value in the first
        linesearch (if None, defaults to 1.0)
    verbose : bool
        If True - prints the progress of minimize. (default is True)
    concise : bool
        If True - returns concise convergence info, only the minimium function
        value (necessary when optimizing a large number of parameters)
        (default is False)
    Return
    ------
    Xs : numpy array - Shape : (D, 1)
        The found solution.
    convergence : numpy array - Shape : (i, D+1)
        Convergence information. The first column is the function values
        returned by the function being minimized. The next D columns are
        the guesses of X during the minimization process.
        If concise = True, convergence information is only the minimum
        function value. This is necessary only when optimizing a large number
        of parameters.
    i : int
        Number of line searches or function evaluations depending on which
        was selected.
    The function returns when either its length is up, or if no further progress
     can be made (ie, we are at a (local) minimum, or so close that due to
     numerical problems, we cannot get any closer)
     Copyright (C) 2001 - 2006 by Carl Edward Rasmussen (2006-09-08).
     Converted to python by David Lines (2019-23-08)
    """
    INT = 0.1  # don't reevaluate within 0.1 of the limit of the current bracket
    EXT = 3.0  # extrapolate maximum 3 times the current step size
    MAX = 20  # max 20 function evaluations per line search
    RATIO = 10  # maximum allowed slope ratio
    SIG = 0.1
    RHO = SIG / 2
    # SIG and RHO control the Wolfe-Powell conditions
    # SIG is the maximum allowed absolute ratio between
    # previous and new slopes (derivatives in the search direction), thus setting
    # SIG to low (positive) values forces higher precision in the line-searches.
    # RHO is the minimum allowed fraction of the expected (from the slope at the
    # initial point in the linesearch). Constants must satisfy 0 < RHO < SIG < 1.
    # Tuning of SIG (depending on the nature of the function to be optimized) may
    # speed up the minimization; it is probably not worth playing much with RHO.

    # print("Minimizing %s ..." % f)

    if reduction is None:
        red = 1.0
    else:
        red = reduction

    S = 'Linesearch' if length > 0 else 'Function evaluation'

    i = 0  # run length counter
    ls_failed = 0  # no previous line search has failed
    f0, df0 = f(X, fargs)  # get initial function value and gradient
    df0 = df0.reshape(-1, 1)
    fX = [];
    fX.append(f0)
    Xd = [];
    Xd.append(X)
    i += (length < 0)  # count epochs
    s = -df0  # get column vec
    d0 = -s.T @ s  # initial search direction (steepest) and slope
    x3 = red / (1 - d0)  # initial step is red/(|s|+1)

    while i < abs(length):  # while not finished
        i += (length > 0)  # count iterations

        X0 = X;
        F0 = f0;
        dF0 = df0  # copy current vals
        M = MAX if length > 0 else min(MAX, -length - i)

        while 1:  # extrapolate as long as necessary
            x2 = 0;
            f2 = f0;
            d2 = d0;
            f3 = f0;
            df3 = df0
            success = False

            while not success and M > 0:
                try:
                    M -= 1;
                    i += (length < 0)  # count epochs
                    f3, df3 = f(X + x3 * s, *list(*fargs))
                    df3 = df3.reshape(-1, 1)
                    if np.isnan(f3) or np.isinf(f3) or np.any(np.isnan(df3) + np.isinf(df3)):
                        raise Exception('Either nan or inf in function eval or gradients')
                    success = True
                except:  # catch any error occuring in f
                    x3 = (x2 + x3) / 2  # bisect and try again

            if f3 < F0:
                X0 = X + x3 * s;
                F0 = f3;
                dF0 = df3  # keep best values

            d3 = df3.T @ s  # new slope
            if d3 > SIG * d0 or f3 > f0 + x3 * RHO * d0 or M == 0:
                break  # finished extrapolating

            x1 = x2;
            f1 = f2;
            d1 = d2  # move point 2 to point 1
            x2 = x3;
            f2 = f3;
            d2 = d3  # move point 3 to point 2
            A = 6 * (f1 - f2) + 3 * (d2 + d1) * (x2 - x1)  # make cubic extrapolation
            B = 3 * (f2 - f1) - (2 * d1 + d2) * (x2 - x1)
            x3 = x1 - d1 * (x2 - x1) ** 2 / (B + np.sqrt(B * B - A * d1 * (x2 - x1)))  # num. error possible, ok!

            if np.iscomplex(x3) or np.isnan(x3) or np.isinf(x3) or x3 < 0:  # num prob | wrong sign
                x3 = x2 * EXT
            elif x3 > x2 * EXT:
                x3 = x2 * EXT
            elif x3 < x2 + INT * (x2 - x1):
                x3 = x2 + INT * (x2 - x1)

        while (abs(d3) > -SIG * d0 or f3 > f0 + x3 * RHO * d0) and M > 0:  # keep interpolating

            if d3 > 0 or f3 > f0 + x3 * RHO * d0:  # choose subinterval
                x4 = x3;
                f4 = f3;
                d4 = d3  # move point 3 to point 4
            else:
                x2 = x3;
                f2 = f3;
                d2 = d3  # move point 3 to point 2

            if f4 > f0:
                x3 = x2 - (0.5 * d2 * (x4 - x2) ** 2) / (f4 - f2 - d2 * (x4 - x2))  # quadratic interpolation
            else:
                A = 6 * (f2 - f4) / (x4 - x2) + 3 * (d4 + d2)  # cubic interpolation
                B = 3 * (f4 - f2) - (2 * d2 + d4) * (x4 - x2)
                x3 = x2 + (np.sqrt(B * B - A * d2 * (x4 - x2) ** 2) - B) / A  # num. error possible, ok!

            if np.isnan(x3) or np.isinf(x3):
                x3 = (x2 + x4) / 2  # if we had a numerical problem then bisect

            x3 = max(min(x3, x4 - INT * (x4 - x2)), x2 + INT * (x4 - x2))  # don't accept too close
            f3, df3 = f(X + x3 * s, *list(fargs))
            df3 = df3.reshape(-1, 1)

            if f3 < F0:
                X0 = X + x3 * s;
                F0 = f3;
                dF0 = df3  # keep best values

            M -= 1;
            i += (length < 0)  # count epochs?!
            d3 = df3.T @ s  # new slope

        if abs(d3) < -SIG * d0 and f3 < f0 + x3 * RHO * d0:  # if line search succeeded
            X = X + x3 * s;
            f0 = f3;
            fX.append(f0);
            Xd.append(X)  # update variables
            if verbose:
                print('%s %6i;  Value %4.6e\r' % (S, i, f0))
            s = (df3.T @ df3 - df0.T @ df3) / (df0.T @ df0) * s - df3  # Polack-Ribiere CG direction
            df0 = df3  # swap derivatives
            d3 = d0;
            d0 = df0.T @ s
            if d0 > 0:  # new slope must be negative
                s = -df0.reshape(-1, 1);
                d0 = -s.T @ s  # otherwise use steepest direction
            x3 = x3 * min(RATIO, d3 / (d0 - np.finfo(np.double).tiny))  # slope ratio but max RATIO
            ls_failed = False  # this line search did not fail
        else:
            X = X0;
            f0 = F0;
            df0 = dF0  # restore best point so far
            if ls_failed or i > abs(length):  # line search failed twice in a row
                break  # or we ran out of time, so we give up
            s = -df0.reshape(-1, 1);
            d0 = -s.T @ s  # try steepest
            x3 = 1 / (1 - d0)
            ls_failed = True  # this line search failed

    if concise:
        convergence = fX[-1]  # return only the minimum function value
    else:
        convergence = np.hstack((np.array(fX).reshape(-1, 1), np.array(Xd)[:, :, 0]))  # bundle convergence info

    Xs = X  # solution

    return Xs, convergence, i

###############################################################################
