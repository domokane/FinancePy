###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

# THIS IS IN PROGRESS
from numba import njit, boolean, int64, float64, vectorize
import numpy as np

# ******NOTICE***************
# optimize.py module by Travis E. Oliphant
#
# You may copy and use this module as you see fit with no
# guarantee implied provided you keep this notice in all copies.
# *****END NOTICE************

# A collection of optimization algorithms. Version 0.5
# CHANGES
#  Added fminbound (July 2001)
#  Added brute (Aug. 2002)
#  Finished line search satisfying strong Wolfe conditions (Mar. 2004)
#  Updated strong Wolfe conditions line search to use
#  cubic-interpolation (Mar. 2004)

###############################################################################

Inf = np.inf
_epsilon = 1e-20

def fmin_cg(f, x0, fprime=None, fargs=(), gtol=1e-5, norm=Inf, epsilon=_epsilon,
            maxiter=None, full_output=0, disp=1, retall=0, callback=None):
    """
    Minimize a function using a nonlinear conjugate gradient algorithm.
    Parameters
    ----------
    f : callable, ``f(x, *args)``
        Objective function to be minimized. Here `x` must be a 1-D array of
        the variables that are to be changed in the search for a minimum, and
        `args` are the other (fixed) parameters of `f`.
    x0 : ndarray
        A user-supplied initial estimate of `xopt`, the optimal value of `x`.
        It must be a 1-D array of values.
    fprime : callable, ``fprime(x, *args)``, optional
        A function that returns the gradient of `f` at `x`. Here `x` and `args`
        are as described above for `f`. The returned value must be a 1-D array.
        Defaults to None, in which case the gradient is approximated
        numerically (see `epsilon`, below).
    args : tuple, optional
        Parameter values passed to `f` and `fprime`. Must be supplied whenever
        additional fixed parameters are needed to completely specify the
        functions `f` and `fprime`.
    gtol : float, optional
        Stop when the norm of the gradient is less than `gtol`.
    norm : float, optional
        Order to use for the norm of the gradient
        (``-np.Inf`` is min, ``np.Inf`` is max).
    epsilon : float or ndarray, optional
        Step size(s) to use when `fprime` is approximated numerically. Can be a
        scalar or a 1-D array. Defaults to ``sqrt(eps)``, with eps the
        floating point machine precision.  Usually ``sqrt(eps)`` is about
        1.5e-8.
    maxiter : int, optional
        Maximum number of iterations to perform. Default is ``200 * len(x0)``.
    full_output : bool, optional
        If True, return `fopt`, `func_calls`, `grad_calls`, and `warnflag` in
        addition to `xopt`.  See the Returns section below for additional
        information on optional return values.
    disp : bool, optional
        If True, return a convergence message, followed by `xopt`.
    retall : bool, optional
        If True, add to the returned values the results of each iteration.
    callback : callable, optional
        An optional user-supplied function, called after each iteration.
        Called as ``callback(xk)``, where ``xk`` is the current value of `x0`.
    Returns
    -------
    xopt : ndarray
        Parameters which minimize f, i.e., ``f(xopt) == fopt``.
    fopt : float, optional
        Minimum value found, f(xopt). Only returned if `full_output` is True.
    func_calls : int, optional
        The number of function_calls made. Only returned if `full_output`
        is True.
    grad_calls : int, optional
        The number of gradient calls made. Only returned if `full_output` is
        True.
    warnflag : int, optional
        Integer value with warning status, only returned if `full_output` is
        True.
        0 : Success.
        1 : The maximum number of iterations was exceeded.
        2 : Gradient and/or function calls were not changing. May indicate
            that precision was lost, i.e., the routine did not converge.
        3 : NaN result encountered.
    allvecs : list of ndarray, optional
        List of arrays, containing the results at each iteration.
        Only returned if `retall` is True.
    See Also
    --------
    minimize : common interface to all `scipy.optimize` algorithms for
               unconstrained and constrained minimization of multivariate
               functions. It provides an alternative way to call
               ``fmin_cg``, by specifying ``method='CG'``.
    Notes
    -----
    This conjugate gradient algorithm is based on that of Polak and Ribiere
    [1]_.
    Conjugate gradient methods tend to work better when:
    1. `f` has a unique global minimizing point, and no local minima or
       other stationary points,
    2. `f` is, at least locally, reasonably well approximated by a
       quadratic function of the variables,
    3. `f` is continuous and has a continuous gradient,
    4. `fprime` is not too large, e.g., has a norm less than 1000,
    5. The initial guess, `x0`, is reasonably close to `f` 's global
       minimizing point, `xopt`.
    References
    ----------
    .. [1] Wright & Nocedal, "Numerical Optimization", 1999, pp. 120-122.
    Examples
    --------
    Example 1: seek the minimum value of the expression
    ``a*u**2 + b*u*v + c*v**2 + d*u + e*v + f`` for given values
    of the parameters and an initial guess ``(u, v) = (0, 0)``.
    >>> args = (2, 3, 7, 8, 9, 10)  # parameter values
    >>> def f(x, *args):
    ...     u, v = x
    ...     a, b, c, d, e, f = args
    ...     return a*u**2 + b*u*v + c*v**2 + d*u + e*v + f
    >>> def gradf(x, *args):
    ...     u, v = x
    ...     a, b, c, d, e, f = args
    ...     gu = 2*a*u + b*v + d     # u-component of the gradient
    ...     gv = b*u + 2*c*v + e     # v-component of the gradient
    ...     return np.asarray((gu, gv))
    >>> x0 = np.asarray((0, 0))  # Initial guess.
    >>> from scipy import optimize
    >>> res1 = optimize.fmin_cg(f, x0, fprime=gradf, args=args)
    Optimization terminated successfully.
             Current function value: 1.617021
             Iterations: 4
             Function evaluations: 8
             Gradient evaluations: 8
    >>> res1
    array([-1.80851064, -0.25531915])
    Example 2: solve the same problem using the `minimize` function.
    (This `myopts` dictionary shows all of the available options,
    although in practice only non-default values would be needed.
    The returned value will be a dictionary.)
    >>> opts = {'maxiter' : None,    # default value.
    ...         'disp' : True,    # non-default value.
    ...         'gtol' : 1e-5,    # default value.
    ...         'norm' : np.inf,  # default value.
    ...         'eps' : 1.4901161193847656e-08}  # default value.
    >>> res2 = optimize.minimize(f, x0, jac=gradf, args=args,
    ...                          method='CG', options=opts)
    Optimization terminated successfully.
            Current function value: 1.617021
            Iterations: 4
            Function evaluations: 8
            Gradient evaluations: 8
    >>> res2.x  # minimum found
    array([-1.80851064, -0.25531915])
    """
    opts = {'gtol': gtol,
            'norm': norm,
            'eps': epsilon,
            'disp': disp,
            'maxiter': maxiter,
            'return_all': retall}

    res = _minimize_cg(f, x0, fargs, fprime, callback=callback, **opts)

    if full_output:
        retlist = res['x'], res['fun'], res['nfev'], res['njev'], res['status']
        if retall:
            retlist += (res['allvecs'], )
        return retlist
    else:
        if retall:
            return res['x'], res['allvecs']
        else:
            return res['x']

def _check_unknown_options(unknown_options):
    if unknown_options:
        msg = ", ".join(map(str, unknown_options.keys()))
        # Stack level 4: this is called from _minimize_*, which is
        # called from another function in SciPy. Level 4 is the first
        # level in user code.
        print("Unknown solver options")
        

from scipy.optimize._differentiable_functions import ScalarFunction

def _prepare_scalar_function(fun, x0, jac=None, args=(), bounds=None,
                             epsilon=None, finite_diff_rel_step=None,
                             hess=None):
    """
    Creates a ScalarFunction object for use with scalar minimizers
    (BFGS/LBFGSB/SLSQP/TNC/CG/etc).
    Parameters
    ----------
    fun : callable
        The objective function to be minimized.
            ``fun(x, *args) -> float``
        where ``x`` is an 1-D array with shape (n,) and ``args``
        is a tuple of the fixed parameters needed to completely
        specify the function.
    x0 : ndarray, shape (n,)
        Initial guess. Array of real elements of size (n,),
        where 'n' is the number of independent variables.
    jac : {callable,  '2-point', '3-point', 'cs', None}, optional
        Method for computing the gradient vector. If it is a callable, it
        should be a function that returns the gradient vector:
            ``jac(x, *args) -> array_like, shape (n,)``
        If one of `{'2-point', '3-point', 'cs'}` is selected then the gradient
        is calculated with a relative step for finite differences. If `None`,
        then two-point finite differences with an absolute step is used.
    args : tuple, optional
        Extra arguments passed to the objective function and its
        derivatives (`fun`, `jac` functions).
    bounds : sequence, optional
        Bounds on variables. 'new-style' bounds are required.
    eps : float or ndarray
        If `jac is None` the absolute step size used for numerical
        approximation of the jacobian via forward differences.
    finite_diff_rel_step : None or array_like, optional
        If `jac in ['2-point', '3-point', 'cs']` the relative step size to
        use for numerical approximation of the jacobian. The absolute step
        size is computed as ``h = rel_step * sign(x0) * max(1, abs(x0))``,
        possibly adjusted to fit into the bounds. For ``method='3-point'``
        the sign of `h` is ignored. If None (default) then step is selected
        automatically.
    hess : {callable,  '2-point', '3-point', 'cs', None}
        Computes the Hessian matrix. If it is callable, it should return the
        Hessian matrix:
            ``hess(x, *args) -> {LinearOperator, spmatrix, array}, (n, n)``
        Alternatively, the keywords {'2-point', '3-point', 'cs'} select a
        finite difference scheme for numerical estimation.
        Whenever the gradient is estimated via finite-differences, the Hessian
        cannot be estimated with options {'2-point', '3-point', 'cs'} and needs
        to be estimated using one of the quasi-Newton strategies.
    Returns
    -------
    sf : ScalarFunction
    """
    if callable(jac):
        grad = jac
#    elif jac in FD_METHODS:
#        # epsilon is set to None so that ScalarFunction is made to use
#        # rel_step
#        epsilon = None
#        grad = jac
    else:
        # default (jac is None) is to do 2-point finite differences with
        # absolute step size. ScalarFunction has to be provided an
        # epsilon value that is not None to use absolute steps. This is
        # normally the case from most _minimize* methods.
        grad = '2-point'
        epsilon = epsilon

    if hess is None:
        # ScalarFunction requires something for hess, so we give a dummy
        # implementation here if nothing is provided, return a value of None
        # so that downstream minimisers halt. The results of `fun.hess`
        # should not be used.
        def hess(x, *args):
            return None

    if bounds is None:
        bounds = (-np.inf, np.inf)

    # ScalarFunction caches. Reuse of fun(x) during grad
    # calculation reduces overall function evaluations.
    sf = ScalarFunction(fun, x0, args, grad, hess,
                        finite_diff_rel_step, bounds, epsilon=epsilon)

    return sf


def vecnorm(x, ord=2):
    if ord == Inf:
        return np.amax(np.abs(x))
    elif ord == -Inf:
        return np.amin(np.abs(x))
    else:
        return np.sum(np.abs(x)**ord, axis=0)**(1.0 / ord)

class OptimizeResult(dict):
    """ Represents the optimization result.
    Attributes
    ----------
    x : ndarray
        The solution of the optimization.
    success : bool
        Whether or not the optimizer exited successfully.
    status : int
        Termination status of the optimizer. Its value depends on the
        underlying solver. Refer to `message` for details.
    message : str
        Description of the cause of the termination.
    fun, jac, hess: ndarray
        Values of objective function, its Jacobian and its Hessian (if
        available). The Hessians may be approximations, see the documentation
        of the function in question.
    hess_inv : object
        Inverse of the objective function's Hessian; may be an approximation.
        Not available for all solvers. The type of this attribute may be
        either np.ndarray or scipy.sparse.linalg.LinearOperator.
    nfev, njev, nhev : int
        Number of evaluations of the objective functions and of its
        Jacobian and Hessian.
    nit : int
        Number of iterations performed by the optimizer.
    maxcv : float
        The maximum constraint violation.
    Notes
    -----
    There may be additional attributes not listed above depending of the
    specific solver. Since this class is essentially a subclass of dict
    with attribute accessors, one can see which attributes are available
    using the `keys()` method.
    """

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError as e:
            raise AttributeError(name) from e

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __repr__(self):
        if self.keys():
            m = max(map(len, list(self.keys()))) + 1
            return '\n'.join([k.rjust(m) + ': ' + repr(v)
                              for k, v in sorted(self.items())])
        else:
            return self.__class__.__name__ + "()"

    def __dir__(self):
        return list(self.keys())


def _line_search_wolfe12(f, fprime, xk, pk, gfk, old_fval, old_old_fval,
                         **kwargs):
    """
    Same as line_search_wolfe1, but fall back to line_search_wolfe2 if
    suitable step length is not found, and raise an exception if a
    suitable step length is not found.
    Raises
    ------
    _LineSearchError
        If no suitable step size is found
    """

    extra_condition = kwargs.pop('extra_condition', None)

    ret = line_search_wolfe1(f, fprime, xk, pk, gfk,
                             old_fval, old_old_fval,
                             **kwargs)

    if ret[0] is not None and extra_condition is not None:
        xp1 = xk + ret[0] * pk
        if not extra_condition(ret[0], xp1, ret[3], ret[5]):
            # Reject step if extra_condition fails
            ret = (None,)

    if ret[0] is None:
        # line search failed: try different one.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', LineSearchWarning)
            kwargs2 = {}
            for key in ('c1', 'c2', 'amax'):
                if key in kwargs:
                    kwargs2[key] = kwargs[key]
            ret = line_search_wolfe2(f, fprime, xk, pk, gfk,
                                     old_fval, old_old_fval,
                                     extra_condition=extra_condition,
                                     **kwargs2)

    if ret[0] is None:
        raise _LineSearchError()

    return ret

def _minimize_cg(fun, x0, args=(), jac=None, callback=None,
                 gtol=1e-5, norm=Inf, eps=_epsilon, maxiter=None,
                 disp=False, return_all=False, finite_diff_rel_step=None,
                 **unknown_options):
    """
    Minimization of scalar function of one or more variables using the
    conjugate gradient algorithm.
    Options
    -------
    disp : bool
        Set to True to print convergence messages.
    maxiter : int
        Maximum number of iterations to perform.
    gtol : float
        Gradient norm must be less than `gtol` before successful
        termination.
    norm : float
        Order of norm (Inf is max, -Inf is min).
    eps : float or ndarray
        If `jac is None` the absolute step size used for numerical
        approximation of the jacobian via forward differences.
    return_all : bool, optional
        Set to True to return a list of the best solution at each of the
        iterations.
    finite_diff_rel_step : None or array_like, optional
        If `jac in ['2-point', '3-point', 'cs']` the relative step size to
        use for numerical approximation of the jacobian. The absolute step
        size is computed as ``h = rel_step * sign(x0) * max(1, abs(x0))``,
        possibly adjusted to fit into the bounds. For ``method='3-point'``
        the sign of `h` is ignored. If None (default) then step is selected
        automatically.
    """
    _check_unknown_options(unknown_options)

    retall = return_all

    x0 = np.asarray(x0).flatten()
    if maxiter is None:
        maxiter = len(x0) * 200

    sf = _prepare_scalar_function(fun, x0, jac=jac, args=args, epsilon=eps,
                                  finite_diff_rel_step=finite_diff_rel_step)

    f = sf.fun
    myfprime = sf.grad

    old_fval = f(x0)
    gfk = myfprime(x0)

    if not np.isscalar(old_fval):
        try:
            old_fval = old_fval.item()
        except (ValueError, AttributeError):
            raise ValueError("The user-provided "
                             "objective function must "
                             "return a scalar value.")

    k = 0
    xk = x0
    # Sets the initial step guess to dx ~ 1
    old_old_fval = old_fval + np.linalg.norm(gfk) / 2

    if retall:
        allvecs = [xk]
    warnflag = 0
    pk = -gfk
    gnorm = vecnorm(gfk, ord=norm)

    sigma_3 = 0.01

    while (gnorm > gtol) and (k < maxiter):
        deltak = np.dot(gfk, gfk)

        cached_step = [None]

        def polak_ribiere_powell_step(alpha, gfkp1=None):
            xkp1 = xk + alpha * pk
            if gfkp1 is None:
                gfkp1 = myfprime(xkp1)
            yk = gfkp1 - gfk
            beta_k = max(0, np.dot(yk, gfkp1) / deltak)
            pkp1 = -gfkp1 + beta_k * pk
            gnorm = vecnorm(gfkp1, ord=norm)
            return (alpha, xkp1, pkp1, gfkp1, gnorm)

        def descent_condition(alpha, xkp1, fp1, gfkp1):
            # Polak-Ribiere+ needs an explicit check of a sufficient
            # descent condition, which is not guaranteed by strong Wolfe.
            #
            # See Gilbert & Nocedal, "Global convergence properties of
            # conjugate gradient methods for optimization",
            # SIAM J. Optimization 2, 21 (1992).
            cached_step[:] = polak_ribiere_powell_step(alpha, gfkp1)
            alpha, xk, pk, gfk, gnorm = cached_step

            # Accept step if it leads to convergence.
            if gnorm <= gtol:
                return True

            # Accept step if sufficient descent condition applies.
            return np.dot(pk, gfk) <= -sigma_3 * np.dot(gfk, gfk)

        try:
            alpha_k, fc, gc, old_fval, old_old_fval, gfkp1 = \
                     _line_search_wolfe12(f, myfprime, xk, pk, gfk, old_fval,
                                          old_old_fval, c2=0.4, amin=1e-100, amax=1e100,
                                          extra_condition=descent_condition)
        except _LineSearchError:
            # Line search failed to find a better solution.
            warnflag = 2
            break

        # Reuse already computed results if possible
        if alpha_k == cached_step[0]:
            alpha_k, xk, pk, gfk, gnorm = cached_step
        else:
            alpha_k, xk, pk, gfk, gnorm = polak_ribiere_powell_step(alpha_k, gfkp1)

        if retall:
            allvecs.append(xk)
        if callback is not None:
            callback(xk)
        k += 1

    fval = old_fval
    if warnflag == 2:
        msg = _status_message['pr_loss']
    elif k >= maxiter:
        warnflag = 1
        msg = _status_message['maxiter']
    elif np.isnan(gnorm) or np.isnan(fval) or np.isnan(xk).any():
        warnflag = 3
        msg = _status_message['nan']
    else:
        msg = _status_message['success']

    if disp:
        print("%s%s" % ("Warning: " if warnflag != 0 else "", msg))
        print("         Current function value: %f" % fval)
        print("         Iterations: %d" % k)
        print("         Function evaluations: %d" % sf.nfev)
        print("         Gradient evaluations: %d" % sf.ngev)

    result = OptimizeResult(fun=fval, jac=gfk, nfev=sf.nfev,
                            njev=sf.ngev, status=warnflag,
                            success=(warnflag == 0), message=msg, x=xk,
                            nit=k)
    if retall:
        result['allvecs'] = allvecs
    return result
