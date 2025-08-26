import warnings
import numpy as np
from scipy.optimize.linesearch import (
    LineSearchWarning,
    line_search_wolfe1,
    line_search_wolfe2,
)

Inf = np.inf
EPSILON = 1e-20

_status_message = {
    "success": "Optimization terminated successfully.",
    "maxiter": "Maximum number of iterations exceeded.",
    "pr_loss": "Desired error not necessarily achieved due to precision loss.",
    "nan": "NaN result encountered.",
}


class _LineSearchError(RuntimeError):
    """Line search failed to find a suitable step length."""

    pass


class ScalarFunction:
    """Minimal wrapper for function, gradient, and Hessian caching."""

    def __init__(
        self,
        fun,
        x0,
        args=(),
        grad=None,
        hess=None,
        epsilon=1.4901161193847656e-8,
        bounds=None,
        finite_diff_rel_step=None,
    ):
        self.fun = lambda x: fun(x, *args)
        if callable(grad):
            self.grad = lambda x: grad(x, *args)
        elif grad == "2-point":
            self.grad = self._finite_diff_grad
        else:
            raise ValueError("Unsupported grad type")
        self.hess = hess or (lambda x: None)
        self.nfev = 0
        self.ngev = 0
        self.epsilon = epsilon

    def __call__(self, x):
        self.nfev += 1
        return self.fun(x)

    def _finite_diff_grad(self, x):
        self.ngev += 1
        grad = np.zeros_like(x, dtype=float)
        fx = self.fun(x)
        eps = self.epsilon
        for i in range(len(x)):
            x_eps = np.array(x, dtype=float)
            x_eps[i] += eps
            grad[i] = (self.fun(x_eps) - fx) / eps
        return grad


def vecnorm(x, ord=2):
    if ord == Inf:
        return np.amax(np.abs(x))
    elif ord == -Inf:
        return np.amin(np.abs(x))
    else:
        return np.sum(np.abs(x) ** ord, axis=0) ** (1.0 / ord)


class OptimizeResult(dict):
    """Represents the optimization result."""

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
            return "\n".join(
                [k.rjust(m) + ": " + repr(v) for k, v in sorted(self.items())]
            )
        else:
            return self.__class__.__name__ + "()"

    def __dir__(self):
        return list(self.keys())


def _prepare_scalar_function(fun, x0, jac=None, args=(), epsilon=None):
    """Prepare scalar function for optimization."""
    grad = jac if callable(jac) else "2-point"
    if epsilon is None:
        epsilon = EPSILON
    sf = ScalarFunction(fun, x0, args=args, grad=grad, epsilon=epsilon)
    return sf


def _line_search_wolfe12(f, fprime, xk, pk, gfk, old_fval, old_old_fval, **kwargs):
    """Wolfe line search with fallback to line_search_wolfe2."""
    extra_condition = kwargs.pop("extra_condition", None)
    ret = line_search_wolfe1(f, fprime, xk, pk, gfk, old_fval, old_old_fval, **kwargs)

    if ret[0] is not None and extra_condition is not None:
        xp1 = xk + ret[0] * pk
        if not extra_condition(ret[0], xp1, ret[3], ret[5]):
            ret = (None,)

    if ret[0] is None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", LineSearchWarning)
            kwargs2 = {k: kwargs[k] for k in ("c1", "c2", "amax") if k in kwargs}
            ret = line_search_wolfe2(
                f,
                fprime,
                xk,
                pk,
                gfk,
                old_fval,
                old_old_fval,
                extra_condition=extra_condition,
                **kwargs2,
            )

    if ret[0] is None:
        raise _LineSearchError()
    return ret


def _minimize_cg(
    fun,
    x0,
    args=(),
    jac=None,
    callback=None,
    gtol=1e-5,
    norm=Inf,
    eps=EPSILON,
    maxiter=None,
    disp=False,
    return_all=False,
    **unknown_options,
):

    x0 = np.asarray(x0).flatten()
    if maxiter is None:
        maxiter = len(x0) * 200

    sf = _prepare_scalar_function(fun, x0, jac=jac, args=args, epsilon=eps)

    f = sf.fun
    myfprime = sf.grad

    old_fval = f(x0)
    gfk = myfprime(x0)
    if not np.isscalar(old_fval):
        try:
            old_fval = old_fval.item()
        except Exception:
            raise ValueError("Objective function must return a scalar.")

    k = 0
    xk = x0
    old_old_fval = old_fval + np.linalg.norm(gfk) / 2
    pk = -gfk
    gnorm = vecnorm(gfk, ord=norm)

    if return_all:
        allvecs = [xk]
    warnflag = 0
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
            return alpha, xkp1, pkp1, gfkp1, gnorm

        def descent_condition(alpha, xkp1, fp1, gfkp1):
            cached_step[:] = polak_ribiere_powell_step(alpha, gfkp1)
            alpha, xk, pk, gfk, gnorm = cached_step
            return gnorm <= gtol or np.dot(pk, gfk) <= -sigma_3 * np.dot(gfk, gfk)

        try:
            alpha_k, fc, gc, old_fval, old_old_fval, gfkp1 = _line_search_wolfe12(
                f,
                myfprime,
                xk,
                pk,
                gfk,
                old_fval,
                old_old_fval,
                c2=0.4,
                amin=1e-100,
                amax=1e100,
                extra_condition=descent_condition,
            )
        except _LineSearchError:
            warnflag = 2
            break

        if alpha_k == cached_step[0]:
            alpha_k, xk, pk, gfk, gnorm = cached_step
        else:
            alpha_k, xk, pk, gfk, gnorm = polak_ribiere_powell_step(alpha_k, gfkp1)

        if return_all:
            allvecs.append(xk)
        if callback is not None:
            callback(xk)
        k += 1

    fval = old_fval
    if warnflag == 2:
        msg = _status_message["pr_loss"]
    elif k >= maxiter:
        warnflag = 1
        msg = _status_message["maxiter"]
    elif np.isnan(gnorm) or np.isnan(fval) or np.isnan(xk).any():
        warnflag = 3
        msg = _status_message["nan"]
    else:
        msg = _status_message["success"]

    if disp:
        prefix = "Warning: " if warnflag != 0 else ""
        print(f"{prefix}{msg}")
        print(f"         Current function value: {fval:.6f}")
        print(f"         Iterations: {k}")
        print(f"         Function evaluations: {sf.nfev}")
        print(f"         Gradient evaluations: {sf.ngev}")

    result = OptimizeResult(
        fun=fval,
        jac=gfk,
        nfev=sf.nfev,
        njev=sf.ngev,
        status=warnflag,
        success=(warnflag == 0),
        message=msg,
        x=xk,
        nit=k,
    )
    if return_all:
        result["allvecs"] = allvecs
    return result


def fmin_cg(
    f,
    x0,
    fprime=None,
    fargs=(),
    gtol=1e-5,
    norm=Inf,
    epsilon=EPSILON,
    maxiter=None,
    full_output=0,
    disp=1,
    retall=0,
    callback=None,
):
    """Minimize a function using nonlinear conjugate gradient (Polak-Ribiere)."""
    opts = {
        "gtol": gtol,
        "norm": norm,
        "eps": epsilon,
        "disp": disp,
        "maxiter": maxiter,
        "return_all": retall,
    }
    res = _minimize_cg(f, x0, args=fargs, jac=fprime, callback=callback, **opts)

    if full_output:
        retlist = res["x"], res["fun"], res["nfev"], res["njev"], res["status"]
        if retall:
            retlist += (res["allvecs"],)
        return retlist
    else:
        if retall:
            return res["x"], res["allvecs"]

        return res["x"]
