"""Hagan-West monotone-convex interpolation (Hagan and West, 2008, §6)."""
import numpy as np
from ...utils.error import FinError


def _section(g0, g1, x):
    """Normalized instantaneous-forward correction and its integral."""
    if abs(g0) + abs(g1) < 1e-14:
        return 0., 0.
    a, b = g1 + 2*g0, g0 + 2*g1
    if (a < 0 <= b) or (a > 0 >= b):
        return ((3*(g0+g1)*x-2*a)*x+g0,
                x*(1-x)*(g0*(1-x)-g1*x))
    if (g0 < 0 <= a) or (g0 > 0 >= a):
        eta = a/(g1-g0)
        z = max(0., (x-eta)/(1-eta))
        return (g0+(g1-g0)*z*z, g0*x+(g1-g0)*(1-eta)*z**3/3)
    if (g1 <= 0 < b) or (g1 >= 0 > b):
        eta = 3*g1/(g1-g0)
        z = max(0., (eta-min(x,eta))/eta) if eta > 0 else 0.
        return (g1+(g0-g1)*z*z, g1*x+eta*(g0-g1)*(1-z**3)/3)
    eta = g1/(g0+g1)
    middle = -.5*(eta*g0+(1-eta)*g1)
    zl = (eta-min(x,eta))/eta if eta > 0 else 0.
    zr = (max(x,eta)-eta)/(1-eta) if eta < 1 else 0.
    return (middle+(g0-middle)*zl**2+(g1-middle)*zr**2,
            middle*x+(g0-middle)*eta*(1-zl**3)/3
            +(g1-middle)*(1-eta)*zr**3/3)


class HaganWestMonotoneConvex:
    """Build exact-knot discount factors from continuous interpolated forwards.

    Input times include zero, and the first discount factor equals one.
    Negative rates are permitted; extrapolation uses a flat terminal
    instantaneous forward.
    """

    def __init__(self, times, dfs):
        t=np.asarray(times,dtype=float)
        d=np.asarray(dfs,dtype=float)
        if (t.ndim != 1 or d.shape != t.shape or t.size < 2
                or not np.isfinite(t).all() or not np.isfinite(d).all()
                or t[0] != 0 or abs(d[0]-1) > 1e-10
                or np.any(np.diff(t)<=0) or np.any(d<=0)):
            raise FinError("Need strictly increasing [0,...] times and positive discount factors (df(0)=1).")
        self.times=t
        self.dfs=d
        self.lengths=np.diff(t)
        self.discrete=-np.diff(np.log(d))/self.lengths
        n=len(self.discrete)
        self.nodes=np.empty(n+1)
        if n==1:
            self.nodes[:]=self.discrete[0]
        else:
            self.nodes[1:-1]=(self.lengths[:-1]*self.discrete[1:]
                              +self.lengths[1:]*self.discrete[:-1])/(self.lengths[:-1]+self.lengths[1:])
            self.nodes[0]=self.discrete[0]-.5*(self.nodes[1]-self.discrete[0])
            self.nodes[-1]=self.discrete[-1]-.5*(self.nodes[-2]-self.discrete[-1])

    def _single(self, time, is_forward):
        if not np.isfinite(time) or time<0:
            raise FinError("Interpolation time must be finite and non-negative.")
        if time==0:
            return float(self.nodes[0]) if is_forward else 1.
        if time>=self.times[-1]:
            if is_forward:
                return float(self.nodes[-1])
            return float(self.dfs[-1]*np.exp(-self.nodes[-1]*(time-self.times[-1])))
        i=max(0,int(np.searchsorted(self.times,time))-1)
        x=(time-self.times[i])/self.lengths[i]
        fd=self.discrete[i]
        g, integral=_section(self.nodes[i]-fd,self.nodes[i+1]-fd,x)
        if is_forward:
            return float(fd+g)
        return float(self.dfs[i]*np.exp(-self.lengths[i]*(fd*x+integral)))

    def _evaluate(self, time, is_forward):
        v=np.asarray(time,dtype=float)
        if v.ndim==0:
            return self._single(float(v),is_forward)
        return np.array([self._single(float(t),is_forward) for t in v.flat]).reshape(v.shape)

    def forward_rate(self, time):
        """Continuous instantaneous forward rate."""
        return self._evaluate(time,True)

    def discount_factor(self,time):
        """Discount factor consistent with market DF knots."""
        return self._evaluate(time,False)

