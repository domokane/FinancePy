##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


### DO NOT USE!!!!!!!!!!!!!!!!!

from numba import njit, float64, int64
from scipy import exp, log
from numpy import sqrt
import numpy as np

from finutils.FinMath import norminvcdf
from FinProcess import FinProcess
from ...finutils.FinHelperFunctions import labelToString

################################################################################
# Heston Process
# dS = rS dt + sqrt(V) * S * dz
# dV = kappa(theta-V) dt + sigma sqrt(V) dz
# corr(dV,dS) = rho dt
# Rewritten as
# dS = rS dt + sqrt(V) * S * (rhohat dz1 + rho dz2)
# dV = kappa(theta-V) dt + sigma sqrt(V) dz2
# where rhohat = sqrt(1-rho*rho)
################################################################################
# TODO - DECIDE WHETHER TO OO MODEL
# TODO - NEEDS CHECKING FOR MC CONVERGENCE
################################################################################

class FinHestonProcess(FinProcess):
    
    def getPathsAssets(self, 
                 t,
                 mus,
                 stockPrices,
                 volatilities,
                 betas,
                 seed, 
                 fast = FinFastNumericalApproach.NUMBA):

          if self._numTimeSteps == 1:
              paths = getGBMAssetsNUMBA(self._numAssets,
                                             self._numPaths,
                                             t,mus,stockPrices,
                                             volatilities,betas,seed)
          else:
              paths = getGBMAssetsPathsNUMBA(self._numAssets,
                                             self._numPaths,
                                             self._numTimeSteps,
                                             t,mus,stockPrices,
                                             volatilities,betas,seed)
              
          return paths


from enum import Enum
class FinHestonScheme(Enum):
    EULER = 1
    EULERLOG = 2
    QUADEXP = 3

################################################################################

@njit(float64[:,:](float64,float64,float64,float64,float64,float64,float64,
      float64,float64,float64,int64,int64,int64),fastmath=True)
def getPaths(s0,r,q,v0,kappa,theta,sigma,rho,t,dt,numPaths,seed,scheme):

    np.random.seed(seed)
    numSteps = int(t/dt)
    sPaths = np.zeros(shape=(numPaths,numSteps))
    sPaths[:,0] = s0
    sdt = sqrt(dt)
    rhohat = sqrt(1.0 - rho*rho)
    sigma2 = sigma*sigma

    if scheme == FinHestonNumericalScheme.EULER.value:
        # Basic scheme to first order with truncation on variance
        for iPath in range(0,numPaths):
            s = s0
            v = v0
            for iStep in range(1,numSteps):
                z1 = np.random.normal(0.0,1.0)*sdt
                z2 = np.random.normal(0.0,1.0)*sdt
                zV = z1
                zS = rho*z1 + rhohat*z2
                vplus = max(v,0.0)
                rtvplus = sqrt(vplus)
                v += kappa*(theta-vplus)*dt + sigma*rtvplus*zV + 0.25*sigma2*(zV*zV-dt)
                s += (r-q)*s*dt + rtvplus*s*zS + 0.5*s*vplus*(zV*zV-dt)
                sPaths[iPath,iStep] = s

    elif scheme == FinHestonNumericalScheme.EULERLOG.value:
        # Basic scheme to first order with truncation on variance
        for iPath in range(0,numPaths):
            x = log(s0)
            v = v0
            for iStep in range(1,numSteps):
                zV = np.random.normal(0.0,1.0)*sdt
                zS = rho*zV + rhohat*np.random.normal(0.0,1.0)*sdt
                vplus = max(v,0.0)
                rtvplus = sqrt(vplus)
                x += (r-q-0.5*vplus)*dt + rtvplus*zS
                v += kappa*(theta-vplus)*dt + sigma*rtvplus*zV + sigma2*(zV*zV-dt)/4.0
                sPaths[iPath,iStep] = exp(x)

    elif scheme == FinHestonNumericalScheme.QUADEXP.value:
        # Due to Leif Andersen(2006) 
        Q = exp(-kappa*dt)
        psic = 1.50
        gamma1 = 0.50
        gamma2 = 0.50
        K0 = -rho*kappa*theta*dt/sigma
        K1 = gamma1*dt*(kappa*rho/sigma-0.5)-rho/sigma
        K2 = gamma2*dt*(kappa*rho/sigma-0.5)+rho/sigma
        K3 = gamma1*dt*(1.0-rho*rho)
        K4 = gamma2*dt*(1.0-rho*rho)
        A = K2 + 0.5*K4
        mu = (r-q)
        c1 = sigma2*Q*(1.0-Q)/kappa
        c2 = theta*sigma2*((1.0-Q)**2)/2.0/kappa

        for iPath in range(0,numPaths):
            x = log(s0)
            vn = v0            
            for iStep in range(1,numSteps):
                zV = np.random.normal(0,1)
                zS = rho*zV + rhohat*np.random.normal(0,1)
                m = theta + (vn-theta)*Q
                s2 = c1*vn + c2
                psi = s2/(m*m)
                u = np.random.uniform(0.0,1.0)

                if psi <= psic:
                    b2 = 2.0/psi - 1.0 + sqrt((2.0/psi)*(2.0/psi-1.0))
                    a = m/(1.0+b2)
                    b = sqrt(b2)
                    zV = norminvcdf(u)
                    vnp = a*((b+zV)**2)
                    d = (1.0-2.0*A*a)
                    M = exp((A*b2*a)/d)/sqrt(d)
                    K0 = -log(M) - (K1+0.5*K3)*vn
                else:
                    p = (psi-1.0)/(psi+1.0)
                    beta = (1.0-p)/m
                    
                    if u <= p:
                        vnp = 0.0
                    else:
                        vnp = log((1.0-p)/(1.0-u))/beta

                    M = p + beta*(1.0-p)/(beta-A)
                    K0 = -log(M) - (K1+0.5*K3)*vn
                    
                x += mu*dt + K0 + (K1*vn + K2*vnp) + sqrt(K3*vn+K4*vnp)*zS
                sPaths[iPath,iStep] = exp(x)
                vn = vnp
    else:
        raise FinError("Unknown FinHestonNumericalSchme")

    return sPaths

################################################################################