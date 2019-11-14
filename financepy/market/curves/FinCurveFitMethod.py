import numpy as np
from scipy.interpolate import splev


###############################################################################


class FinCurveFitMethod():
    pass

###############################################################################


class FinCurveFitMethodPolynomial():

    def __init__(self, power=3):
        self._parentType = FinCurveFitMethod
        self._power = power
        self._coeffs = []

    def interpolatedYield(self, t):
        yld = np.polyval(self._coeffs, t)
        return yld

###############################################################################


class FinCurveFitMethodNelsonSiegel():

    def __init__(self, tau=None):
        self._parentType = FinCurveFitMethod
        self._beta1 = None
        self._beta2 = None
        self._beta3 = None
        self._tau = tau

    def interpolatedYield(self, t,
                          beta1=None, beta2=None, beta3=None, tau=None):

        t = np.maximum(t, 1e-10)

        if beta1 is None:
            beta1 = self._beta1

        if beta2 is None:
            beta2 = self._beta2

        if beta3 is None:
            beta3 = self._beta3

        if tau is None:
            tau = self._tau

        theta = t / tau
        expTerm = np.exp(-theta)
        yld = beta1
        yld += beta2 * (1.0 - expTerm) / theta
        yld += beta3 * ((1.0 - expTerm) / theta - expTerm)
        return yld

###############################################################################


class FinCurveFitMethodNelsonSiegelSvensson():

    def __init__(self, tau1=None, tau2=None):
        self._parentType = FinCurveFitMethod
        self._beta1 = None
        self._beta2 = None
        self._beta3 = None
        self._beta4 = None
        self._tau1 = tau1
        self._tau2 = tau2

    def interpolatedYield(self, t,
                          beta1=None, beta2=None, beta3=None, beta4=None,
                          tau1=None, tau2=None):

        # Careful if we get a time zero point
        t = np.maximum(t, 1e-10)

        if beta1 is None:
            beta1 = self._beta1

        if beta2 is None:
            beta2 = self._beta2

        if beta3 is None:
            beta3 = self._beta3

        if beta4 is None:
            beta4 = self._beta4

        if tau1 is None:
            tau1 = self._tau1

        if tau2 is None:
            tau2 = self._tau2

        theta1 = t / tau1
        theta2 = t / tau2
        expTerm1 = np.exp(-theta1)
        expTerm2 = np.exp(-theta2)
        yld = beta1
        yld += beta2 * (1.0 - expTerm1) / theta1
        yld += beta3 * ((1.0 - expTerm1) / theta1 - expTerm1)
        yld += beta4 * ((1.0 - expTerm2) / theta2 - expTerm2)
        return yld

###############################################################################

class FinCurveFitMethodBSpline():

    def __init__(self, power=3, knots=[1,3,5,10]):
        self._parentType = FinCurveFitMethod
        self._power = power
        self._knots = knots

    def interpolatedYield(self, t):
        t = np.maximum(t, 1e-10)
        yld = splev(t, self._spline)
        return yld
