# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 09:26:27 2016

@author: Dominic O'Kane
"""

import numpy as np
import matplotlib.pyplot as plt

from ...finutils.FinDate import FinDate
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinMath import scale
from .FinCurveFitMethod import FinCurveFitMethodPolynomial
from .FinCurveFitMethod import FinCurveFitMethodNelsonSiegel
from .FinCurveFitMethod import FinCurveFitMethodNelsonSiegelSvensson
from .FinCurveFitMethod import FinCurveFitMethodBSpline

from scipy.optimize import curve_fit
from scipy.interpolate import splrep

##########################################################################
# TODO: Inherit from FinCurve and add df method
# TODO: Add fitting optimiser to take in bonds and do a best fit
##########################################################################


class FinFitBondYieldCurve():
    ''' Class to do fitting of the yield curve and to enable interpolation of 
    yields. Because yields assume a flat term structure for each bond, this 
    class does not allow discounting to be done and so does not inherit from 
    FinCurve. It should only be used for visualisation and simple 
    interpolation but not for full term-structure-consistent pricing. '''

    def __init__(self, settlementDate, bonds, ylds, curveFitMethod):
        ''' Fit the curve to a set of bond yields using the type of curve 
        specified. Weights can be provided to deal with illiquid bonds. '''

        self._settlementDate = settlementDate
        self._bonds = bonds
        self._ylds = np.array(ylds)
        self._curveFitMethod = curveFitMethod

        fitType = type(self._curveFitMethod)

        yearsToMaturities = []
        for bond in bonds:
            maturityYears = (bond._maturityDate-settlementDate)/gDaysInYear
            yearsToMaturities.append(maturityYears)
        self._yearsToMaturity = np.array(yearsToMaturities)

        if fitType is FinCurveFitMethodPolynomial:

            d = curveFitMethod._power
            coeffs = np.polyfit(self._yearsToMaturity, self._ylds, deg=d)
            self._curveFitMethod._coeffs = coeffs

        elif fitType is FinCurveFitMethodNelsonSiegel:

            xdata = self._yearsToMaturity
            ydata = self._ylds
#            bnds = [(0.0, -1.0, -1.0, 0.0), (1.0, 1.0, 1.0, 100.0)]
            popt, pcov = curve_fit(self._curveFitMethod.interpolatedYield,
                                   xdata, ydata)

            self._curveFitMethod._beta1 = popt[0]
            self._curveFitMethod._beta2 = popt[1]
            self._curveFitMethod._beta3 = popt[2]
            self._curveFitMethod._tau = popt[3]

        elif fitType is FinCurveFitMethodNelsonSiegelSvensson:

            xdata = self._yearsToMaturity
            ydata = self._ylds
            bnds = [(0., -1, -1, -1, 0, 1), (1, 1, 1, 1, 10, 100)]
            popt, pcov = curve_fit(self._curveFitMethod.interpolatedYield,
                                   xdata, ydata, bounds=bnds)

            self._curveFitMethod._beta1 = popt[0]
            self._curveFitMethod._beta2 = popt[1]
            self._curveFitMethod._beta3 = popt[2]
            self._curveFitMethod._beta4 = popt[3]
            self._curveFitMethod._tau1 = popt[4]
            self._curveFitMethod._tau2 = popt[5]

        elif fitType is FinCurveFitMethodBSpline:

            xdata = self._yearsToMaturity
            ydata = self._ylds

            ''' Cubic splines as k=3 '''
            knots = [1, 3, 5, 10]
            spline = splrep(xdata, ydata, k=3, t=knots)
            self._curveFitMethod._spline = spline

        else:
            raise ValueError("Unrecognised curve fit type.")

##########################################################################

    def interpolatedYield(self, maturityDate):

        if type(maturityDate) is FinDate:
            t = (maturityDate - self._settlementDate) / gDaysInYear
        elif type(maturityDate) is list:
            t = maturityDate
        elif type(maturityDate) is np.ndarray:
            t = maturityDate
        else:
            raise ValueError("Unknown date type.")

        fit = self._curveFitMethod

        if type(fit) == FinCurveFitMethodPolynomial:
            yld = np.polyval(fit._coeffs, t)
        elif type(fit) == FinCurveFitMethodNelsonSiegel:
            yld = fit.interpolatedYield(t,
                                        self._curveFitMethod._beta1,
                                        self._curveFitMethod._beta2,
                                        self._curveFitMethod._beta3,
                                        self._curveFitMethod._tau)

        elif type(fit) == FinCurveFitMethodNelsonSiegelSvensson:
            yld = fit.interpolatedYield(t,
                                        self._curveFitMethod._beta1,
                                        self._curveFitMethod._beta2,
                                        self._curveFitMethod._beta3,
                                        self._curveFitMethod._beta4,
                                        self._curveFitMethod._tau1,
                                        self._curveFitMethod._tau2)
        elif type(fit) == FinCurveFitMethodBSpline:
            yld = fit.interpolatedYield(t)

        return yld

##########################################################################

    def display(self, title):
        ''' Calculation of forward rates. This function can return a vector
        of forward rates given a vector of times. '''

        plt.figure(figsize=(12, 6))
        plt.title(title)
        bond_ylds_scaled = scale(self._ylds, 100.0)
        plt.plot(self._yearsToMaturity, bond_ylds_scaled, 'o')
        plt.xlabel('Time to Maturity (years)')
        plt.ylabel('Yield To Maturity (%)')

        tmax = np.max(self._yearsToMaturity)
        t = np.linspace(0.0, int(tmax+0.5), 100)

        yld = self.interpolatedYield(t)
        yld = scale(yld, 100.0)
        plt.plot(t, yld, label=str(self._curveFitMethod))
        plt.legend(loc='lower right')
        plt.ylim((min(yld)-0.3, max(yld)*1.1))
        plt.grid(True)
