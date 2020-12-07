###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################


from math import exp, sqrt, fabs, log
from numba import njit, boolean, int64, float64, vectorize
import numpy as np
from .FinError import FinError

PI = 3.14159265358979323846
INVROOT2PI = 0.3989422804014327

ONE_MILLION = 1000000
TEN_MILLION = 10000000
ONE_BILLION = 1000000000

###############################################################################
# TODO: Move this somewhere else.
###############################################################################


@njit(fastmath=True, cache=True)
def accruedInterpolator(tset: float,  # Settlement time in years
                        couponTimes: np.ndarray,
                        couponAmounts: np.ndarray):
    ''' Fast calulation of accrued interest using an Actual/Actual type of
    convention. This does not calculate according to other conventions. '''

    numCoupons = len(couponTimes)

    for i in range(1, numCoupons):

        pct = couponTimes[i - 1]
        nct = couponTimes[i]
        denom = (nct-pct)
       
        if tset >= pct and tset < nct:
            accdFrac = (tset-pct) / denom
            accdCpn = accdFrac * couponAmounts[i]
            return accdCpn

    # TODO: NEED TO REVISIT THIS TODO
    return 0.0
    print("t", tset)
    print("CPN TIMES", couponTimes)
    print("CPN AMNTS", couponAmounts)
 
    raise FinError("Failed to calculate accrued")

###############################################################################


@njit(boolean(int64), fastmath=True, cache=True)
def isLeapYear(y: int):
    ''' Test whether year y is a leap year - if so return True, else False '''
    leapYear = ((y % 4 == 0) and (y % 100 != 0) or (y % 400 == 0))
    return leapYear

###############################################################################


@njit(float64[:](float64[:], float64), fastmath=True, cache=True)
def scale(x: np.ndarray,
          factor: float):
    ''' Scale all of the elements of an array by the same amount factor. '''
    x_scale = np.empty(len(x))
    for i in range(0, len(x)):
        x_scale[i] = x[i] * factor
    return x_scale

###############################################################################


@njit(boolean(float64[:]), fastmath=True, cache=True)
def testMonotonicity(x: np.ndarray):
    ''' Check that an array of doubles is monotonic and strictly increasing.'''
    for i in range(1, len(x)):
        if x[i] <= x[i-1]:
            return False
    return True

###############################################################################


@njit(fastmath=True, cache=True)
def testRange(x: np.ndarray,
              lower: float,
              upper: float):
    ''' Check that all of the values of an array fall between a lower and
    upper bound. '''
    for i in range(0, len(x)):
        if x[i] < lower:
            raise FinError("Value below lower.")
        if x[i] > upper:
            raise FinError("Value above upper.")

###############################################################################


@njit(fastmath=True, cache=True)
def maximum(a: np.ndarray,
            b: np.ndarray):
    ''' Determine the array in which each element is the maximum of the
    corresponding element in two equally length arrays a and b. '''

    n = len(a)
    out = [0.0] * n

    for i in range(0, n):
        if a[i] > b[i]:
            out[i] = a[i]
        else:
            out[i] = b[i]

    return out

###############################################################################


@njit(float64[:](float64[:, :]), fastmath=True, cache=True)
def maxaxis(s: np.ndarray):
    ''' Perform a search for the vector of maximum values over an axis of a
    2D Numpy Array '''

    shp = s.shape

    maxVector = np.empty(shp[0])

    for i in range(0, shp[0]):
        xmax = s[i, 0]
        for j in range(1, shp[1]):
            x = s[i, j]
            if x > xmax:
                xmax = x

        maxVector[i] = xmax

    return maxVector

###############################################################################


@njit(float64[:](float64[:, :]), fastmath=True, cache=True)
def minaxis(s: np.ndarray):
    ''' Perform a search for the vector of minimum values over an axis of a
    2D Numpy Array '''
    shp = s.shape

    minVector = np.empty(shp[0])

    for i in range(0, shp[0]):
        xmin = s[i, 0]
        for j in range(1, shp[1]):
            x = s[i, j]
            if x < xmin:
                xmin = x

        minVector[i] = xmin

    return minVector

###############################################################################


@njit(fastmath=True, cache=True)
def covar(a: np.ndarray,
          b: np.ndarray):
    ''' Calculate the Covariance of two arrays of numbers.
    TODO: check that this works well for Numpy Arrays and add NUMBA function
    signature to code. Do test of timings against Numpy. '''

    n = len(a)
    ma = 0.0
    mb = 0.0
    mab = 0.0
    ma2 = 0.0
    mb2 = 0.0

    for i in range(0, n):
        ma = ma + a[i]
        mb = mb + b[i]
        ma2 = ma2 + a[i]**2
        mb2 = mb2 + b[i]**2
        mab = mab + a[i] * b[i]

    ma /= n
    mb /= n
    ma2 /= n
    mb2 /= n
    mab /= n

    caa = ma2 - ma * ma
    cab = mab - ma * mb
    cbb = mb2 - mb * mb

    m = [[0.0, 0.0], [0.0, 0.0]]
    m[0][0] = caa
    m[1][0] = cab
    m[0][1] = cab
    m[1][1] = cbb
    return m

###############################################################################


@njit(float64(float64, float64), fastmath=True, cache=True)
def pairGCD(v1: float,
            v2: float):
    ''' Determine the Greatest Common Divisor of two integers using Euclid's
    algorithm. TODO - compare this with math.gcd(a,b) for speed. Also examine
    to see if I should not be declaring inputs as integers for NUMBA. '''

    if v1 == 0 or v2 == 0:
        return 0

    while v2 != 0:
        temp = v2
        factor = v1 / v2
        v2 = v1 - factor * v2
        v1 = temp

    pairGCD = abs(v1)
    return pairGCD

###############################################################################


@njit(fastmath=True, cache=True)
def nprime(x: float):
    '''Calculate the first derivative of the Cumulative Normal CDF which is
    simply the PDF of the Normal Distribution '''

    InvRoot2Pi = 0.3989422804014327
    return np.exp(-x * x / 2.0) * InvRoot2Pi

###############################################################################


@njit(fastmath=True, cache=True)
def heaviside(x: float):
    ''' Calculate the Heaviside function for x '''
    if x >= 0.0:
        return 1.0
    return 0.0

###############################################################################


@njit(fastmath=True, cache=True)
def frange(start: int,
           stop: int,
           step: int):
    ''' Calculate a range of values from start in steps of size step. Ends as
    soon as the value equals or exceeds stop. '''
    x = []
    while start <= stop:
        x.append(start)
        start += step

    return x

###############################################################################


@njit(fastmath=True, cache=True)
def normpdf(x: float):
    ''' Calculate the probability density function for a Gaussian (Normal)
    function at value x'''
    InvRoot2Pi = 0.3989422804014327
    return np.exp(-x * x / 2.0) * InvRoot2Pi

###############################################################################

@njit(float64(float64), fastmath=True, cache=True)
def N(x):
    ''' Fast Normal CDF function based on Hull OFAODS  4th Edition Page 252. 
    This function is accurate to 6 decimal places. '''

    a1 = 0.319381530
    a2 = -0.356563782
    a3 = 1.781477937
    a4 = -1.821255978
    a5 = 1.330274429
    g = 0.2316419

    k = 1.0 / (1.0 + g * fabs(x))
    k2 = k * k
    k3 = k2 * k
    k4 = k3 * k
    k5 = k4 * k

    if x >= 0.0:
        c = (a1 * k + a2 * k2 + a3 * k3 + a4 * k4 + a5 * k5)
        phi = 1.0 - c * exp(-x*x/2.0) * INVROOT2PI
    else:
        phi = 1.0 - N(-x)

    return phi

###############################################################################

@vectorize([float64(float64)], fastmath=True, cache=True)
def NVect(x):
    return N(x)

###############################################################################

@vectorize([float64(float64)], fastmath=True, cache=True)
def NPrimeVect(x):
    return nprime(x)

###############################################################################

@njit(float64(float64), fastmath=True, cache=True)
def normcdf_integrate(x: float):
    ''' Calculation of Normal Distribution CDF by simple integration
    which can become exact in the limit of the number of steps tending
    towards infinity. This function is used for checking as it is slow
    since the number of integration steps is currently hardcoded to 10,000.'''
    lower = -6.0
    upper = x
    numSteps = 10000
    dx = (upper - lower) / numSteps
    InvRoot2Pi = 0.3989422804014327

    x = lower
    fx = exp(-x * x / 2.0)
    integral = fx / 2.0

    for _ in range(0, numSteps - 1):
        x = x + dx
        fx = exp(-x * x / 2.0)
        integral += fx

    x = x + dx
    fx = exp(-x * x / 2.0)
    integral += fx / 2.0
    integral *= InvRoot2Pi * dx
    return integral

###############################################################################


@njit(float64(float64), fastmath=True, cache=True)
def normcdf_slow(z: float):
    '''  Calculation of Normal Distribution CDF accurate to 1d-15. This
    method is faster than integration but slower than other approximations.
    Reference: J.L. Schonfelder, Math Comp 32(1978), pp 1232-1240. '''

    a = [0.0] * 25
    bp = 0.0

    RTWO = 1.4142135623731

    a[0] = 0.6101430819232
    a[1] = -0.434841272712578
    a[2] = 0.176351193643605
    a[3] = -6.07107956092494E-02
    a[4] = 1.77120689956941E-02
    a[5] = -4.32111938556729E-03
    a[6] = 8.54216676887099E-04
    a[7] = -1.27155090609163E-04
    a[8] = 1.12481672436712E-05
    a[9] = 3.13063885421821E-07
    a[10] = -2.70988068537762E-07
    a[11] = 3.07376227014077E-08
    a[12] = 2.51562038481762E-09
    a[13] = -1.02892992132032E-09
    a[14] = 2.99440521199499E-11
    a[15] = 2.60517896872669E-11
    a[16] = -2.63483992417197E-12
    a[17] = -6.43404509890636E-13
    a[18] = 1.12457401801663E-13
    a[19] = 1.72815333899861E-14
    a[20] = -4.26410169494238E-15
    a[21] = -5.45371977880191E-16
    a[22] = 1.58697607761671E-16
    a[23] = 2.0899837844334E-17
    a[24] = -5.900526869409E-18

    xa = abs(z) / RTWO

    if xa > 100:
        p = 0
    else:
        t = (8 * xa - 30) / (4 * xa + 15)
        bm = 0.0
        b = 0.0

        for i in range(0, 25):
            bp = b
            b = bm
            bm = t * b - bp + a[24 - i]

        p = exp(-xa * xa) * (bm - bp) / 4

    if z > 0:
        p = 1.0 - p

    return p

###############################################################################


# @njit([float64(float64)], fastmath=True, cache=True)
# def normcdf(x: float):
#     ''' This is the Normal CDF function which forks to one of three of the
#     implemented approximations. This is based on the choice of the fast flag
#     variable. A value of 1 is the fast routine, 2 is the slow and 3 is the
#     even slower integration scheme. '''

#     return normcdf_fast(x)

#     # if fastFlag == 1:
#     #     return normcdf_fast(x)
#     # elif fastFlag == 2:
#     #     return normcdf_slow(x)
#     # elif fastFlag == 3:
#     #     return normcdf_integrate(x)
#     # else:
#     #     return 0.0

###############################################################################


# @vectorize([float64(float64)], fastmath=True, cache=True)
# def N(x: float):
#     ''' This is the shortcut to the default Normal CDF function and currently
#     is hardcoded to the fastest of the implemented routines. This is the most
#     widely used way to access the Normal CDF. '''
#     return normcdf_fast(x)

###############################################################################


@njit(fastmath=True, cache=True)
def phi3(b1: float,
         b2: float,
         b3: float,
         r12: float,
         r13: float,
         r23: float):
    ''' Bivariate Normal CDF function to upper limits $b1$ and $b2$ which uses
    integration to perform the innermost integral. This may need further
    refinement to ensure it is optimal as the current range of integration is
    from -7 and the integration steps are dx = 0.001. This may be excessive.'''

    dx = 0.001
    lowerLimit = -7
    upperlimit = b1
    numPoints = int((b1 - lowerLimit) / dx)
    dx = (upperlimit - lowerLimit) / numPoints
    x = lowerLimit

    r12p = sqrt(1.0 - r12 * r12)
    r13p = sqrt(1.0 - r13 * r13)
    r123 = (r23 - r12 * r13) / r12p / r13p

    v = 0.0

    for _ in range(1, numPoints + 1):
        dp = N(x + dx) - N(x)
        h = (b2 - r12 * x) / r12p
        k = (b3 - r13 * x) / r13p
        bivariate = M(h, k, r123)
        v = v + bivariate * dp
        x += dx

    return v

###############################################################################


@njit(fastmath=True, cache=True)
def norminvcdf(p):
    '''  This algorithm computes the inverse Normal CDF and is based on the
    algorithm found at (http:#home.online.no/~pjacklam/notes/invnorm/)
    which is by John Herrero (3-Jan-03) '''

    # Define coefficients in rational approximations
    a1 = -39.6968302866538
    a2 = 220.946098424521
    a3 = -275.928510446969
    a4 = 138.357751867269
    a5 = -30.6647980661472
    a6 = 2.50662827745924

    b1 = -54.4760987982241
    b2 = 161.585836858041
    b3 = -155.698979859887
    b4 = 66.8013118877197
    b5 = -13.2806815528857

    c1 = -7.78489400243029E-03
    c2 = -0.322396458041136
    c3 = -2.40075827716184
    c4 = -2.54973253934373
    c5 = 4.37466414146497
    c6 = 2.93816398269878

    d1 = 7.78469570904146E-03
    d2 = 0.32246712907004
    d3 = 2.445134137143
    d4 = 3.75440866190742

    inverseCDF = 0.0

    # Define break-points
    p_low = 0.02425
    p_high = 1.0 - p_low

    # If argument out of bounds, raise error
    if p < 0.0 or p > 1.0:
        raise FinError("p must be between 0.0 and 1.0")

    if p == 0.0:
        p = 1e-10

    if p == 1.0:
        p = 1.0 - 1e-10

    if p < p_low:
        # Rational approximation for lower region
        q = sqrt(-2.0 * log(p))
        inverseCDF = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5)
                    * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0)
    elif p <= p_high:
        # Rational approximation for lower region
        q = p - 0.5
        r = q * q
        inverseCDF = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * \
            q / (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1.0)
    elif p < 1.0:
        # Rational approximation for upper region
        q = sqrt(-2.0 * log(1 - p))
        inverseCDF = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5)
                       * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0)

    return inverseCDF

###############################################################################
# This is used for consistency with Haug and its conciseness. Consider renaming
# phi2 to M
@njit(fastmath=True, cache=True)
def M(a, b, c):
    return phi2(a, b, c)

###############################################################################


@njit(float64(float64, float64, float64), fastmath=True, cache=True)
def phi2(h1, hk, r):
    ''' Drezner and Wesolowsky implementation of bi-variate normal '''

#    if abs(r) > 0.9999999:
#        raise FinError("Phi2: |Correlation| > 1")

    x = [0.0, 0.0, 0.0, 0.0, 0.0]
    w = [0.0, 0.0, 0.0, 0.0, 0.0]

    x[0] = 0.04691008
    x[1] = 0.23076534
    x[2] = 0.5
    x[3] = 0.76923466
    x[4] = 0.95308992

    w[0] = 0.018854042
    w[1] = 0.038088059
    w[2] = 0.0452707394
    w[3] = 0.038088059
    w[4] = 0.018854042

    h2 = hk
    h12 = (h1 * h1 + h2 * h2) * 0.5
    bv = 0.0

    if fabs(r) < 0.7 or fabs(h1) > 35 or fabs(h2) > 35:

        h3 = h1 * h2

        for i in range(0, 5):
            r1 = r * x[i]
            rr2 = 1.0 - r1 * r1
            bv = bv + w[i] * exp((r1 * h3 - h12) / rr2) / sqrt(rr2)

        bv = N(h1) * N(h2) + r * bv
    else:
        r2 = 1.0 - r * r
        r3 = sqrt(r2)

        if r < 0.0:
            h2 = -h2

        h3 = h1 * h2
        h7 = exp(-h3 * 0.5)

        if r2 != 0.0:
            h6 = abs(h1 - h2)
            h5 = h6 * h6 * 0.5
            h6 = h6 / r3
            aa = 0.5 - h3 * 0.125
            ab = 3.0 - 2.0 * aa * h5
            bv = 0.13298076 * h6 * ab * \
                N(-h6) - exp(-h5 / r2) * (ab + aa * r2) * 0.053051647

            for i in range(0, 5):
                r1 = r3 * x[i]
                rr = r1 * r1
                r2 = sqrt(1.0 - rr)
                bv = bv - w[i] * exp(-h5 / rr) * \
                    (exp(-h3 / (1.0 + r2)) / r2 / h7 - 1.0 - aa * rr)

        if r > 0.0:
            bv = bv * r3 * h7 + N(min(h1, h2))
        else:
            if h1 < h2:
                bv = -bv * r3 * h7
            else:
                bv = -bv * r3 * h7 + N(h1) + N(hk) - 1.0

    return bv

###############################################################################


@njit(float64[:, :](float64[:, :]), cache=True, fastmath=True)
def cholesky(rho):
    ''' Numba-compliant wrapper around Numpy cholesky function. '''
    chol = np.linalg.cholesky(rho)
    return chol

###############################################################################


@njit(fastmath=True, cache=True)
def corrMatrixGenerator(rho, n):
    ''' Utility function to generate a full rank n x n correlation matrix with
    a flat correlation structure and value rho. '''

    corrMatrix = np.zeros(shape=(n, n))
    for i in range(0, n):
        corrMatrix[i, i] = 1.0
        for j in range(0, i):
            corrMatrix[i, j] = rho
            corrMatrix[j, i] = rho

    return corrMatrix

###############################################################################
