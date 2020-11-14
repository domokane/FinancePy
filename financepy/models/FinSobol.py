"""
BSD 3-Clause License

Copyright (c) 2019, Ghifari Adam Faza
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import os
import numpy as np
from numba import njit

from financepy.finutils.FinMath import norminvcdf

###############################################################################
# This code loads sobol coefficients from binary numpy file and allocates
# contents to static global variables.
###############################################################################

dirname = os.path.abspath(os.path.dirname(__file__))
path = os.path.join(dirname, "sobolcoeff.npz")

with np.load(path, mmap_mode='r') as f:
    sArr = np.array(f['sa'][0])
    aArr = np.array(f['sa'][1])
    m_i = f['c']

###############################################################################


@njit(cache=True)
def getGaussianSobol(numPoints, dimension):
    ''' Sobol Gaussian quasi random points generator based on graycode order.
    The generated points follow a normal distribution. '''
    points = getUniformSobol(numPoints, dimension)

    for i in range(numPoints):
        for j in range(dimension):
            points[i, j] = norminvcdf(points[i, j])
    return points

###############################################################################


@njit(cache=True)
def getUniformSobol(numPoints, dimension):
    ''' Sobol uniform quasi random points generator based on graycode order. 
    This function returns a 2D Numpy array of values where the number of rows
    is the number of draws and the number of columns is the number of 
    dimensions of the random values. Each dimension has the same number of 
    random draws. Each column of random numbers is ordered so as not to
    correlate, i.e be independent from any other column.'''
    
    global sArr
    global aArr
    global m_i

    # ll = number of bits needed
    ll = int(np.ceil(np.log(numPoints+1)/np.log(2.0)))

    # c[i] = index from the right of the first zero bit of i
    c = np.zeros(numPoints, dtype=np.int64)
    c[0] = 1
    for i in range(1, numPoints):
        c[i] = 1
        value = i
        while value & 1:
            value >>= 1
            c[i] += 1

    # points initialization
    points = np.zeros((numPoints, dimension))

    # ----- Compute the first dimension -----
    # Compute direction numbers v[1] to v[L], scaled by 2**32
    v = np.zeros(ll+1)
    for i in range(1, ll+1):
        v[i] = 1 << (32-i)
        v[i] = int(v[i])

    #  Evalulate x[0] to x[N-1], scaled by 2**32
    x = np.zeros(numPoints+1)
    for i in range(1, numPoints+1):
        x[i] = int(x[i-1]) ^ int(v[c[i-1]])
        points[i-1, 0] = x[i]/(2**32)

    # ----- Compute the remaining dimensions -----
    for j in range(1, dimension):
        # read parameters from file
        s = sArr[j-1]
        a = aArr[j-1]
        mm = m_i[j-1]
        m = np.concatenate((np.zeros(1), mm))

        # Compute direction numbers V[1] to V[L], scaled by 2**32
        v = np.zeros(ll+1)
        if ll <= s:
            for i in range(1, ll+1):
                v[i] = int(m[i]) << (32-i)

        else:
            for i in range(1, s+1):
                v[i] = int(m[i]) << (32-i)

            for i in range(s+1, ll+1):
                v[i] = int(v[i-s]) ^ (int(v[i-s]) >> s)
                for k in range(1, s):
                    v[i] = int(v[i]) ^ (((int(a) >> int(s-1-k)) & 1)
                                        * int(v[i-k]))

        # Evalulate X[0] to X[N-1], scaled by pow(2,32)
        x = np.zeros(numPoints+1)
        for i in range(1, numPoints+1):
            x[i] = int(x[i-1]) ^ int(v[c[i-1]])
            points[i-1, j] = x[i]/(2**32)

    return points

###############################################################################
