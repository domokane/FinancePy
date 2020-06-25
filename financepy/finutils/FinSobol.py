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
from functools import partial

from numba import njit, objmode

sArr = np.empty(0, dtype='int64')
aArr = np.empty(0, dtype='int64')
m_i = np.empty((0, 0), dtype='int64')
dirname = os.path.abspath(os.path.dirname(__file__))
path = os.path.join(dirname, "sobolcoeff.csv")
del dirname

def bytesToIntArray(l, s):
    arr = s.split(b' ')
    npArr = np.array(arr, dtype='int64')
    result = np.pad(npArr, (0, l - len(npArr)))
    return result

def loadSobolCoeff(dimension):
    global sArr
    global aArr
    global m_i
    global path

    rowsLoaded = len(sArr)
    if (dimension - 1 > rowsLoaded):
        sArrNew, aArrNew = np.loadtxt(path, delimiter=',', dtype='int64', skiprows=rowsLoaded+1, max_rows=dimension-1, usecols=(0,1), unpack=True)
        maxLen = sArrNew[-1]
        m_iNew = np.loadtxt(path, delimiter=',', dtype='int64', skiprows=rowsLoaded+1, max_rows=dimension-1, usecols=(2,), converters={2: partial(bytesToIntArray, maxLen)})
        
        sArr = np.concatenate((sArr, sArrNew))
        aArr = np.concatenate((aArr, aArrNew))
        padAmount = m_iNew.shape[1] - m_i.shape[1]
        m_i = np.pad(m_i, ((0,0), (0,padAmount)))
        m_i = np.concatenate((m_i, m_iNew))
    return sArr, aArr, m_i


@njit
def generateSobol(numPoints, dimension):
    """
    Sobol points generator based on graycode order.
    This function is translated from the original c++ program.
    Original c++ program: https://web.maths.unsw.edu.au/~fkuo/sobol/
    Args:
         numPoints (int): number of points (cannot be greater than 2^32)
         dimension (int): number of dimensions
     Return:
         point (nparray): 2-dimensional array with row as the point and column as the dimension.
    """

    with objmode(sArr='int64[:]', aArr='int64[:]', m_i='int64[:,:]'):
        sArr, aArr, m_i = loadSobolCoeff(dimension)

    # ll = number of bits needed
    ll = int(np.ceil(np.log(numPoints)/np.log(2.0)))

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
    x = np.zeros(numPoints)
    x[0] = 0
    for i in range(1, numPoints):
        x[i] = int(x[i-1]) ^ int(v[c[i-1]])
        points[i, 0] = x[i]/(2**32)

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
                    v[i] = int(v[i]) ^ (((int(a) >> int(s-1-k)) & 1) * int(v[i-k]))

        # Evalulate X[0] to X[N-1], scaled by pow(2,32)
        x = np.zeros(numPoints)
        x[0] = 0
        for i in range(1, numPoints):
            x[i] = int(x[i-1]) ^ int(v[c[i-1]])
            points[i, j] = x[i]/(2**32)

    return points
