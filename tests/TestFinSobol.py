from financepy.finutils.FinSobol import getSobol, getGaussianSobol

import time
from numba import jit
import numpy as np
def test_FinSobol():
    
    numPoints = 1000
    dimensions = 3
    
    points = getSobol(numPoints, dimensions)

    for d in range(dimensions):
        av = 0.0
        var = 0.0
        
        for point in points[:, d]:
            av += point
            var += point ** 2
            
        av /= numPoints
        var /= numPoints
        
        avError = abs(av - (1/2))
        varError = abs(var - (1/3))
        assert(avError < 0.002)
        assert(varError < 0.002)
      
    numRepeats = 100
    numDimensions = 10

    start = time.time()
    for _ in range(numRepeats):
        getSobol(1000, numDimensions)
    end = time.time()
    print("Average time taken", (end - start) / numRepeats)

    start = time.time()    
    for _ in range(numRepeats):
        getGaussianSobol(1000, numDimensions)
    end = time.time()
    print("Average time taken", (end - start) / numRepeats)
    
"""
from functools import partial
from numba import objmode

def bytesToIntArray(l, s):
    arr = s.split(' ')
    npArr = np.array(arr, dtype='int64')
    result = np.pad(npArr, (0, l - len(npArr)))
    return result

path = "C:\\Users\\ferga\\Repositories\\financepy\\financepy\\finutils\\sobolcoeff.csv"
start = time.time()
sa = np.loadtxt(path, delimiter=',', dtype='int64', max_rows=1000, skiprows=1, usecols=(0,1), unpack=True)
print(time.time() - start)
maxLen = sa[0][-1]
c = np.loadtxt(path, delimiter=',', dtype='int64', skiprows=1, max_rows=1000, usecols=(2,), converters={2: partial(bytesToIntArray, maxLen)}, encoding='latin1')

print(time.time() - start)
@jit(nopython=False)
def pythonFunction():
    return noPythonFunction()


path = ("C:\\Users\\ferga\\Repositories\\financepy\\financepy\\finutils\\saved.npz")
print(sa)
print(c[1])
np.savez(path, sa=sa, c=c)

start = time.time()
with np.load(path) as f:
    #print(f['sa'])
    pass
print(time.time() - start)


@jit(nopython=True)
def noPythonFunction():
    global numbers
    print(numbers)
    return 1
"""

@jit(cache=True, nopython=True)    
def cachedFunction():
    print(getSobol(2,2))
    return getSobol(2,2)


test_FinSobol()
cachedFunction()