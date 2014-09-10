"""
Defines sympy functions that return Hodgkin-Huxley network equations
(connected by gap junctions).
"""

import sympy

# "helper" functions
def rate(V, Vt, Vs, A, B, C):
    u = (V - Vt)/Vs
    den = np.exp(u) - B
    testA = np.abs(den) < 10.0**(-15.0)  # test for small denominators
    testB = B == 1
    testC = C == 0
    testind = np.logical_and(testA,testB,testC)
    not_testind = np.array(-1*(testind-1), dtype='bool')
    
    ratio = np.zeros(V.shape[0])
    ratio[not_testind] = (A[not_testind]*(V[not_testind]-Vt[not_testind]) \
                          + C[not_testind]) / den[not_testind]
    ratio[testind] = A[testind]*Vs[testind]
    return ratio

def rate(V, Vt, Vs, A, B, C):
