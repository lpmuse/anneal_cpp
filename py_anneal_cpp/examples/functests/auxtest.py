"""
Testing how the C++ code generator handles "auxiliary" functions.
"""

import sympy as sym

# auxiliary functions
def func1(x):
    return x

def func2(x):
    return 5*x

def func3(x):
    if x < 0:
        return 5*x
    else:
        return 10*x

################################################################################

# problem name
problem_name = 'auxtest'

# symbols for stimuli, state variables, and parameters
stims = []

syms = []
X_list = ['x']
NX = len(X_list)
for X in X_list:
    symbol_ptr = sym.symbols(X)
    tmp = '{0} = symbol_ptr'.format(X)
    exec tmp
    syms.append(symbol_ptr)

# define vector field symbolically
dxdt_list = [func2(x)]
VF = sym.Matrix(1, NX, dxdt_list)
