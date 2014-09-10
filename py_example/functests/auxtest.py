"""
Testing how the C++ code generator handles "auxiliary" functions.
"""

import sympy as sym

# auxiliary functions
def func1(x):
    return x

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
dxdt_list = [func1(x)]
VF = sym.Matrix(1, NX, dxdt_list)
