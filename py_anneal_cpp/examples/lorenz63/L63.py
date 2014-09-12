# Model equations file for Lorenz 63.  This is specifically meant to be fed
# into minAzero, ALGLIB L-BFGS minimizer, etc.

import sympy as sym

# define the problem name
problem_name = 'lorenz63'

# no stimuli; define as empty tuple
stims = []

# define lists of symbols for state variables and stimuli
syms = []
X_list = ['x', 'y', 'z', 'sigma', 'r', 'b']
NX = len(X_list)
for X in X_list:
    symbol_ptr = sym.symbols(X)
    tmp = '{0} = symbol_ptr'.format(X)
    exec tmp
    syms.append(symbol_ptr)

# now define the vector field symbolically
dxdt_list = [sigma * (y - x),
             x * (r - z) - y,
             x * y - b * z,
             0, 0, 0]
VF = sym.Matrix(1, NX, dxdt_list)

