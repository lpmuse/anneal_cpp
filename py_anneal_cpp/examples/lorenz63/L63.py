# Model equations file for Lorenz 63.  This is specifically meant to be fed 
# into minAzero, ALGLIB L-BFGS minimizer, etc.

import sympy as sym

# define the problem name
problem_name = 'lorenz63'

# no stimuli; define as empty tuple
stims = []

# define lists of symbols for state variables and stimuli
syms = []
#X_list = ['x', 'y', 'z', 'sigma', 'r', 'b']
X_list = ['x', 'y', 'z']
NX = len(X_list)
for X in X_list:
    symbol_ptr = sym.symbols(X)
    tmp = '{0} = symbol_ptr'.format(X)
    exec tmp
    syms.append(symbol_ptr)

# define fixed parameters
sigma = 10.0
r = 28.0
b = 8.0/3.0

# now define the vector field symbolically
# note that each variable is going to be rescaled btwn 0 and 1
xmin, xmax = (-20.0, 20.0)
ymin, ymax = (-25.0, 25.0)
zmin, zmax = (0.0, 50.0)

xdiff = xmax - xmin
ydiff = ymax - ymin
zdiff = zmax - zmin

X = x*xdiff + xmin
Y = y*ydiff + ymin
Z = z*zdiff + zmin

dxdt_list = [(sigma * (Y - X)) / xdiff,
             (X * (r - Z) - Y) / ydiff,
             (X * Y - b * Z) / zdiff]
VF = sym.Matrix(1, NX, dxdt_list)

