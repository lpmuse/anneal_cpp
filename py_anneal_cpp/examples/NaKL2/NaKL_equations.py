import sympy as sym

# the problem name
problem_name = 'NaKL'

# define symbols for stimuli, state variables, and parameters
stims = []
STIMS_list = ['Iinj']
NS = len(STIMS_list)
for STIM in STIMS_list:
    symbol_ptr = sym.symbols(STIM)
    tmp = '{0} = symbol_ptr'.format(STIM)
    exec tmp
    stims.append(symbol_ptr)

syms = []
X_list = ['VV', 'nn', 'mm', 'hh', 'pinj', 'gNa', 'ENa', 'gK', 'EK', 'gM', 'Erest',
          'anV1', 'anV2', 't0n', 'epsn', 'amV1', 'amV2', 't0m', 'epsm', 'ahV1',
          'ahV2', 't0h', 'epsh']
NX = len(X_list)
for X in X_list:
    symbol_ptr = sym.symbols(X)
    tmp = '{0} = symbol_ptr'.format(X)
    exec tmp
    syms.append(symbol_ptr)

# define the vector field symbolically
dxdt_list = [ pinj * Iinj + gNa *pow(mm ,3)* hh *(ENa - VV) \
              + gK *pow(nn ,4)*(EK - VV) + gM *( Erest - VV ),
              0.5 * (1 + sym.tanh(( VV - anV1 )/ anV2 ) - 2 * nn ) \
              / ( t0n + epsn * (1 - sym.tanh(( VV - anV1 )/ anV2 ) \
                                * sym.tanh( VV - anV1 )/ anV2 )),
              0.5 * (1 + sym.tanh(( VV - amV1 )/ amV2 ) - 2 * mm ) \
              / ( t0m + epsm * (1 - sym.tanh(( VV - amV1 )/ amV2 ) \
                                * sym.tanh( VV - amV1 )/ amV2 )),
              0.5 * (1 + sym.tanh(( VV - ahV1 )/ ahV2 ) - 2 * hh) \
              / ( t0h + epsh * (1 - sym.tanh(( VV - ahV1 )/ ahV2 ) \
                                * sym.tanh( VV - ahV1 )/ ahV2 )),
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
VF = sym.Matrix(1, NX, dxdt_list)
