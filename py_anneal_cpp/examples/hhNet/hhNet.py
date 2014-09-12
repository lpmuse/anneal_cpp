"""
Defines sympy functions that return Hodgkin-Huxley network equations
(connected by gap junctions).
"""

import sympy as sym
import sympy.mpmath as mpm

# "helper" functions
def rate(V, Vt, Vs, A, B, C):
    """
    Generic form for the alpha and beta gating variable rate functions.
    """
    return (A*(V - Vt) + C) / (sym.exp((V - Vt)/Vs) - B)

def rate_exp(V, Vt, Vs, A, B, C):
    """
    The rate function when C = 0 and B = 1 needs to be expanded for small 
    (V - Vt)/Vs.  Since the sympy --> C++ code conversion can't handle
    conditional expressions, the expansion is done to very high order and
    used for ALL x!
    """
    u = (V - Vt)/Vs
    coef = mpm.taylor(den, 0, 25)
    den_poly = sym.expand(mpm.polyval(coef[::-1], u) / u)
    return A / den_poly

def den(u):
    return sym.exp(u) - 1

def x_inf(V, Vta, Vsa, Aa, Ba, Ca, Vtb, Vsb, Ab, Bb, Cb):
    alpha = rate(V, Vta, Vsa, Aa, Ba, Ca)
    beta = rate(V, Vtb, Vsb, Ab, Bb, Cb)
    return alpha / (alpha + beta)

################################################################################
# set up some information about the network
NN = 1  # number of neurons

# the problem name
problem_name = 'hhNet'

# define symbols for stimuli, state variables, and parameters (to be estimated)
stims = []
S_list = []
for i in range(NN):
    S_list.append('Iinj_{0}'.format(i+1))
NS = len(S_list)
for S in S_list:
    symbol_ptr = sym.symbols(S)
    tmp = '{0} = symbol_ptr'.format(S)
    exec tmp
    stims.append(symbol_ptr)

syms = []
X_list = []
for i in range(NN):
    X_list.append('VV'.format(i+1))
    X_list.append('mm'.format(i+1))
    X_list.append('hh'.format(i+1))
    X_list.append('nn'.format(i+1))
NX = len(X_list)
for X in X_list:
    symbol_ptr = sym.symbols(X)
    tmp = '{0} = symbol_ptr'.format(X)
    exec tmp
    syms.append(symbol_ptr)

# define the TOTAL number of POSSIBLE fixed parameters
NP = 37
# define all of the fixed parameters
NP_FIXED = 37
# injected current surface area factor
pinj = 1.0
# maximal conductances
gNa = 120.0
gK = 36.0
gL = 0.3
# reversal potentials
ENa = 45.0
EK = -82.0
EL = -59.387
# alpha rate parameters, m
Vtam = -45.0
Vsam = -10.0
Aam = -0.1
Bam = 1.0
Cam = 0.0
# alpha rate parameters, h
Vtah = -70.0
Vsah = 20.0
Aah = 0.0
Bah = 0.0
Cah = 0.07
# alpha rate parameters, n
Vtan = -60.0
Vsan = -10.0
Aan = -0.01
Ban = 1.0
Can = 0.0
# beta rate parameters, m
Vtbm = -70.0
Vsbm = 18.0
Abm = 0.0
Bbm = 0.0
Cbm = 4.0
# beta rate parameters, h
Vtbh = -40.0
Vsbh = -10.0
Abh = 0.0
Bbh = -1.0
Cbh = 1.0
# beta rate parameters, n
Vtbn = -70.0
Vsbn = 80.0
Abn = 0.0
Bbn = 0.0
Cbn = 0.125

# define the vector field symbolically
# NOTE: V is rescaled to be between 0 and 1.  This means applying the transformations:
#             V --> (V - V_min) / (V_max - V_min)
#             dVdt --> dVdt / (V_max - V_min)
# This means applying the following replacements:
#             V replaced with V*(V_max - V_min) + V_min
#             dVdt replaced with dVdt*(V_max - V_min)
Vmax = 20.0  # mV
Vmin = -80.0  # mV
Vdiff = Vmax - Vmin
# define U as a placeholder for transformed V
UU = sym.symbols('UU')
UU = VV * Vdiff + Vmin

dxdt_list = []
dxdt_list.append((gNa*pow(mm,3)*hh*(ENa - UU) + gK*pow(nn,4)*(EK - UU) + gL*(EL - UU) \
                  + pinj*Iinj_1) / Vdiff)
GV_list = [mm, hh, nn]
Vkin_list = [Vtam, Vsam, Vtbm, Vsbm,
             Vtah, Vsah, Vtbh, Vsbh,
             Vtan, Vsan, Vtbn, Vsbn]
GV_coef_list = [Aam, Bam, Cam, Abm, Bbm, Cbm,
                Aah, Bah, Cah, Abh, Bbh, Cbh,
                Aan, Ban, Can, Abn, Bbn, Cbn]
for i,gv in enumerate(GV_list):
    Aa, Ba, Ca, Ab, Bb, Cb = GV_coef_list[6*i:6*i+6]
    print Aa, Ba, Ca, Ab, Bb, Cb
    Vta, Vsa, Vtb, Vsb = Vkin_list[4*i:4*i+4]
    if Aa != 0 and Ba == 1 and Ca == 0:
        alpha = rate_exp(UU, Vta, Vsa, Aa, Ba, Ca)
    else:
        alpha = rate(UU, Vta, Vsa, Aa, Ba, Ca)
    if Ab != 0 and Bb == 1 and Cb == 0:
        beta = rate_exp(UU, Vtb, Vsb, Ab, Bb, Cb)
    else:
        beta = rate(UU, Vtb, Vsb, Ab, Bb, Cb)
    
    dxdt_list.append(alpha*(1 - gv) - beta*gv)

for i in range(NP-NP_FIXED):
    dxdt_list.append(0)

VF = sym.Matrix(1, NX, dxdt_list)
