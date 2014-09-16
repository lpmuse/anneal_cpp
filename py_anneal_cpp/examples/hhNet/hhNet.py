"""
Defines sympy functions that return Hodgkin-Huxley network equations
(connected by gap junctions).
"""

import numpy as np
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
    return A / (Vs*den_poly)

def den(u):
    return sym.exp(u) - 1

def x_inf(V, Vta, Vsa, Aa, Ba, Ca, Vtb, Vsb, Ab, Bb, Cb):
    alpha = rate(V, Vta, Vsa, Aa, Ba, Ca)
    beta = rate(V, Vtb, Vsb, Ab, Bb, Cb)
    return alpha / (alpha + beta)

################################################################################
# set up some information about the network
NN = 4  # number of neurons

# the problem name
problem_name = 'hhNet'

# declare a dict to hold ALL symbols for this model
Sd = {}

# define symbols for stimuli, state variables, and parameters (to be estimated)
stims = []
for i in range(NN):
    sstr = 'Iinj_{0}'.format(i)
    Sd[sstr] = sym.symbols(sstr)
    stims.append(Sd[sstr])
NS = len(stims)

syms = []
for i in range(NN):
    Sd['VV_{0}'.format(i)] = sym.symbols('VV_{0}'.format(i))
    syms.append(Sd['VV_%d'%i])
    Sd['mm_{0}'.format(i)] = sym.symbols('mm_{0}'.format(i))
    syms.append(Sd['mm_%d'%i])
    Sd['hh_{0}'.format(i)] = sym.symbols('hh_{0}'.format(i))
    syms.append(Sd['hh_%d'%i])
    Sd['nn_{0}'.format(i)] = sym.symbols('nn_{0}'.format(i))
    syms.append(Sd['nn_%d'%i])
NX = len(syms)

# define the TOTAL number of POSSIBLE fixed parameters
NP = 37
# define all of the fixed parameters
NP_FIXED = 37
# load in parameter values
p = np.loadtxt('parameters.dat')
for i in range(NN):
    # injected current surface area factor
    Sd['pinj_{0}'.format(i)] = p[0+4*i]
    # maximal conductances
    Sd['gNa_{0}'.format(i)] = p[1+4*i]
    Sd['gK_{0}'.format(i)] = p[2+4*i]
    Sd['gL_{0}'.format(i)] = p[3+4*i]
    # reversal potentials
    Sd['ENa_{0}'.format(i)] = p[4+4*i]
    Sd['EK_{0}'.format(i)] = p[5+4*i]
    Sd['EL_{0}'.format(i)] = p[6+4*i]
    # alpha rate parameters, m
    Sd['Vtam_{0}'.format(i)] = p[7+4*i]
    Sd['Vsam_{0}'.format(i)] = p[8+4*i]
    Sd['Aam_{0}'.format(i)] = p[9+4*i]
    Sd['Bam_{0}'.format(i)] = p[10+4*i]
    Sd['Cam_{0}'.format(i)] = p[11+4*i]
    # alpha rate parameters, h
    Sd['Vtah_{0}'.format(i)] = p[12+4*i]
    Sd['Vsah_{0}'.format(i)] = p[13+4*i]
    Sd['Aah_{0}'.format(i)] = p[14+4*i]
    Sd['Bah_{0}'.format(i)] = p[15+4*i]
    Sd['Cah_{0}'.format(i)] = p[16+4*i]
    # alpha rate parameters, n
    Sd['Vtan_{0}'.format(i)] = p[17+4*i]
    Sd['Vsan_{0}'.format(i)] = p[18+4*i]
    Sd['Aan_{0}'.format(i)] = p[19+4*i]
    Sd['Ban_{0}'.format(i)] = p[20+4*i]
    Sd['Can_{0}'.format(i)] = p[21+4*i]
    # beta rate parameters, m
    Sd['Vtbm_{0}'.format(i)] = p[22+4*i]
    Sd['Vsbm_{0}'.format(i)] = p[23+4*i]
    Sd['Abm_{0}'.format(i)] = p[24+4*i]
    Sd['Bbm_{0}'.format(i)] = p[25+4*i]
    Sd['Cbm_{0}'.format(i)] = p[26+4*i]
    # beta rate parameters, h
    Sd['Vtbh_{0}'.format(i)] = p[27+4*i]
    Sd['Vsbh_{0}'.format(i)] = p[28+4*i]
    Sd['Abh_{0}'.format(i)] = p[29+4*i]
    Sd['Bbh_{0}'.format(i)] = p[30+4*i]
    Sd['Cbh_{0}'.format(i)] = p[31+4*i]
    # beta rate parameters, n
    Sd['Vtbn_{0}'.format(i)] = p[32+4*i]
    Sd['Vsbn_{0}'.format(i)] = p[33+4*i]
    Sd['Abn_{0}'.format(i)] = p[34+4*i]
    Sd['Bbn_{0}'.format(i)] = p[35+4*i]
    Sd['Cbn_{0}'.format(i)] = p[36+4*i]

# load in coupling strengths
G = np.loadtxt('G.dat')
if NN == 1:
    G = sym.Matrix(np.array([G]))
else:
    G = sym.Matrix(G)

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
for i in range(NN):
    Sd['UU_{0}'.format(i)] = Sd['VV_{0}'.format(i)] * Vdiff + Vmin

dxdt_list = []
# make a list of all the voltages for interneuron coupling
allU = []
for i in range(NN):
    allU.append(Sd['UU_%d'%i])
allU = sym.Matrix(1, NN, allU)
# vector of 1's for part of coupling term
ONES = sym.ones(1, NN)
# sum over cols of G
if NN == 1:
    Gsum = sym.Matrix([G[0] * 1.0])
else:
    Gsum = G.dot(ONES)
# take dot product with voltages
if NN == 1:
    Gdot = sym.Matrix([G.dot(allU)])
else:
    Gdot = G.dot(allU)

for i in range(NN):
    dxdt_list.append((Sd['gNa_%d'%i]*pow(Sd['mm_%d'%i],3) \
        * Sd['hh_%d'%i]*(Sd['ENa_%d'%i] - Sd['UU_%d'%i]) \
        + Sd['gK_%d'%i]*pow(Sd['nn_%d'%i],4)*(Sd['EK_%d'%i] - Sd['UU_%d'%i]) \
        + Sd['gL_%d'%i]*(Sd['EL_%d'%i] - Sd['UU_%d'%i]) \
        + Sd['pinj_%d'%i]*Sd['Iinj_%d'%i] \
        + Gdot[i] - Gsum[i]*Sd['UU_%d'%i]) / Vdiff)
    GV_list = [Sd['mm_%d'%i], Sd['hh_%d'%i], Sd['nn_%d'%i]]
    Vkin_list = [Sd['Vtam_%d'%i], Sd['Vsam_%d'%i], Sd['Vtbm_%d'%i], Sd['Vsbm_%d'%i],
                 Sd['Vtah_%d'%i], Sd['Vsah_%d'%i], Sd['Vtbh_%d'%i], Sd['Vsbh_%d'%i],
                 Sd['Vtan_%d'%i], Sd['Vsan_%d'%i], Sd['Vtbn_%d'%i], Sd['Vsbn_%d'%i]]
    GV_coef_list = [Sd['Aam_%d'%i], Sd['Bam_%d'%i], Sd['Cam_%d'%i],
                    Sd['Abm_%d'%i], Sd['Bbm_%d'%i], Sd['Cbm_%d'%i],
                    Sd['Aah_%d'%i], Sd['Bah_%d'%i], Sd['Cah_%d'%i],
                    Sd['Abh_%d'%i], Sd['Bbh_%d'%i], Sd['Cbh_%d'%i],
                    Sd['Aan_%d'%i], Sd['Ban_%d'%i], Sd['Can_%d'%i],
                    Sd['Abn_%d'%i], Sd['Bbn_%d'%i], Sd['Cbn_%d'%i]]
    
    for j,gv in enumerate(GV_list):
        Sd['Aa_%d'%i], Sd['Ba_%d'%i], Sd['Ca_%d'%i], Sd['Ab_%d'%i], \
                     Sd['Bb_%d'%i], Sd['Cb_%d'%i] = GV_coef_list[6*j:6*j+6]
        Sd['Vta_%d'%i], Sd['Vsa_%d'%i], Sd['Vtb_%d'%i], \
                     Sd['Vsb_%d'%i] = Vkin_list[4*j:4*j+4]
        if Sd['Aa_%d'%i] != 0 and Sd['Ba_%d'%i] == 1 and Sd['Ca_%d'%i] == 0:
            alpha = rate_exp(Sd['UU_%d'%i], Sd['Vta_%d'%i], Sd['Vsa_%d'%i], \
                             Sd['Aa_%d'%i], Sd['Ba_%d'%i], Sd['Ca_%d'%i])
        else:
            alpha = rate(Sd['UU_%d'%i], Sd['Vta_%d'%i], Sd['Vsa_%d'%i],
                         Sd['Aa_%d'%i], Sd['Ba_%d'%i], Sd['Ca_%d'%i])
        if Sd['Ab_%d'%i] != 0 and Sd['Bb_%d'%i] == 1 and Sd['Cb_%d'%i] == 0:
            beta = rate_exp(Sd['UU_%d'%i], Sd['Vtb_%d'%i], Sd['Vsb_%d'%i],
                            Sd['Ab_%d'%i], Sd['Bb_%d'%i], Sd['Cb_%d'%i])
        else:
            beta = rate(Sd['UU_%d'%i], Sd['Vtb_%d'%i], Sd['Vsb_%d'%i],
                        Sd['Ab_%d'%i], Sd['Bb_%d'%i], Sd['Cb_%d'%i])
        dxdt_list.append(alpha*(1 - gv) - beta*gv);
    
    for i in range(NP-NP_FIXED):
        dxdt_list.append(0);

VF = sym.Matrix(1, NX, dxdt_list)
