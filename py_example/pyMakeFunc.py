"""
Takes a Python module with defined vector field (and auxiliary functions),
computes the Jacobian symbolically, then writes a C++ source file which is
used by the ALGLIB minimizers to find local minimizers of the action.

The following elements are required in the equations Python module:
  problem_name - name of the problem (string_
  NX - number of state variables (int)
  NS - number of stimuli (int)
  syms - list of symbols, one for each state variable and parameter
  stims - list of symbols, one for each stimulus
"""

import importlib
import os, sys
import sympy as sym
import argparse as agp

# Open the Python module containing the ODE's.
# First, read in the name of the module file.
parser = agp.ArgumentParser()
parser.add_argument('-e', '--equations', action='store', dest='equations', default='equations.py')
args = parser.parse_args()
equations = args.equations

# Import the equations module
cwd = os.getcwd().rstrip('/')
module_loc = '{0}/{1}'.format(cwd, equations)
directory, module_name = os.path.split(module_loc)
module_name = os.path.splitext(module_name)[0]
path = list(sys.path)
sys.path.insert(0, directory)
try:
  eq = __import__(module_name)
finally:
  sys.path[:] = path
#import equations

# Get the problem name
try:
  Problem = eq.problem_name
except:
  print('ERROR: No problem name specified in equations file.')
  exit(1)

# Get number of variables & stimuli
try:
  NX = len(eq.syms)
except NameError:
  print('ERROR: Equations file must specify syms.')
  exit(1)
except:
  print('ERROR: Uknown error with syms.')
  exit(1)

try:
  NS = len(eq.stims)
except NameError:
  print('ERROR: Equations file must specify stims. Make it an empty list () if '
        + 'you don\'t have any stimuli.')
  exit(1)
except:
  print('ERROR: Uknown error with stims.')
  exit(1)

# Compute the Jacobian of the vector field
Jac = eq.VF.jacobian(eq.syms)

# Convert the vector field and Jacobian into C-code format
outVF = []
for i in range(NX):
  tmp = '%s' % sym.ccode(eq.VF[i])
  for j in range(NX):
    tmp = tmp.replace('%s'%eq.syms[j], 'x[%d]'%j)
  for j in range(NS):
    tmp = tmp.replace('%s'%eq.stims[j], 'sti[it][%d]'%j)
  outVF.append(tmp)

outDF=[]
for i in range(NX):
  for j in range(NX):
    tmp = '%s' % sym.ccode(Jac[i,j])
    for k in range(NX):
      tmp=tmp.replace('%s'%eq.syms[k], 'x[%d]'%k)
    for k in range(NS):
      tmp=tmp.replace('%s'%eq.stims[k], 'sti[it][%d]'%k)
    outDF.append(tmp)

################################################################################
# This is where we start writing the C++ source file itself.
f = open('func.cpp','w')
f.write('#define NX %d\t// dim of state variable + number of parameters\n\
#define NS %d\t// number of stimulus\n\
using namespace alglib;\n\
void func_origin(real_1d_array &x, int it, real_2d_array &sti, real_1d_array &dxdt);\n\
void func_DF(real_1d_array &x, int it, real_2d_array &sti, real_2d_array &Jac);\n'% (NX,NS))

f.write('\n\nvoid func_origin(real_1d_array &x, int it, real_2d_array &sti, real_1d_array &dxdt){\n\t//%s vector field\n'% Problem)
for i in range(NX):
  f.write('\tdxdt[%d]=%s;\n' %(i,outVF[i]))
f.write('\n}\n\n\
void func_DF(real_1d_array &x, int it, real_2d_array &sti, real_2d_array &Jac){\n\t//%s Jacobian matrix\n' % Problem )
for i in range(NX):
  for j in range(NX):
    f.write('\tJac[%d][%d]=%s;\n' %(i,j,outDF[i*NX+j]))
f.write('\n}\n')
print 'Generated func.cpp sucessfully\n'
f.close()
