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
parser.add_argument('-e', '--equations', action='store', dest='equations',
                    default='equations.py')
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
outVF = []  # vector field
for i in range(NX):
  Fstr = '{0}'.format(sym.ccode(eq.VF[i]))
  for j in range(NX):
    Fstr = Fstr.replace('{0}'.format(eq.syms[j]),
                       'x[{0}]'.format(j))
  for j in range(NS):
    Fstr = Fstr.replace('{0}'.format(eq.stims[j]),
                        'sti[it][{0}]'.format(j))
  outVF.append(Fstr)

outDF = []  # Jacobian of the vector field
for i in range(NX):
  for j in range(NX):
    Fstr = '{0}'.format(sym.ccode(Jac[i,j]))
    for k in range(NX):
      Fstr = Fstr.replace('{0}'.format(eq.syms[k]),
                          'x[{0}]'.format(k))
    for k in range(NS):
      Fstr = Fstr.replace('{0}'.format(eq.stims[k]),
                          'sti[it][{0}]'.format(k))
    outDF.append(Fstr)

################################################################################
# This is where we start writing the C++ source file itself.
f = open('func.cpp','w')
f.write('#define NX {0}  // dim of state variable + number of parameters\n'.format(NX))
f.write('#define NS {0}  // number of stimuli\n'.format(NS))
f.write('\n')
f.write('using namespace alglib;\n')
f.write('\n')
f.write('void func_origin(real_1d_array &x, int it, real_2d_array &sti, real_1d_array &dxdt);\n')
f.write('void func_DF(real_1d_array &x, int it, real_2d_array &sti, real_2d_array &Jac);\n')
f.write('\n')
f.write('// {0} vector field\n'.format(Problem))
f.write('void func_origin(real_1d_array &x, int it, real_2d_array &sti, real_1d_array &dxdt)\n{\n')
for i in range(NX):
  f.write('    dxdt[{0}] = {1};\n'.format(i,outVF[i]))
f.write('}\n')
f.write('\n')
f.write('//{0} Jacobian matrix\n'.format(Problem))
f.write('void func_DF(real_1d_array &x, int it, real_2d_array &sti, real_2d_array &Jac)\n{\n')
for i in range(NX):
  for j in range(NX):
    f.write('    Jac[{0}][{1}] = {2};\n'.format(i,j,outDF[i*NX+j]))
f.write('}\n')

f.close()
print 'Generated func.cpp sucessfully!'
