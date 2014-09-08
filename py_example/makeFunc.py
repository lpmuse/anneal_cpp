import sympy as sym 
from sympy import *

#  Opening and reading the text file with the vector field information
print "Importing equations.txt\n"
file = open('equations.txt','r')
temp=[]  # Array to hold equations.txt information
for line in file:
  if line.startswith('#'): # Pound used as comment in text file
     continue
  elif line.startswith('\\'): # In case file has UTF-8 markers
     continue
  else:
     temp.append(line)
file.close()

h=[]  # Array to hold unformatted equations.txt information
for i in range(len(temp)):
  temp1=temp[i].rstrip( )
  h.append(temp1)

# Problem name
Problem = h[0]
tmp = h[1].split(',', 1 );
NX = int(tmp[0])
#ND = int(tmp[1])
NS = int(tmp[1])
print 'Processing ',Problem, 'model\n'
mytemp = 'syms=['
for i in range(NX):
    mytemp2 = "%s=Symbol('%s')" % (h[2+NX+i],h[2+NX+i])
    exec mytemp2
    mytemp = '%s %s,' % (mytemp, h[2+NX+i])
mytemp='%s]' % mytemp[:-1]
exec mytemp

if NS>0:
	mytemp = 'stims=['
	for i in range(NS):
    		mytemp2 = "%s=Symbol('%s')" % (h[2+2*NX+i],h[2+2*NX+i])
    		exec mytemp2
    		mytemp = '%s %s,' % (mytemp, h[2+2*NX+i])
	mytemp='%s]' % mytemp[:-1]
	exec mytemp

mytemp = 'VF=Matrix(1,%d,[' % NX
for i in range(NX):
    mytemp = '%s %s,' % (mytemp, h[2+i])
mytemp='%s])' % mytemp[:-1]
exec mytemp
Jac = VF.jacobian(syms)



outVF=[]
for i in range(NX):
    tmp = '%s' % ccode(VF[i])
    for j in range(NX):
        tmp=tmp.replace('%s'%syms[j], 'x[%d]'%j)
    for j in range(NS):
        tmp=tmp.replace('%s'%stims[j], 'sti[it][%d]'%j)
    outVF.append(tmp)

    
outDF=[]
for i in range(NX):
    for j in range(NX):
        tmp = '%s' % ccode(Jac[i,j])
        for k in range(NX):
            tmp=tmp.replace('%s'%syms[k], 'x[%d]'%k)
        for k in range(NS):
            tmp=tmp.replace('%s'%stims[k], 'sti[it][%d]'%k)
        outDF.append(tmp)
        
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
