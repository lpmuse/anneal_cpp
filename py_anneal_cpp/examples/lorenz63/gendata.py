# generate data from Lorenz 63 model
import numpy as np
from scipy.integrate import ode

import matplotlib.pyplot as plt

# define the model ODE's
def L63(t, x, dxdt, s, r, b):
    dxdt[0] = s * (x[1] - x[0])
    dxdt[1] = x[0] * (r - x[2]) - x[1]
    dxdt[2] = x[0] * x[1] - b * x[2]
    return dxdt

# run out the transient
# integration parameters
t0 = 0.0
tf = 40.0
dt = 0.05
NT = (tf - t0) / dt + 1
times = np.linspace(t0, tf, NT)

x0 = np.array([1,1,1])
dxdt = np.zeros(3)

sol = np.zeros((NT,3))
sol[0] = x0

s = 10.0
r = 28.0
b = 8.0/3.0

itg = ode(L63)
itg.set_integrator('dop853')
itg.set_initial_value(x0, t0)
itg.set_f_params(dxdt, s, r, b)

for i,tnext in enumerate(times[1:]):
    itg.integrate(tnext)
    sol[i+1] = itg.y

# now integrate and save the data
t0 = tf
tf += 10.0
dt = 0.05
NT = (tf - t0) / dt + 1
print NT
times = np.linspace(t0, tf, NT)

x0 = sol[-1]
dxdt = L63(t0, x0, dxdt, s, r, b)

sol = np.zeros((NT,3))
sol[0] = x0

for i,tnext in enumerate(times[1:]):
    itg.integrate(tnext)
    sol[i+1] = itg.y


fig,ax = plt.subplots(3,1,sharex=True)
fig.set_tight_layout(True)
ax[0].plot(times, sol[:,0])
ax[1].plot(times, sol[:,1])
ax[2].plot(times, sol[:,2])
plt.show()

# add some noise
xmin, xmax = (-20.0, 20.0)
ymin, ymax = (-25.0, 25.0)
zmin, zmax = (0.0, 50.0)

xdiff = xmax - xmin
ydiff = ymax - ymin
zdiff = zmax - zmin

std = 0.01*np.array([xdiff, ydiff, zdiff])

sol[:,0] += std[0]*np.random.randn(NT)
sol[:,1] += std[1]*np.random.randn(NT)
sol[:,2] += std[2]*np.random.randn(NT)

fig,ax = plt.subplots(3,1,sharex=True)
fig.set_tight_layout(True)
ax[0].plot(times, sol[:,0])
ax[1].plot(times, sol[:,1])
ax[2].plot(times, sol[:,2])
plt.show()

# save the data to disk
np.savetxt('twin_data_original.dat', sol)
