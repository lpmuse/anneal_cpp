import numpy as np
import matplotlib.pyplot as plt

NT = 1001
D = 4

path = np.loadtxt('path/D4_M1_PATH0.dat')

plt.plot(path[28,1::4])
plt.ylim(0,1)
plt.show()
