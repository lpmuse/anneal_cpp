# rescale data according to the transformation
#     x --> (x - x_min) / (x_max - x_min)
# where x_min and x_max are set by the dynamical voltage range (determined
# roughly from the data).

import numpy as np

# import the original data
original_data = np.loadtxt('twin_data_original.dat')

# set up the rescaling
xmin, xmax = (-20.0, 20.0)
ymin, ymax = (-25.0, 25.0)
zmin, zmax = (0.0, 50.0)

rescaled_data = np.copy(original_data)
rescaled_data[:,0] = (rescaled_data[:,0] - xmin) / (xmax - xmin)
rescaled_data[:,1] = (rescaled_data[:,1] - xmin) / (xmax - xmin)
rescaled_data[:,2] = (rescaled_data[:,2] - xmin) / (xmax - xmin)

np.savetxt('twin_data.dat', rescaled_data)
