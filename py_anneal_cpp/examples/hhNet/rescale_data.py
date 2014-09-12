# rescale voltage data using the transformation
#     V --> (V - V_min) / (V_max - V_min)
# where V_min and V_max are set by the dynamical voltage range (determined
# roughly from the data).

import numpy as np

# import the original data
original_data = np.loadtxt('twin_data_original.dat')

# set up the rescaling
Vmin = -80.0
Vmax = 20.0

rescaled_data = (original_data - Vmin) / (Vmax - Vmin)

np.savetxt('twin_data.dat', rescaled_data)
