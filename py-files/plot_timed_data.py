import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('../data/psi_total_out_1p2.txt_01')
data = np.resize(data, [65, 33])
plt.contourf(data)
plt.colorbar()
plt.show()
