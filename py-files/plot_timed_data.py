import os
import numpy as np
import matplotlib.pyplot as plt

files = [f for f in os.listdir() if os.path.isfile(f)]
for file in files:
    if file[0:14] == 'psi_total_out_':
        print(file)
        data = np.loadtxt(file)
        data = np.resize(data, [65, 33])
        plt.contourf(data)
        plt.colorbar()
        plt.show()
