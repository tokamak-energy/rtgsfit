
import os
import ctypes
import numpy as np
import matplotlib.pyplot as plt
import pytest
from pytest_cases import fixture
import c_librtgsfit
import c_libconstants
from findiff import FinDiff
from scipy.io import loadmat
import sys
sys.path.append("../utility")
from py_to_ctypes import run_c_func

#filename = '../data/12001000_RUN01_for_python_fj.mat'
#inputsname = '../data/12001000_RUN01_inputs_fj.mat'

filename = '../data/12001000_RUN01_for_c.mat'
inputsname = '../data/12001000_RUN01_inputs.mat'

data = loadmat(filename, squeeze_me=True)
inputs = loadmat(inputsname, squeeze_me=True)

n_z = data['n_z']
n_r = data['n_r']
r_grid = data['r_grid']
z_grid = data['z_grid']
r_vec = data['r_vec']
z_vec = data['z_vec']
meas = inputs['measurements']
coil_curr = inputs['coil_current'].astype(float)   
n_meas = data['n_meas']
n_coef = data['n_coef']
n_grid = data['n_grid']
n_lcfs_max = data['n_lcfs_max']
n_lcfs_max = data['n_lcfs_max']

mask = np.ones((n_grid, ), dtype=np.int64)  
psi_norm =  1 - np.exp(-10*(z_grid**2)) * np.exp(-10*(r_grid - 0.5)**2);

#plt.imshow(np.reshape(psi_norm, (65, 33)))
#plt.show()
meas_no_coil = np.zeros(n_meas,)
basis = np.zeros(n_coef*n_grid)
g_coef_meas = np.zeros(n_coef*n_meas)
curr_dens = np.zeros(n_grid,)
flux_pls = np.zeros(n_grid,)
flux_coil = np.zeros(n_grid,)
psi_total = np.zeros(psi_norm.shape)
error = np.zeros((1,))
r_lcfs = np.zeros(n_lcfs_max, )
z_lcfs = np.zeros(n_lcfs_max, )
coef = np.zeros(n_coef, )
n_lcfs = np.zeros((1,), dtype=np.int64)

fig, ax = plt.subplots(1, 4, sharex=True, sharey=True)
ax = ax.flatten()

meas_orig = meas.copy()
coil_curr_orig = coil_curr.copy()

np.savetxt('../data/flux_norm.txt', psi_norm)
np.savetxt('../data/mask.txt', mask, fmt='%d')
np.savetxt('../data/psi_total.txt', psi_total)
np.savetxt('../data/error.txt', error)
np.savetxt('../data/meas.txt', meas_orig)
np.savetxt('../data/coil_curr.txt', coil_curr_orig)

for ii in range(16):
#    print(ii)
    meas, coil_curr, psi_norm, mask, psi_total, error, r_lcfs, z_lcfs, n_lcfs, coef = run_c_func( 
                c_librtgsfit.rtgsfit, meas_orig, coil_curr_orig, psi_norm, mask, psi_total, error, r_lcfs, z_lcfs, n_lcfs, coef) 
                
    meas = meas.astype(float)
    coil_curr = coil_curr.astype(float)
    psi_norm = psi_norm.astype(float)
    mask = mask.astype(np.int64)
    psi_total = psi_total.astype(float) 
    error = error.astype(float)   
    r_lcfs = r_lcfs.astype(float)
    z_lcfs = z_lcfs.astype(float)
    n_lcfs = n_lcfs.astype(np.int64)
    coef = coef.astype(float)
        
    if ii % 5 == 0:
        jj = int(ii/5)
        print(ii, coef)
        out = ax[jj].contourf(r_vec, z_vec, np.reshape(psi_total, (n_z, n_r)), 20)
        ax[jj].plot(r_lcfs[:n_lcfs[0]], z_lcfs[:n_lcfs[0]], 'k.')    
        ax[jj].set_aspect('equal')
        ax[jj].set_title(f"{ii}: {error[0]:.5e}")

plt.colorbar(out, ax = ax.ravel().tolist(), shrink=0.95)

plt.show()

