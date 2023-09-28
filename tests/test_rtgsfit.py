
import os
import ctypes
import numpy as np
import matplotlib.pyplot as plt
import pytest
from pytest_cases import fixture
import c_rtgsfit
import c_constants
from findiff import FinDiff
from scipy.io import loadmat
import sys
sys.path.append("../utility")
from py_to_ctypes import run_c_func

filename = '../data/12001000_RUN01_for_python.mat'
inputsname = '../data/12001000_RUN01_inputs.mat'

data = loadmat(filename, squeeze_me=True)
inputs = loadmat(inputsname, squeeze_me=True)

n_z = data['n_z']
n_r = data['n_r']
r_grid = data['r_grid']
z_grid = data['z_grid']
meas = inputs['measurements']
coil_curr = inputs['coil_current'].astype(float)   
n_meas = data['n_meas']
n_coef = data['n_coef']
n_grid = data['n_grid']

mask = np.ones((n_grid, ), dtype=np.int64)  
psi_norm =  1 - np.exp(-10*(z_grid**2)) * np.exp(-10*(r_grid - 0.5)**2);
meas_no_coil = np.zeros(n_meas,)
basis = np.zeros(n_coef*n_grid)
g_coef_meas = np.zeros(n_coef*n_meas)
curr_dens = np.zeros(n_grid,)
flux_pls = np.zeros(n_grid,)
flux_coil = np.zeros(n_grid,)
psi_total = np.zeros(psi_norm.shape)
error = np.zeros((1,))

fig, ax = plt.subplots(3,6, sharex=True, sharey=True)
ax = ax.flatten()

meas_orig = meas.copy()
coil_curr_orig = coil_curr.copy()

np.savetxt('../data/flux_norm.txt', psi_norm)
np.savetxt('../data/mask.txt', mask, fmt='%d')
np.savetxt('../data/psi_total.txt', psi_total)
np.savetxt('../data/error.txt', error)
np.savetxt('../data/meas.txt', meas_orig)
np.savetxt('../data/coil_curr.txt', coil_curr_orig)

for ii in range(18):
    print(ii)
    meas, coil_curr, psi_norm, mask, psi_total, error = run_c_func( 
                c_rtgsfit.rtgsfit, meas_orig, coil_curr_orig, psi_norm, mask, psi_total, error) 
                
    meas = meas.astype(float)
    coil_curr = coil_curr.astype(float)
    psi_norm = psi_norm.astype(float)
    mask = mask.astype(np.int64)
    psi_total = psi_total.astype(float) 
    error = error.astype(float)   
      
    out = ax[ii].contourf(np.reshape(psi_total, (n_z, n_r)), 20)
    plt.colorbar(out, ax = ax[ii])    
    ax[ii].set_title(f"{ii}: {error[0]:.5e}")


plt.show()

