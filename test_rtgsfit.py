
import ctypes
import numpy as np
import matplotlib.pyplot as plt
import pytest
from pytest_cases import fixture
from py_to_ctypes import run_c_func
import c_rtgsfit
from findiff import FinDiff

def test_rm_coil_from_meas(n_meas=20, n_coil=5):
    

    g_meas_coil = np.random.random_sample((n_meas, n_coil))

    coil_curr = np.random.random_sample((n_coil, ))

    meas = np.random.random_sample((n_meas, ))
    
    meas_no_coil = np.zeros((n_meas, ))

    meas_no_coil_gt = meas - g_meas_coil @ coil_curr
    
    # https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-0/cblas-gemm-001.html
    n_meas, n_coil, g_meas_coil, coil_curr, meas,  meas_no_coil = run_c_func( 
                c_rtgsfit.rm_coil_from_meas, n_meas, n_coil, g_meas_coil,
                coil_curr, meas,  meas_no_coil)    
     
    assert np.all(np.isclose(meas_no_coil, meas_no_coil_gt))
    

def test_make_basis(n_row=30, n_col=15):

    n_grid = n_row*n_col
    d_row = np.random.random_sample(1)[0]
    psi_norm = np.random.random_sample((n_row, n_col))
    r_grid =  np.random.random_sample((n_grid, ))
    inv_r_mu0 =  np.random.random_sample((n_grid, ))
    basis = np.zeros((3, n_grid))                
   
    grad_row = FinDiff((0, d_row, 1), acc=2)

    dpsi_dz = grad_row(psi_norm)
    
    basis_gt = np.vstack(((1-psi_norm.flatten())*r_grid, (1-psi_norm.flatten())*inv_r_mu0, dpsi_dz.flatten()))
                
    n_row, n_col, n_grid, d_row, psi_norm, r_grid, inv_r_mu0, basis = run_c_func( 
                c_rtgsfit.make_basis, n_row, n_col, n_grid, d_row, psi_norm, r_grid, inv_r_mu0, basis)   
    
    assert np.all(np.isclose(basis, basis_gt))
                   
                   
if __name__ == "__main__":
    test_rm_coil_from_meas()
    test_make_basis()
