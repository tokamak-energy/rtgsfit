
import ctypes
import numpy as np
import matplotlib.pyplot as plt
import pytest
from pytest_cases import fixture
from py_to_ctypes import run_c_func
import c_rtgsfit
import c_constants
from findiff import FinDiff

def test_rm_coil_from_meas(n_meas=20, n_coil=5):
    
    print(n_meas, n_coil)
    coil_curr = np.random.random_sample((n_coil, ))

    meas = np.random.random_sample((n_meas, ))
    
    meas_no_coil = np.zeros((n_meas, ))

#    g_meas_coil = np.ctypeslib.as_array(c_constants.G_MEAS_COIL, (n_meas, n_coil))
#    
#    meas_no_coil_gt = meas - g_meas_coil @ coil_curr
    
    # https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-0/cblas-gemm-001.html
    coil_curr, meas,  meas_no_coil = run_c_func(c_rtgsfit.rm_coil_from_meas, coil_curr, meas,  meas_no_coil)    
    
    breakpoint()
#     
#    assert np.all(np.isclose(meas_no_coil, meas_no_coil_gt))
    

def test_make_basis(n_z=30, n_r=15, n_coef=3):

    assert n_coef == 3
    n_grid = n_z*n_r
    psi_norm = np.random.random_sample((n_z, n_r))
    basis = np.zeros((n_coef, n_grid))                
   
    grad_row = FinDiff((0, d_row, 1), acc=2)

    dpsi_dz = grad_row(psi_norm)
    
    basis_gt = np.vstack(((1-psi_norm.flatten())*r_grid, (1-psi_norm.flatten())*inv_r_mu0, dpsi_dz.flatten()))
                
    n_z, n_r, n_grid, d_row, psi_norm, r_grid, inv_r_mu0, basis = run_c_func( 
                c_rtgsfit.make_basis, psi_norm, basis)   
    
    assert np.all(np.isclose(basis, basis_gt))
                   
                   
if __name__ == "__main__":

    n_meas = c_constants.N_MEAS.value
    n_coil = c_constants.N_COIL.value
    
    n_z = c_constants.N_Z.value
    n_r = c_constants.N_R.value   
    mylib = ctypes.CDLL('./constants.so')
    double_array_type = (ctypes.c_double * n_r)
    double_array_ptr = ctypes.cast(mylib.R_VEC, ctypes.POINTER(double_array_type))
    print(list(double_array_ptr.contents))
    n_coef = c_constants.N_COEF.value

    test_rm_coil_from_meas(n_meas, n_coil)
#    test_make_basis(n_z, n_r, n_coef)
