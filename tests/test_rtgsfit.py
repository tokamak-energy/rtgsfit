
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

filename = '../data/12001000_RUN01_for_python_works.mat'
inputsname = '../data/12001000_RUN01_inputs.mat'

#@fixture(scope="module")
#def data():
#    return loadmat(filename, squeeze_me=True)


#def test_rm_coil_from_meas(data):
#    
#    n_coil = data['n_coil']
#    n_meas = data['n_meas']
#    g_meas_coil = data['g_meas_coil']
#    coil_curr = np.random.random_sample((n_coil, ))

#    meas = np.random.random_sample((n_meas, ))
#    
#    meas_no_coil = np.zeros((n_meas, ))
#    
#    meas_no_coil_gt = meas - g_meas_coil @ coil_curr    
#    
#    # https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-0/cblas-gemm-001.html
#    coil_curr, meas,  meas_no_coil = run_c_func(c_rtgsfit.rm_coil_from_meas, coil_curr, meas,  meas_no_coil)    
#         
#    assert np.all(np.isclose(meas_no_coil, meas_no_coil_gt))
#    

#def test_make_basis(data ):

#    n_z = data['n_z']
#    n_r = data['n_r']
#    n_coef = data['n_coef']
#    n_grid = data['n_grid']
#    dz = data['dz']
#    r_grid = data['r_grid']
#    inv_r_mu0 = data['inv_r_mu0']
#    
#    psi_norm = np.random.random_sample((n_z, n_r))
#    mask = np.ones((n_z, n_r), dtype=np.int64)
#    basis = np.zeros((n_coef, n_grid))                
#   
#    grad_row = FinDiff((0, dz, 1), acc=2)

#    dpsi_dz = grad_row(psi_norm)
#    
#    basis_gt = np.vstack(((1-psi_norm.flatten())*r_grid.flatten(), (1-psi_norm.flatten())*inv_r_mu0.flatten(), dpsi_dz.flatten()))
#                
#    psi_norm, mask, basis = run_c_func( 
#                c_rtgsfit.make_basis, psi_norm, mask, basis)   
#    
#    
#    assert np.all(np.isclose(basis, basis_gt))


#def test_rtgsfit(data):

data = loadmat(filename, squeeze_me=True)
inputs = loadmat(inputsname, squeeze_me=True)

n_z = data['n_z']
n_r = data['n_r']
r_grid = data['r_grid']
z_grid = data['z_grid']
meas = data['measurements']
coil_curr = data['coil_current'].astype(float)   
#flux = data['flux'] 
n_meas = data['n_meas']
n_coef = data['n_coef']
n_grid = data['n_grid']
#g_grid_meas = data['g_grid_meas']
#r_mu0_dz2 = data['r_mu0_dz2']
#n_coil = data['n_coil']
#g_grid_coil = data['g_grid_coil']
#weight = data['weight']
#n_xpt_max = data['n_xpt_max']
#n_lcfs_max = data['n_lcfs_max']


mask = np.ones((n_z, n_r), dtype=np.int64)  
psi_norm =  1 - np.exp(-10*(z_grid**2)) * np.exp(-10*(r_grid - 0.5)**2);
meas_no_coil = np.zeros(n_meas,)
basis = np.zeros(n_coef*n_grid)
g_coef_meas = np.zeros(n_coef*n_meas)
curr_dens = np.zeros(n_grid,)
flux_pls = np.zeros(n_grid,)
flux_coil = np.zeros(n_grid,)

meas, coil_curr, psi_norm, mask = run_c_func( 
            c_rtgsfit.rtgsfit, meas, coil_curr, psi_norm, mask)   
            
plt.imshow(psi_norm)
plt.show()

#coil_curr, meas,  meas_no_coil = run_c_func(c_rtgsfit.rm_coil_from_meas, coil_curr, meas,  meas_no_coil)    
# 
#psi_norm, mask, basis = run_c_func(c_rtgsfit.make_basis, psi_norm, mask, basis)   


#                
#_, _, _, n_coef, n_meas, n_grid, _, basis, n_grid, g_grid_meas, n_meas, _, g_coef_meas, n_meas = run_c_func( 
#            c_cblas.cblas_dgemm, c_cblas.CblasRowMajor, c_cblas.CblasNoTrans, 
#            c_cblas.CblasNoTrans, n_coef, n_meas, n_grid, 1.0, basis, 
#            n_grid, g_grid_meas, n_meas, 0.0, g_coef_meas, n_meas)   

#meas_no_coil = meas_no_coil * weight  
#coef_est = meas_no_coil.copy()      

#g_meas_coef = np.reshape(g_coef_meas,(n_coef, n_meas)).T

#coef_gt,_,_,_ =np.linalg.lstsq(g_meas_coef, meas_no_coil)


#    
#g_coef_meas, coef_est = run_c_func(c_dgelss.dgelss2, g_coef_meas, meas_no_coil)

#_,_,_,_,_, basis, _, coef_est,_,_, curr_dens, _ = run_c_func(c_cblas.cblas_dgemv, 
#        c_cblas.CblasColMajor, c_cblas.CblasNoTrans, n_grid, n_coef, 1.0, basis, 
#        n_grid, coef_est[:n_coef], 1, 0.0, curr_dens, 1)      
#        
#        

##plt.imshow(np.reshape(curr_dens, (n_z, n_r)))
##plt.show()

#source = -curr_dens * r_mu0_dz2.flatten()

#source, flux_pls = run_c_func(c_poisson_solver.poisson_solver, source, flux_pls)    

#_, _,_,_, _, _,_, coil_curr, _,_, flux_coil, _ =  run_c_func(
#        c_cblas.cblas_dgemv, c_cblas.CblasRowMajor, c_cblas.CblasNoTrans, n_grid, 
#        n_coil, 1.0, g_grid_coil, n_coil, coil_curr, 1, 0.0, flux_coil, 1)   

#flux_total = flux_coil + 2 * np.pi * flux_pls;

#print(flux_total[::10])

#opt_n = np.array([0])
#opt_r = np.zeros(n_xpt_max)
#opt_z = np.zeros(n_xpt_max)
#opt_flux = np.zeros(n_xpt_max)
#xpt_n = np.array([0])
#xpt_r = np.zeros(n_xpt_max)
#xpt_z = np.zeros(n_xpt_max)
#xpt_flux = np.zeros(n_xpt_max)


#flux_total, opt_r, opt_z, opt_flux, opt_n, xpt_r, xpt_z, xpt_flux, xpt_n = \
#        run_c_func(c_find_x_point.find_null_in_gradient, flux_total, opt_r, opt_z, 
#        opt_flux, opt_n, xpt_r, xpt_z, xpt_flux, xpt_n)


        
        
#    printf("j\n");    
#    // find x point & opt
#    find_null_in_gradient(flux_total, opt_r, opt_z, opt_flux, &opt_n, 
#            xpt_r, xpt_z, xpt_flux, &xpt_n);

#    printf("k\n");    
#    // select opt
#    i_opt = max_idx(opt_n, opt_flux);    
#    axis_flux = opt_flux[i_opt];
#    axis_r = opt_r[i_opt];
#    axis_z = opt_z[i_opt];  

#    printf("l\n");
#    printf("%d\n", xpt_n);
#/*    for (i_grid=0; i_grid<xpt_n; i_grid++)*/
#/*    {*/
#/*        printf("%f\n", xpt_flux[i_grid]);*/
#/*    }*/
#    // select xpt      
#    i_xpt = max_idx(xpt_n, xpt_flux);
#    printf("%d\n", i_xpt);
#    xpt_flux_max = xpt_flux[i_xpt];      
#    lcfs_flux = FRAC * xpt_flux_max + (1-FRAC)*axis_flux;

#    printf("m\n");    
#    // extract LCFS
#    find_lcfs_rz(flux_total, lcfs_flux, lcfs_r, lcfs_z, &lcfs_n);         

#    printf("m\n");            
#    // extract inside of LCFS
#    inside_lcfs(axis_r, axis_z, lcfs_r, lcfs_z, lcfs_n, mask);

#    printf("n\n");    
#    // normalise total psi                                
#    normalise_flux(flux_norm, flux_total, lcfs_flux, axis_flux, mask);

#plt.contour(np.reshape(flux_pls, (n_z, n_r)))
#plt.show()
 
#plt.plot(g_meas_coef @ coef_gt / meas_no_coil, 'x')
#plt.plot(meas_no_coil/meas_no_coil, '+')
#plt.plot(g_meas_coef @ coef_est[:n_coef] / meas_no_coil, '.')
#plt.show()
                            
#    print(f"{meas.shape}, {coil_curr.shape}, {psi_norm.shape}")
#    breakpoint()
#    psi_norm, mask, basis = run_c_func( 
#                c_rtgsfit.rtgsfit, meas, coil_curr, psi_norm, mask)   
                
#    plt.contour(psi_norm)

#if __name__ == "__main__":

#    n_meas = c_constants.N_MEAS.value
#    n_coil = c_constants.N_COIL.value
#    
#    n_z = c_constants.N_Z.value
#    n_r = c_constants.N_R.value   
#    mylib = ctypes.CDLL('./constants.so')
#    double_array_type = (ctypes.c_double * n_r)
#    double_array_ptr = ctypes.cast(mylib.R_VEC, ctypes.POINTER(double_array_type))
#    print(list(double_array_ptr.contents))
#    n_coef = c_constants.N_COEF.value

#    test_rm_coil_from_meas(n_meas, n_coil)
##    test_make_basis(n_z, n_r, n_coef)
