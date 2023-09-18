import ctypes
import numpy as np
import matplotlib.pyplot as plt
import pytest
import c_solve_tria
from scipy.io import loadmat
from curr_dens import curr_rand
from grad_shafranov import gs_finite_diff, plu, plu_band
from py_to_ctypes import run_c_func




@pytest.mark.parametrize("seed", np.arange(0, 2))            
def test_permute(seed, datafile='12001000_RUN01_for_python.mat', thresh=1e-10, show=False):
  
    data = loadmat(datafile, squeeze_me=True)
    
    n_r = data['n_r']
    n_z = data['n_z']
    lower = data['lower']
    perm = data['perm']
    upper = data['upper']
  
    b_vec = curr_rand(n_r, n_z, seed).flatten()    
    truth = perm @ b_vec 
    estimate = np.zeros(n_r*n_z)    

    b_vec, estimate = run_c_func(c_solve_tria.permute, b_vec, estimate) 

    if show:
        truth = np.reshape(truth, (n_z, n_r))
        estimate = np.reshape(estimate, (n_z, n_r))   
         
        fig, ax = plt.subplots(2, 2)
        ax[0, 0].imshow(truth)
        ax[0, 1].plot(truth, '+')
        ax[0, 1].set_prop_cycle(None)
        ax[0, 1].plot(estimate, 'x')
        ax[1, 0].imshow(estimate)
        ax[1, 1].imshow(np.abs(truth-estimate))
        plt.show()
 
    for est, tru in zip(estimate.flatten(), truth.flatten()):
        assert(np.abs(est-tru) < thresh)
        
@pytest.mark.parametrize("seed", np.arange(0, 2))            
def test_back_sub_upper_short(seed, datafile='12001000_RUN01_for_python.mat', thresh=1e-10, show=False):
    data = loadmat(datafile, squeeze_me=True)
    
    n_r = data['n_r']
    n_z = data['n_z']
    lower = data['lower']
    perm = data['perm']
    upper = data['upper']
  
    b_vec = curr_rand(n_r, n_z, seed).flatten()    
    truth = np.linalg.solve(upper, b_vec) 
    estimate = np.zeros(n_r*n_z)    

    b_vec, estimate = run_c_func(c_solve_tria.back_sub_upper_short, b_vec, estimate) 

    if show:
        truth = np.reshape(truth, (n_z, n_r))
        estimate = np.reshape(estimate, (n_z, n_r))   
         
        fig, ax = plt.subplots(2, 2)
        ax[0, 0].imshow(truth)
        ax[0, 1].plot(truth, '+')
        ax[0, 1].set_prop_cycle(None)
        ax[0, 1].plot(estimate, 'x')
        ax[1, 0].imshow(estimate)
        ax[1, 1].imshow(np.abs(truth-estimate))
        plt.show()
 
    for est, tru in zip(estimate.flatten(), truth.flatten()):
        assert(np.abs(est-tru) < thresh)
        
@pytest.mark.parametrize("seed", np.arange(0, 2))            
def test_back_sub_lower_short(seed, datafile='12001000_RUN01_for_python.mat', thresh=1e-10, show=False):
    data = loadmat(datafile, squeeze_me=True)
    
    n_r = data['n_r']
    n_z = data['n_z']
    lower = data['lower']
    perm = data['perm']
    upper = data['upper']
  
    b_vec = curr_rand(n_r, n_z, seed).flatten()    
    truth = np.linalg.solve(lower, b_vec) 
    estimate = np.zeros(n_r*n_z)    

    b_vec, estimate = run_c_func(c_solve_tria.back_sub_lower_short, b_vec, estimate) 

    if show:
        truth = np.reshape(truth, (n_z, n_r))
        estimate = np.reshape(estimate, (n_z, n_r))   
         
        fig, ax = plt.subplots(2, 2)
        ax[0, 0].imshow(truth)
        ax[0, 1].plot(truth, '+')
        ax[0, 1].set_prop_cycle(None)
        ax[0, 1].plot(estimate, 'x')
        ax[1, 0].imshow(estimate)
        ax[1, 1].imshow(np.abs(truth-estimate))
        plt.show()
 
    for est, tru in zip(estimate.flatten(), truth.flatten()):
        assert(np.abs(est-tru) < thresh)
        
@pytest.mark.parametrize("seed", np.arange(0, 2))            
def test_solve_tria(seed, datafile='12001000_RUN01_for_python.mat', thresh=1e-10, show=False):
  
    data = loadmat(datafile, squeeze_me=True)
    
    n_r = data['n_r']
    n_z = data['n_z']
    lower = data['lower']
    perm = data['perm']
    upper = data['upper']
  
    b_vec = curr_rand(n_r, n_z, seed).flatten()    
    truth = np.linalg.solve(lower, perm @ b_vec) 
    truth = np.linalg.solve(upper, truth)   
    estimate = np.zeros(n_r*n_z)    

    b_vec, estimate = run_c_func(c_solve_tria.solve_tria, b_vec, estimate) 

    if show:
        truth = np.reshape(truth, (n_z, n_r))
        estimate = np.reshape(estimate, (n_z, n_r))   
         
        fig, ax = plt.subplots(2, 2)
        ax[0, 0].imshow(truth)
        ax[0, 1].plot(truth, '+')
        ax[0, 1].set_prop_cycle(None)
        ax[0, 1].plot(estimate, 'x')
        ax[1, 0].imshow(estimate)
        ax[1, 1].imshow(np.abs(truth-estimate))
        plt.show()
 
    for est, tru in zip(estimate.flatten(), truth.flatten()):
        assert(np.abs(est-tru) < thresh)

if __name__ == '__main__':
    test_permute(0, show=True)
    test_back_sub_upper_short(0, show=True)
    test_back_sub_lower_short(0, show=True)    
    test_solve_tria(0, show=True)    
    
#@pytest.mark.parametrize("n_r", 2**np.arange(3, 6))
#@pytest.mark.parametrize("n_z", 2**np.arange(3, 7))
#@pytest.mark.parametrize("seed", np.arange(0, 2))            
#def test_solve_tria(n_r, n_z, seed, r_min=0.2, r_max=1.0, z_min=-1.0, z_max=1.0,
#        thresh=1.0e-10, show=False):
#        
#    r_vec = np.linspace(r_min, r_max, n_r)
#    z_vec = np.linspace(z_min, z_max, n_z)
#   
#    sys_mat = gs_finite_diff(r_vec, z_vec)
#    perm, lower, upper = plu(sys_mat)
#    perm_idx, lower_band, upper_band = plu_band(perm, lower, upper, n_r)
#    b_vec = curr_rand(n_r, n_z, seed).flatten()
#    
#    truth = np.linalg.solve(lower, perm @ b_vec) 
#    truth = np.linalg.solve(upper, truth)
#    
#    estimate = np.zeros(n_r*n_z)    
#    n_grid = n_r*n_z

#    n_grid, n_r, lower_band, upper_band, b_vec, perm_idx, estimate = \
#            run_c_func(c_solve_tria.solve_tria, n_grid, n_r, lower_band, 
#            upper_band, b_vec, perm_idx, estimate) 

#    if show:
#        truth = np.reshape(truth, (n_z, n_r))
#        estimate = np.reshape(estimate, (n_z, n_r))   
#         
#        fig, ax = plt.subplots(2, 2)
#        ax[0, 0].imshow(truth)
#        ax[0, 1].plot(truth, '+')
#        ax[0, 1].set_prop_cycle(None)
#        ax[0, 1].plot(estimate, 'x')
#        ax[1, 0].imshow(estimate)
#        ax[1, 1].imshow(np.abs(truth-estimate))
#        plt.show()
# 
#    for est, tru in zip(estimate.flatten(), truth.flatten()):
#        assert(np.abs(est-tru) < thresh)           
