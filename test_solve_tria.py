import ctypes
import numpy as np
import matplotlib.pyplot as plt
import pytest
import c_solve_tria
from curr_dens import curr_rand
from grad_shafranov import gs_finite_diff, plu, plu_band
from py_to_ctypes import run_c_func


@pytest.mark.parametrize("n_r", 2**np.arange(3, 6))
@pytest.mark.parametrize("n_z", 2**np.arange(3, 7))
@pytest.mark.parametrize("seed", np.arange(0, 2))            
def test_solve_tria(n_r, n_z, seed, r_min=0.2, r_max=1.0, z_min=-1.0, z_max=1.0,
        thresh=1.0e-10, show=False):
        
    r_vec = np.linspace(r_min, r_max, n_r)
    z_vec = np.linspace(z_min, z_max, n_z)
   
    sys_mat = gs_finite_diff(r_vec, z_vec)
    perm, lower, upper = plu(sys_mat)
    perm_idx, lower_band, upper_band = plu_band(perm, lower, upper, n_r)
    b_vec = curr_rand(n_r, n_z, seed).flatten()
    
    truth = np.linalg.solve(lower, perm @ b_vec) 
    truth = np.linalg.solve(upper, truth)
    
    estimate = np.zeros(n_r*n_z)    
    n_grid = n_r*n_z

    n_grid, n_r, lower_band, upper_band, b_vec, perm_idx, estimate = \
            run_c_func(c_solve_tria.solve_tria, n_grid, n_r, lower_band, 
            upper_band, b_vec, perm_idx, estimate) 

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
           
