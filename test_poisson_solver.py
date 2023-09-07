import ctypes
import numpy as np
import matplotlib.pyplot as plt
import c_poisson_solver
import pytest
from scipy.constants import mu_0
from findiff import FinDiff
from grad_shafranov import gs_finite_diff, plu, plu_band    
from curr_dens import curr_ellipse
from greens import greens_grid, greens_bound
from py_to_ctypes import run_c_func

    
def hagenow(lower, upper, perm, source, inv_r_mu0, g_bound, n_z, n_r, dz, dr):

    psi_zero = np.linalg.solve(lower, perm @ source) 
    psi_zero = np.linalg.solve(upper, psi_zero)

    psi_zero = np.reshape(psi_zero, (n_z, n_r))

    grad_row = FinDiff((0, dz, 1), acc=2)
    truth_row = grad_row(psi_zero)

    grad_col = FinDiff((1, dr, 1), acc=2)
    truth_col = grad_col(psi_zero)        

    dpsi_bound = np.concatenate([truth_col[:, 0]*dz, truth_row[0, :]*dr, 
            -truth_col[:, -1]*dz, -truth_row[-1, :]*dr])
            
    dpsi_hagen = inv_r_mu0 * dpsi_bound

    psi_bound = g_bound @ dpsi_hagen
    
    return psi_bound


@pytest.mark.parametrize("n_r", 2**np.arange(4, 7))
@pytest.mark.parametrize("n_z", 2**np.arange(4, 7))
def test_hagenow(n_r, n_z, r_min=0.11, r_max=1.07, z_min=-0.96, z_max=0.96, 
        r0 = 0.5, z0 = 0.0, a_minor = 0.2, sigma_r = 1.0, sigma_z = 2.0, 
        thresh=1.0e-10, show=False):    

    n_ele = n_r*n_z
    n_bound = 2*(n_r + n_z)

    r_vec = np.linspace(r_min, r_max, n_r)
    z_vec = np.linspace(z_min, z_max, n_z)
    r_grid, z_grid = np.meshgrid(r_vec, z_vec)     
    dz = z_vec[1] - z_vec[0]
    dr = r_vec[1] - r_vec[0]   
    r_ltrb = np.concatenate((r_grid[:, 0], r_grid[0, :], r_grid[:, -1], r_grid[-1, :]))
    z_ltrb = np.concatenate((z_grid[:, 0], z_grid[0, :], z_grid[:, -1], z_grid[-1, :]))
    inv_r_mu0 = 1.0/(mu_0 * r_ltrb)    
    
    sys_mat = gs_finite_diff(r_vec, z_vec)
    perm, lower, upper = plu(sys_mat)
    perm_idx, lower_band, upper_band = plu_band(perm, lower, upper, n_r)         
    curr_dens = curr_ellipse(r_vec, z_vec, r0, z0, sigma_r, sigma_z, a_minor)

    source = -curr_dens.flatten()*mu_0*r_grid.flatten()*(dz**2)
    
    g_grid = greens_grid(r_grid, z_grid, dr, dz)
    g_bound = greens_bound(r_ltrb, z_ltrb, dr, dz, n_r, n_z)

    truth_psi = g_grid @ curr_dens.flatten() * dr * dz
    truth_psi = np.reshape(truth_psi, (n_z, n_r))
    truth_psi_ltrb = np.concatenate((truth_psi[:, 0], truth_psi[0, :], 
            truth_psi[:, -1], truth_psi[-1, :]))

    py_psi_ltrb = hagenow(lower, upper, perm, source, inv_r_mu0, g_bound, n_z, 
            n_r, dz, dr)

    psi = np.zeros(n_ele, )
    dpsi_ltrb = np.zeros(n_bound, )
    c_psi_ltrb = np.zeros(n_bound, )

    n_r, n_z, n_ele, lower_band, upper_band, source, perm_idx, n_bound, g_bound, \
            inv_r_mu0, dz, dr, psi, psi_bound, dpsi_ltrb = run_c_func( \
            c_poisson_solver.hagenow_bound, n_r, n_z, n_ele, lower_band, \
            upper_band, source, perm_idx, n_bound, g_bound,
            inv_r_mu0, dz, dr, psi, c_psi_ltrb, dpsi_ltrb)
            
    if show:
        fig, ax = plt.subplots(1, 3)
        ax[0].imshow(curr_dens)
        ax[1].imshow(np.reshape(psi,(n_z, n_r)))
        ax[2].plot(truth_psi_ltrb, 'o')
        ax[2].plot(py_psi_ltrb, 'x')
        ax[2].plot(c_psi_ltrb, '+')
        ax[2].legend(['truth', 'py-hag', 'c-hag'])
        plt.show()        
    
    for ii, (est, py_est, tru) in enumerate(zip(psi_bound, py_psi_ltrb, truth_psi_ltrb)):
        assert np.abs((est-py_est)/py_est) < thresh, f"{ii} {est} {py_est} {np.abs((est-py_est)/py_est)} {thresh}"
        assert np.abs((est - tru)/tru) < (2.0/n_z), f"{ii} {est} {tru} {np.abs((est - tru)/tru)} {(1.0/np.sqrt(n_r*n_z))}"
        
        
        

            
