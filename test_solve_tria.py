import ctypes
import numpy as np
from scipy.linalg import lu
import matplotlib.pyplot as plt
import c_solve_tria
import pytest
        
def lu_setup(n_r, n_z, seed, r_min, r_max, z_min, z_max):

    r_vec = np.linspace(r_min, r_max, n_r)
    z_vec = np.linspace(z_min, z_max, n_z)
    
    dr = r_vec[1] - r_vec[0]
    dz = z_vec[1] - z_vec[0]
    
    n_grid = n_r*n_z
    sys_mat = make_finite_diff_matrix(n_r, n_z, n_grid, r_vec, dr, dz)
    perm, lower, upper = lu(sys_mat)
    
    lower_pad = np.pad(lower, ((0, 0), (n_r, 0)))
    upper_pad = np.pad(upper, ((0, 0), (0, n_r+1)))

    lower_band = np.zeros((n_grid, n_r))
    upper_band = np.zeros((n_grid, n_r+2))

    for ii in range(n_grid):
        lower_band[ii, :] = lower_pad[ii, ii:ii+n_r]
        upper_band[ii, :] = upper_pad[ii, ii:ii+n_r+2]
        
    # replace with recipricals 
    upper_band[:, 0] = 1.0/upper_band[:, 0]    
    
    rr, cc = np.where(perm!=0)
    assert(np.all(rr==np.arange(rr.size)))
    
    b_vec = source_vec(n_r, n_z, seed)
    
    output = {
        'sys_mat': sys_mat,
        'lower': lower,
        'upper': upper,
        'perm': perm,
        'lower_band': lower_band,
        'upper_band': upper_band,
        'idx_final': cc,
        'n_row': n_grid,
        'n_rep': n_r,
        'n_r': n_r,
        'n_z': n_z,
        'b_vec': b_vec
        }
        
    return output
        
def source_vec(n_r, n_z, seed):

    n_grid = n_r*n_z
    Rng = np.random.default_rng(seed=seed)
    b_vec = Rng.random(n_grid, )
    
    return b_vec
  

def make_finite_diff_matrix(n_r, n_z, n_grid, r_vec, dr, dz):
        
    a_bar = (dz/dr)**2 * (r_vec / (r_vec+ 0.5*dr))
    a_bar = np.tile(a_bar, (n_z, 1)).flatten()

    c_bar = (dz/dr)**2 * (r_vec / (r_vec - 0.5*dr))   
    c_bar = np.tile(c_bar, (n_z, 1)).flatten()

    b_bar = a_bar + c_bar

    sys_mat = np.diag(-b_bar - 2.0).flatten()

    idx = np.ravel_multi_index((np.arange(n_r, n_grid), 
            np.arange(0, n_grid - n_r)), (n_grid, n_grid)) 
    sys_mat[idx] = 1.0

    idx = np.ravel_multi_index((np.arange(0, n_grid - n_r), 
            np.arange(n_r, n_grid)), (n_grid, n_grid))
    sys_mat[idx] = 1.0

    idx = np.ravel_multi_index((np.arange(0, n_grid-1), 
            np.arange(1, n_grid)), (n_grid, n_grid))
    sys_mat[idx] = a_bar[0:(n_grid-1)]

    idx = np.ravel_multi_index((np.arange(1, n_grid), 
            np.arange(0, n_grid-1)), (n_grid, n_grid))
    sys_mat[idx] = c_bar[1:n_grid] 

    bndry_ind = np.concatenate((np.arange(1, n_r-1), np.arange(0, 
            n_grid, n_r), np.arange(n_r-1, n_grid, 
            n_r), np.arange((n_z-1)*n_r + 1, n_grid -1)))
    
    sys_mat = np.reshape(sys_mat, (n_grid, n_grid))

    sys_mat[bndry_ind, :] = 0.0

    sys_mat = sys_mat.flatten()
    
    idx = np.ravel_multi_index((bndry_ind, bndry_ind), (n_grid, n_grid))

    sys_mat[idx] = 1.0
    
    sys_mat = np.reshape(sys_mat, (n_grid, n_grid))
            
    return sys_mat
        


@pytest.mark.parametrize("n_r", 2**np.arange(3, 6))
@pytest.mark.parametrize("n_z", 2**np.arange(3, 7))
@pytest.mark.parametrize("seed", np.arange(0, 2))            
def test_solve_tria(n_r, n_z, seed, r_min=0.2, r_max=1.0, z_min=-1.0, z_max=1.0, thresh=1.0e-10, show=False):
    
    lu_mtrx = lu_setup(n_r, n_z, seed, r_min=0.2, r_max=1.0, z_min=-1.0, z_max=1.0)
    
    truth = np.linalg.solve(lu_mtrx['lower'], lu_mtrx['perm'] @ lu_mtrx['b_vec']) 
    truth = np.linalg.solve(lu_mtrx['upper'], truth)
    
    c_lower_band = lu_mtrx['lower_band'].ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    c_upper_band = lu_mtrx['upper_band'].ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    c_b_vec =  lu_mtrx['b_vec'].ctypes.data_as(ctypes.POINTER(ctypes.c_double))    
    c_n_row = np.int32(lu_mtrx['n_row'])
    c_n_rep = np.int32(lu_mtrx['n_rep'])
    idx_final = np.asarray(lu_mtrx['idx_final'], dtype=np.int32)
    c_idx_final = idx_final.ctypes.data_as(ctypes.POINTER(ctypes.c_int))   
    
    estimate = np.ascontiguousarray(np.zeros((lu_mtrx['n_row'])))    
    c_estimate = estimate.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    
    c_solve_tria.solve_tria(c_n_row, c_n_rep, c_lower_band, c_upper_band, 
        c_b_vec, c_idx_final, c_estimate)
    estimate = np.ctypeslib.as_array(c_estimate, shape=(lu_mtrx['n_z'], lu_mtrx['n_r']))     

    if show:
        truth = np.reshape(truth, (lu_mtrx['n_z'], lu_mtrx['n_r']))
        estimate = np.reshape(estimate, (lu_mtrx['n_z'], lu_mtrx['n_r']))   
         
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
   
    return lu_mtrx    

            
