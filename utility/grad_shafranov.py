import numpy as np
from scipy.linalg import lu
   
def gs_finite_diff(r_vec, z_vec):

    n_r = r_vec.size
    n_z = z_vec.size
    n_grid = n_r*n_z
    dr = r_vec[1] - r_vec[0]
    dz = z_vec[1] - z_vec[0]
        
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
        



def plu(sys_mat):

    perm, lower, upper = lu(sys_mat)
    
    perm = np.linalg.inv(perm)
    
    assert np.all(np.isclose(perm @ sys_mat, lower @ upper))
        
    return perm, lower, upper


def plu_band(perm, lower, upper, n_rep):
   
    n_grid = perm.shape[0]
   
    lower_pad = np.pad(lower, ((0, 0), (n_rep, 0)))
    upper_pad = np.pad(upper, ((0, 0), (0, n_rep+1)))

    lower_band = np.zeros((n_grid, n_rep))
    upper_band = np.zeros((n_grid, n_rep+2))

    for ii in range(n_grid):
        lower_band[ii, :] = lower_pad[ii, ii:ii+n_rep]
        upper_band[ii, :] = upper_pad[ii, ii:ii+n_rep+2]
        
    # replace with recipricals 
    upper_band[:, 0] = 1.0/upper_band[:, 0]    
    
    rr, final_idx = np.where(perm!=0)
    assert(np.all(rr==np.arange(rr.size)))
        
    return final_idx, lower_band, upper_band
