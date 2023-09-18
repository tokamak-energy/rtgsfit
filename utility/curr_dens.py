import numpy as np
    

def curr_ellipse(r_vec, z_vec, r0, z0, sigma_r, sigma_z, a_minor):
    
    r_grid, z_grid = np.meshgrid(r_vec, z_vec)
    
    curr_dens = ((r_grid - r0)**2)/(sigma_r**2) + ((z_grid - z0)**2)/(sigma_z**2)
    
    curr_bound = ((a_minor - r0)**2)/(sigma_r**2)
    
    curr_dens = curr_bound - curr_dens 
    curr_dens[curr_dens<0.0] = 0.0
    
    return curr_dens
    
    
def curr_rand(n_r, n_z, seed):

    Rng = np.random.default_rng(seed=seed)
    curr_dens = Rng.random((n_z, n_r), )
    
    return curr_dens
