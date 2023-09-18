import numpy as np
from scipy.constants import mu_0
from scipy.integrate import quad
from scipy.special import ellipk, ellipe

def greens(r_c, z_c, r, z):

    # Calculate k^2
    k2 = 4.0 * r * r_c / ((r + r_c) ** 2 + (z - z_c) ** 2)

    # Clip to between 0 and 1 to avoid nans e.g. when coil is on grid point
    k2 = np.clip(k2, 1e-10, 1.0 - 1e-10)
    k = np.sqrt(k2)

    return ((mu_0 / (2.0 * np.pi))
        * np.sqrt(r * r_c)
        * ((2.0 - k2) * ellipk(k2) - 2.0 * ellipe(k2))
        / k)


def greens_grid(r_grid, z_grid, dr, dz):
    """
    Calculate poloidal flux at (R,Z) due to a unit current
    at (Rc,Zc) using Greens function

    """
    
    r_c = np.reshape(r_grid, (-1, 1))
    z_c = np.reshape(z_grid, (-1, 1))
    
    r = r_c.T
    z = z_c.T
    
    g_grid = greens(r_c, z_c, r, z)
            
    g_self = mu_0 * r_grid * (np.log(8 * r_grid/(dr + dz)) - 0.5) / (2*np.pi)
    np.fill_diagonal(g_grid, g_self.flatten())
    
    return g_grid
    
    
def greens_bound(r_ltrb, z_ltrb, dr, dz, n_r, n_z):
            
    r_c = np.reshape(r_ltrb, (-1, 1))
    z_c = np.reshape(z_ltrb, (-1, 1))
    
    r = r_c.T
    z = z_c.T
        
    g_grid = greens(r_c, z_c, r, z)
    
    self_induc = np.zeros(2*(n_r+n_z))
    
    f = lambda x: greens(r_ltrb[0], x, r_ltrb[0], z_ltrb[0])
    self_induc[0:n_z] = quad(f, z_ltrb[0]-0.5*dz, z_ltrb[0]+0.5*dz)[0]/dz

    for jj in range(n_z, n_z+n_r):
        f = lambda x: greens(x,  z_ltrb[jj], r_ltrb[jj], z_ltrb[jj])
        self_induc[jj] = quad(f, r_ltrb[jj]-0.5*dr, r_ltrb[jj]+0.5*dr)[0]/dr  

    f = lambda x: greens(r_ltrb[n_r+n_z], x, r_ltrb[n_r+n_z], z_ltrb[n_r+n_z])
    self_induc[n_r+n_z:n_r+n_z*2] = quad(f, z_ltrb[n_r+n_z]-0.5*dz, 
                z_ltrb[n_r+n_z]+0.5*dz)[0]/dz
                                      
    for jj in range(n_z*2 + n_r, 2*n_z+2*n_r):
        f = lambda x: greens(x,  z_ltrb[jj], r_ltrb[jj], z_ltrb[jj])
        self_induc[jj] = quad(f, r_ltrb[jj]-0.5*dr, r_ltrb[jj]+0.5*dr)[0]/dr                                    

    np.fill_diagonal(g_grid, self_induc)
    
    return g_grid
