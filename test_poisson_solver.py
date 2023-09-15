import ctypes
import numpy as np
import matplotlib.pyplot as plt
import c_poisson_solver
import pytest
from scipy.constants import mu_0
from findiff import FinDiff
from curr_dens import curr_ellipse
from py_to_ctypes import run_c_func
from pytest_cases import fixture
from scipy.io import loadmat

    
def hagenow(lower, upper, perm, source, inv_r_mu0, g_bound, n_z, n_r, dz, dr):

    psi_zero = np.linalg.solve(lower, perm @ source) 
    psi_zero = np.linalg.solve(upper, psi_zero)

    psi_zero = np.reshape(psi_zero, (n_z, n_r))

    grad_row = FinDiff((0, dz, 1), acc=2)
    truth_row = grad_row(psi_zero)

    grad_col = FinDiff((1, dr, 1), acc=2)
    truth_col = grad_col(psi_zero)        

    dpsi_bound = np.concatenate([truth_col[:, 0]*dz, -truth_row[-1, 1:-1]*dr, 
            -truth_col[::-1, -1]*dz, truth_row[0, -2:0:-1]*dr])
            
    dpsi_hagen = inv_r_mu0 * dpsi_bound

    psi_bound = g_bound @ dpsi_hagen
    
    return psi_bound

#n_r_list, n_z_list = np.meshgrid([16, 32, 64], [16,32, 64])

#r_min=0.11
#r_max=1.07
#z_min=-0.96
#z_max=0.96
r0 = 0.5
z0 = 0.0
a_minor = 0.2
sigma_r = 1.0
sigma_z = 2.0
thresh=1.0e-10
show=True
filename = '12001000_RUN01_for_python.mat'

#@fixture(scope="module", params=[16, 32, 64])
#def r_vec(request):
#    return np.linspace(r_min, r_max, request.param)
#    
#    
#@fixture(scope="module", params=[16, 32, 64])
#def z_vec(request):
#    return np.linspace(z_min, z_max, request.param)  
#     

#@fixture(scope="module", unpack_into="perm,lower,upper,perm_idx,lower_band,upper_band")
#def gs_mtrx(r_vec, z_vec):

#    sys_mat= gs_finite_diff(r_vec, z_vec)
#    perm, lower, upper = plu(sys_mat)
#    perm_idx, lower_band, upper_band = plu_band(perm, lower, upper, r_vec.size)         
#    return perm, lower, upper, perm_idx, lower_band, upper_band
   


    

#@fixture(scope="module", unpack_into="g_grid,g_bound,inv_r_mu0")
#def greens(r_vec, z_vec): 

#    n_r, n_z, _, dr, dz, _ = grid_params(r_vec, z_vec)
#    r_grid, z_grid = np.meshgrid(r_vec, z_vec)     
#    r_ltrb = np.concatenate((r_grid[:, 0], r_grid[0, :], r_grid[:, -1], r_grid[-1, :]))
#    z_ltrb = np.concatenate((z_grid[:, 0], z_grid[0, :], z_grid[:, -1], z_grid[-1, :]))    
#    g_grid = greens_grid(r_grid, z_grid, dr, dz)
#    g_bound = greens_bound(r_ltrb, z_ltrb, dr, dz, n_r, n_z)   
#    inv_r_mu0 = 1.0/(mu_0 * r_ltrb)    
#    return g_grid, g_bound, inv_r_mu0
#    

#def grid_params(r_vec, z_vec):
#    n_r = r_vec.size
#    n_z = z_vec.size
#    n_ele = n_r*n_z
#    dz = z_vec[1] - z_vec[0]
#    dr = r_vec[1] - r_vec[0]    
#    n_bound = 2*(n_r + n_z)
#    
#    return n_r, n_z, n_ele, dr, dz, n_bound

@fixture(scope="module")
def data():
    return loadmat(filename, squeeze_me=True)

@fixture(scope="module")
def curr_dens(data):    
    return  curr_ellipse(data['r_vec'], data['z_vec'], r0, z0, sigma_r, sigma_z, a_minor)

def curr_to_source(curr, r_vec, dz):
    source = -curr*mu_0*r_vec.reshape((1, -1))*(dz**2)
    return source.flatten()


def test_hagenow_py(curr_dens, data):
    
    n_r = data['n_r']
    n_z = data['n_z']
    n_grid = data['n_grid']
    dr = data['dr']
    dz = data['dz']
    r_vec = data['r_vec']
    lower = data['lower']
    upper = data['upper']
    perm = data['perm']
    inv_r_mu0 = data['inv_r_ltrb_mu0']
    g_bound = data['g_ltrb']
    n_ltrb = data['n_ltrb']
    
    source = curr_to_source(curr_dens, r_vec, dz)
          
    py_psi_ltrb = hagenow(lower, upper, perm, 
            source, inv_r_mu0, g_bound, n_z, 
            n_r, dz, dr)

    psi = np.zeros(n_grid, )
    c_psi_ltrb = np.zeros(n_ltrb, )

    source, psi, c_psi_ltrb, = run_c_func( \
            c_poisson_solver.hagenow_bound, source, psi, c_psi_ltrb)
            
    if show:
        fig, ax = plt.subplots(1, 3)
        ax[0].imshow(np.reshape(curr_dens,(n_z, n_r)))
        ax[1].imshow(np.reshape(psi, (n_z, n_r)))
        ax[2].plot(py_psi_ltrb, 'x')
        ax[2].plot(c_psi_ltrb, '+')
        ax[2].legend(['py-hag', 'c-hag'])
        plt.show()        
    
    for ii, (est, py_est) in enumerate(zip(c_psi_ltrb, py_psi_ltrb)):
        assert np.abs((est-py_est)/py_est) < thresh, \
                f"{ii} {est} {py_est} {np.abs((est-py_est)/py_est)} {thresh}"

          
    
def test_hagenow_true(curr_dens, data):    

    n_r = data['n_r']
    n_z = data['n_z']
    n_grid = data['n_grid']
    dr = data['dr']
    dz = data['dz']
    r_vec = data['r_vec']
    n_ltrb = data['n_ltrb']
    lower = data['lower']
    upper = data['upper']
    perm = data['perm']
    inv_r_mu0 = data['inv_r_ltrb_mu0']
    g_bound = data['g_ltrb']
    g_grid = data['g_grid_grid']
    
    source = curr_to_source(curr_dens, r_vec, dz)
    
    truth_psi = g_grid @ curr_dens.flatten() * dr * dz
    truth_psi = np.reshape(truth_psi, (n_z, n_r))
    truth_psi_ltrb = np.concatenate((truth_psi[:, 0], truth_psi[-1, 1:-1], 
            truth_psi[::-1, -1], truth_psi[0, -2:0:-1]))
        
    psi = np.zeros(n_grid, )
    c_psi_ltrb = np.zeros(n_ltrb, )

    source, psi, c_psi_ltrb = run_c_func( \
            c_poisson_solver.hagenow_bound, source, psi, c_psi_ltrb)
            
    if show:
        fig, ax = plt.subplots(1, 3)
        ax[0].imshow(curr_dens)
        ax[1].imshow(np.reshape(psi,(n_z, n_r)))
        ax[2].plot(truth_psi_ltrb, 'x')
        ax[2].plot(c_psi_ltrb, '+')
        ax[2].legend(['truth', 'c-hag'])
        plt.show()        
    
    for ii, (est, tru) in enumerate(zip(c_psi_ltrb, truth_psi_ltrb)):
        assert np.abs((est - tru)/tru) < (2.0/n_z), \
                f"{ii} {est} {tru} {np.abs((est - tru)/tru)} {(1.0/np.sqrt(n_r*n_z))}"
          

def test_poisson_true(curr_dens, data):    

    n_r = data['n_r']
    n_z = data['n_z']
    n_grid = data['n_grid']
    dr = data['dr']
    dz = data['dz']
    r_vec = data['r_vec']
    n_ltrb = data['n_ltrb']
    lower = data['lower']
    upper = data['upper']
    perm = data['perm']
    inv_r_mu0 = data['inv_r_ltrb_mu0']
    g_bound = data['g_ltrb']
    g_grid = data['g_grid_grid']

    source =curr_to_source(curr_dens, r_vec, dz)
    
    truth_psi = g_grid @ curr_dens.flatten() * dr * dz
    truth_psi = np.reshape(truth_psi, (n_z, n_r))

    psi = np.zeros(n_grid, )

    source, psi = run_c_func(c_poisson_solver.poisson_solver, source, psi) 

    if show:
        fig, ax = plt.subplots(1, 3)
        ax[0].imshow(np.reshape(truth_psi,(n_z, n_r)))
        ax[1].imshow(np.reshape(psi,(n_z, n_r)))
        ax[2].plot(truth_psi.flatten(), 'x')
        ax[2].plot(psi.flatten(), '+')
        ax[2].legend(['truth', 'c-hag'])
        plt.show()        

    for ii, (est, tru) in enumerate(zip(psi.flatten(), truth_psi.flatten())):
        assert np.abs((est - tru)/tru) < (2.0/n_z), \
                f"{ii} {est} {tru} {np.abs((est - tru)/tru)} {(1.0/np.sqrt(n_r*n_z))}"
    
        

            
