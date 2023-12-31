import ctypes
import numpy as np
import matplotlib.pyplot as plt
import c_libfind_x_point
import pytest
from pytest_cases import fixture
from scipy.constants import mu_0
import sys
sys.path.append("../utility")
from py_to_ctypes import run_c_func
import contourpy
from findiff import FinDiff
from shapely.geometry import Point, Polygon
from scipy.io import loadmat

filename = '../data/12001000_RUN01_for_c.mat'

@fixture(scope="module")
def data():
    return loadmat(filename, squeeze_me=True)

def solovev(r_vec, z_vec, r0=0.4, c0=-1.0, c1=15.0, c2=1.5):

    p_prime = -c1/mu_0;
    ff_prime = -c2 * r0**2;
    [r_grid, z_grid] = np.meshgrid(r_vec, z_vec);

    j_phi = (r_grid * p_prime) + (ff_prime / (mu_0 * r_grid));

    psi_func = lambda r, z: 1 - 0.5*(c2*r0**2 + c0*r**2)*z**2 - 0.125*(c1 - c0) * \
            (r**2 - r0**2)**2

    psi_grid = psi_func(r_grid, z_grid)
            
    z_opt = np.array([0.0])
    r_opt = np.array([r0])
    psi_opt = np.array([psi_func(r_opt, z_opt)])
    
    r_xpt = np.array([1, 1])*np.sqrt(-c2/c0)*r0
    z_xpt = np.array([1, -1])*r0*np.sqrt((c1 - c0)*(c0 + c2)/(2*c0**2))
    psi_xpt = psi_func(r_xpt, z_xpt)  
              
    return psi_grid, r_xpt, z_xpt, psi_xpt, r_opt, z_opt, psi_opt
     
        
def test_find_null_in_gradient(data, thresh=1.e-3):

    r_vec = data['r_vec']
    z_vec = data['z_vec']

    psi, xpt_r_gt, xpt_z_gt, xpt_psi_gt, opt_r_gt, opt_z_gt, opt_psi_gt = solovev(r_vec, z_vec)
    
    opt_r = np.zeros(10)
    opt_z = np.zeros(10)
    opt_psi = np.zeros(10)
    i_opt = np.array([0], dtype=np.int64)
    
    xpt_r = np.zeros(10)
    xpt_z = np.zeros(10)
    xpt_psi = np.zeros(10)
    i_xpt =  np.array([0], dtype=np.int64)
    
    psi, opt_r, opt_z, opt_psi, i_opt, xpt_r, xpt_z, xpt_psi, i_xpt \
            = run_c_func(c_libfind_x_point.find_null_in_gradient, psi, opt_r, 
            opt_z, opt_psi, i_opt, xpt_r, xpt_z, xpt_psi, i_xpt)
        
    i_xpt = i_xpt[0]
    i_opt = i_opt[0]
    
    for ii in range(i_xpt):
        assert np.any(np.isclose(xpt_r[ii], xpt_r_gt, rtol=thresh, atol=np.inf))         
        assert np.any(np.isclose(xpt_z[ii], xpt_z_gt, rtol=thresh, atol=np.inf))         

    for rr, zz in zip(xpt_r_gt, xpt_z_gt):
        assert np.any(np.isclose(rr, xpt_r[:i_xpt], rtol=thresh, atol=np.inf))         
        assert np.any(np.isclose(zz, xpt_z[:i_xpt], rtol=thresh, atol=np.inf))       
        
    for ii in range(i_opt):
        assert np.any(np.isclose(opt_r[ii], opt_r_gt, rtol=thresh, atol=np.inf))         
        assert np.any(np.isclose(opt_z[ii], opt_z_gt, rtol=thresh, atol=np.inf))         
        
    for rr, zz in zip(opt_r_gt, opt_z_gt):
        assert np.any(np.isclose(rr, opt_r[:i_opt], rtol=thresh, atol=np.inf))         
        assert np.any(np.isclose(zz, opt_z[:i_opt], rtol=thresh, atol=np.inf))   


def test_find_null_in_gradient_march(data, thresh=1.e-3):

    r_vec = data['r_vec']
    z_vec = data['z_vec']

    psi, xpt_r_gt, xpt_z_gt, xpt_psi_gt, opt_r_gt, opt_z_gt, opt_psi_gt = solovev(r_vec, z_vec)
    
    opt_r = np.zeros(10)
    opt_z = np.zeros(10)
    opt_psi = np.zeros(10)
    i_opt = np.array([0], dtype=np.int64)
    
    xpt_r = np.zeros(10)
    xpt_z = np.zeros(10)
    xpt_psi = np.zeros(10)
    i_xpt =  np.array([0], dtype=np.int64)
    
    psi, opt_r, opt_z, opt_psi, i_opt, xpt_r, xpt_z, xpt_psi, i_xpt \
            = run_c_func(c_libfind_x_point.find_null_in_gradient_march, psi, opt_r, 
            opt_z, opt_psi, i_opt, xpt_r, xpt_z, xpt_psi, i_xpt)
        
    i_xpt = i_xpt[0]
    i_opt = i_opt[0]
    
    for ii in range(i_xpt):
        assert np.any(np.isclose(xpt_r[ii], xpt_r_gt, rtol=thresh, atol=np.inf))         
        assert np.any(np.isclose(xpt_z[ii], xpt_z_gt, rtol=thresh, atol=np.inf))         

    for rr, zz in zip(xpt_r_gt, xpt_z_gt):
        assert np.any(np.isclose(rr, xpt_r[:i_xpt], rtol=thresh, atol=np.inf))         
        assert np.any(np.isclose(zz, xpt_z[:i_xpt], rtol=thresh, atol=np.inf))       
        
    for ii in range(i_opt):
        assert np.any(np.isclose(opt_r[ii], opt_r_gt, rtol=thresh, atol=np.inf))         
        assert np.any(np.isclose(opt_z[ii], opt_z_gt, rtol=thresh, atol=np.inf))         
        
    for rr, zz in zip(opt_r_gt, opt_z_gt):
        assert np.any(np.isclose(rr, opt_r[:i_opt], rtol=thresh, atol=np.inf))         
        assert np.any(np.isclose(zz, opt_z[:i_opt], rtol=thresh, atol=np.inf))   
        
    plt.contour(r_vec, z_vec, psi, 50)
    plt.plot(opt_r[:i_opt], opt_z[:i_opt], '+')
    plt.plot(xpt_r[:i_xpt], xpt_z[:i_xpt], 'x')    
    plt.show()
    
    
def test_find_lcfs(data):

    r_vec = data['r_vec']
    z_vec = data['z_vec']
    
    psi_gt, xpt_r_gt, xpt_z_gt, xpt_psi_gt, opt_r_gt, opt_z_gt, opt_psi_gt = solovev(r_vec, z_vec)
    
    ContGen = contourpy.contour_generator(r_vec, z_vec, psi_gt)
    
    n_row = z_vec.size
    n_col = r_vec.size
    
    dr = r_vec[1] - r_vec[0]
    dz = z_vec[1] - z_vec[0]
    
    psi_bound = np.max(xpt_psi_gt)*1.001
    intrscts = np.concatenate(ContGen.lines(psi_bound))
   
    r_lcfs = np.zeros(500)
    z_lcfs = np.zeros(500)
    n_lcfs = int(0)
    psi = psi_gt.copy()
    
    psi, psi_bound, r_lcfs, z_lcfs, n_lcfs = run_c_func(c_libfind_x_point.find_lcfs_rz, 
            psi, psi_bound, r_lcfs, z_lcfs, n_lcfs)

    rz_lcfs = np.vstack((r_lcfs[:n_lcfs], z_lcfs[:n_lcfs])).T

    for rz in zip(rz_lcfs):
        assert np.any(np.all(np.isclose(rz, intrscts), axis=1))
        
    for rz in zip(intrscts):
        assert np.any(np.all(np.isclose(rz, rz_lcfs), axis=1))  
        
          
def test_inside_lcfs(data):

    r_vec = data['r_vec']
    z_vec = data['z_vec']
    
    psi_gt, xpt_r_gt, xpt_z_gt, xpt_psi_gt, opt_r_gt, opt_z_gt, opt_psi_gt = solovev(r_vec, z_vec)
    
    ContGen = contourpy.contour_generator(r_vec, z_vec, psi_gt)
    
    n_row = z_vec.size
    n_col = r_vec.size
    
    dr = r_vec[1] - r_vec[0]
    dz = z_vec[1] - z_vec[0]
    
    psi_bound = np.max(xpt_psi_gt)*1.001
    intrscts = ContGen.lines(psi_bound)
    
    is_closed = np.zeros(len(intrscts), dtype=bool)
    count = 0
    for ii in range(len(intrscts)):
        if np.all(np.isclose(intrscts[ii][0, :], intrscts[ii][-1, :])):
            idx = ii
            count = count + 1
    assert count == 1

    intrscts = intrscts[idx]
    
    r_opt = opt_r_gt[0]
    z_opt = opt_z_gt[0]    
    psi_bound = np.max(xpt_psi_gt)*1.001
    
    r_lcfs = intrscts[:, 0]
    z_lcfs = intrscts[:, 1]
    n_lcfs = intrscts.shape[0]
    mask = np.zeros((n_row, n_col), dtype=np.int64)
        
    r_opt, z_opt, r_lcfs, z_lcfs, n_lcfs, mask = run_c_func(c_libfind_x_point.inside_lcfs, r_opt, z_opt, 
            r_lcfs, z_lcfs, n_lcfs, mask)  
    
        
    mask = mask.astype(bool)
    
    poly = Polygon(intrscts)
    mask_gt = np.zeros((n_row, n_col), dtype=bool)
    
    for i_row in range(n_row):
        for i_col in range(n_col):
            mask_gt[i_row, i_col] = poly.contains(Point(r_vec[i_col], z_vec[i_row]))


    assert np.all(mask_gt == mask)


#def test_flux_to_mask(data):

#    r_vec = data['r_vec']
#    z_vec = data['z_vec']
#    
#    psi_gt, xpt_r_gt, xpt_z_gt, xpt_psi_gt, opt_r_gt, opt_z_gt, opt_psi_gt = solovev(r_vec, z_vec)
#    
#    psi = psi_gt.copy()
#    n_row = z_vec.size
#    n_col = r_vec.size
#    
#    dr = r_vec[1] - r_vec[0]
#    dz = z_vec[1] - z_vec[0]
#    
#    hess_row_col = FinDiff((0, dz, 1), (1, dr, 1), acc=2)
#    hess_row_row = FinDiff((0, dz, 2), acc=2)
#    hess_col_col = FinDiff((1, dr, 2), acc=2)
#    grad_row = FinDiff((0, dz, 1), acc=2)
#    grad_col = FinDiff((1, dr, 1), acc=2)
#    
#    grad_r = grad_col(psi)
#    grad_z = grad_row(psi)
#    hess_rr = hess_col_col(psi)
#    hess_zz = hess_row_row(psi)
#    hess_rz = hess_row_col(psi)
#    
#    opt_r = np.zeros(10)
#    opt_z = np.zeros(10)
#    opt_psi = np.zeros(10)
#    i_opt = int(0)
#    
#    xpt_r = np.zeros(10)
#    xpt_z = np.zeros(10)
#    xpt_psi = np.zeros(10)
#    i_xpt = int(0)
#    
#    psi, opt_r, opt_z, opt_psi, i_opt, xpt_r, xpt_z, xpt_psi, i_xpt \
#        = run_c_func(c_libfind_x_point.find_null_in_gradient, psi, opt_r, 
#        opt_z, opt_psi, i_opt, xpt_r, xpt_z, xpt_psi, i_xpt)
#    
#    r_lcfs = np.zeros(n_row*n_col)
#    z_lcfs = np.zeros(n_row*n_col)
#    n_lcfs = int(0)

#    psi_bound = np.max(xpt_psi)*1.001

#    psi, psi_bound, r_lcfs, z_lcfs, n_lcfs = run_c_func(c_libfind_x_point.find_lcfs_rz, 
#            psi, psi_bound, r_lcfs, z_lcfs, n_lcfs) 
#            
#    r_opt = opt_r[0]
#    z_opt = opt_z[0]    
#    idx = np.zeros(n_row*n_col, dtype=np.int64)
#    n_idx = int(0)
#                
#   
#    mask = np.zeros((n_row, n_col), dtype=np.int64)
#        
#    r_opt, z_opt, r_lcfs, z_lcfs, n_lcfs, mask = run_c_func(c_libfind_x_point.inside_lcfs, r_opt, z_opt, 
#            r_lcfs, z_lcfs, n_lcfs, mask)  
#    
#    mask = mask.astype(bool)
#    
#    ContGen = contourpy.contour_generator(r_vec, z_vec, psi_gt)
#    intrscts = ContGen.lines(np.max(xpt_psi_gt)*1.001)
#    
#    is_closed = np.zeros(len(intrscts), dtype=bool)
#    count = 0
#    for ii in range(len(intrscts)):
#        if np.all(np.isclose(intrscts[ii][0, :], intrscts[ii][-1, :])):
#            idx = ii
#            count = count + 1

#    assert count == 1

#    intrscts = intrscts[idx]
#    
#    poly = Polygon(intrscts)
#    mask_gt = np.zeros((n_row, n_col), dtype=bool)
#    
#    for i_row in range(n_row):
#        for i_col in range(n_col):
#            mask_gt[i_row, i_col] = poly.contains(Point(r_vec[i_col], z_vec[i_row]))
#    
#    breakpoint()
#    assert np.all(mask_gt == mask)
               
if __name__ == "__main__":

    data = loadmat(filename, squeeze_me=True)
    test_find_null_in_gradient_march(data)
#    test_flux_to_mask(r_vec, z_vec)
#    psi, r_xpt, z_xpt, psi_xpt, r_opt, z_opt, psi_opt = solovev(r_vec, z_vec)
#    
#    plt.contour(r_vec, z_vec, psi, 50)
#    plt.contour(r_vec, z_vec, psi, [psi_xpt[0]], colors='black')
#    plt.plot(r_xpt, z_xpt, 'x')
#    plt.plot(r_opt, z_opt, '+')
#    plt.show()
