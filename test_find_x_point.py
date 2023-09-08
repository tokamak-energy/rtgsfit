import ctypes
import numpy as np
import matplotlib.pyplot as plt
import c_find_x_point
import pytest
from pytest_cases import fixture
from scipy.constants import mu_0
from py_to_ctypes import run_c_func
import contourpy

def solovev(r_vec, z_vec, r0=0.5, c0=-1.0, c1=15.0, c2=1.5):

    p_prime = -c1/mu_0;
    ff_prime = -c2 * r0**2;
    [r_grid, z_grid] = np.meshgrid(r_vec, z_vec);

    j_phi = (r_grid * p_prime) + (ff_prime / (mu_0 * r_grid));

    psi_func = lambda r, z:  0.5*(c2*r0**2 + c0*r**2)*z**2 + 0.125*(c1 - c0) * \
            (r**2 - r0**2)**2

    psi_grid = psi_func(r_grid, z_grid)
            
    z_opt = 0.0
    r_opt = r0
    psi_opt = psi_func(r_opt, z_opt)
    
    r_xpt = np.array([1, 1])*np.sqrt(-c2/c0)*r0
    z_xpt = np.array([1, -1])*r0*np.sqrt((c1 - c0)*(c0 + c2)/(2*c0**2))
    psi_xpt = psi_func(r_xpt, z_xpt)  
              
    return psi_grid, r_xpt, z_xpt, psi_xpt, r_opt, z_opt, psi_opt
    
    
def test_find_zero_on_edge(n_test = 1000, thresh=1e-10):

    for ii in range(n_test):

        grad_patch = np.random.random_sample([2, 2])*2 - 1

        ContGen = contourpy.contour_generator((0, 1), (0, 1), grad_patch)
        intrscts = ContGen.lines(0.0)

        if len(intrscts) == 0:
            intrscts_count = 0
        else:
            intrscts = np.concatenate(intrscts, axis=0) 
            intrscts_count = intrscts.shape[0]

        count = np.zeros(1, dtype=np.int64)
        cross_row = np.zeros(4)
        cross_col = np.zeros(4)
        
        grad_patch, thresh, cross_row, cross_col, count = run_c_func( 
                c_find_x_point.find_zero_on_edge, grad_patch, thresh, cross_row, 
                cross_col, count)

        assert count[0] == intrscts_count

        for jj in range(count[0]):
            assert np.any(np.isclose(cross_col[jj],  intrscts[:, 0]))
            assert np.any(np.isclose(cross_row[jj], intrscts[:, 1]))
            
            
    
if __name__ == "__main__":


    r_start = 0.14;
    r_end = 1.0;
    n_r = 33;
    z_start = -1.96; 
    z_end = 1.96;
    n_z = 65; 
    r_vec = np.linspace(r_start, r_end, n_r);
    z_vec = np.linspace(z_start, z_end, n_z);
    
    psi, r_xpt, z_xpt, psi_xpt, r_opt, z_opt, psi_opt = solovev(r_vec, z_vec)
    
    plt.contour(r_vec, z_vec, psi, 50)
    plt.contour(r_vec, z_vec, psi, [psi_xpt[0]], colors='black')
    plt.plot(r_xpt, z_xpt, 'x')
    plt.plot(r_opt, z_opt, '+')
    plt.show()
