import ctypes
import numpy as np
import matplotlib.pyplot as plt
import pytest
from pytest_cases import fixture
from py_to_ctypes import run_c_func
import c_cblas


def test_cblas_dgemm(n_row_a=5, n_row_b=3, n_col_b=7):
    
    a_mat = np.random.random_sample((n_row_a, n_row_b))

    b_mat = np.random.random_sample((n_row_b, n_col_b))

    c_mat = np.zeros((n_row_a, n_col_b))

    # https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-0/cblas-gemm-001.html
    aa, bb, cc, n_row_a, n_col_b, n_row_b, dd, a_mat, n_row_a, b_mat, n_row_b, \
                ee, c_mat, n_row_a = run_c_func( 
                c_cblas.cblas_dgemm, c_cblas.CblasRowMajor, c_cblas.CblasNoTrans, 
                c_cblas.CblasNoTrans, n_row_a, n_col_b, n_row_b, 1.0, a_mat, 
                n_row_b, b_mat, n_col_b, 0.0, c_mat, n_col_b)    
                
    c_mat_gt = a_mat @ b_mat
    
    assert np.all(np.isclose(c_mat_gt, c_mat))
    
    
def test_cblas_dgemv(n_meas=70, n_coef=20):

    arr = np.random.random_sample((n_coef, n_meas))
    x_vec = np.random.random_sample((n_coef,))
    
    b_vec = np.zeros((n_meas,))
    
    b_vec_gt = arr.T @ x_vec

    aa, bb, n_meas, n_coef, cc, arr, n_meas, x_vec, dd, ee, b_vec, ff = run_c_func(c_cblas.cblas_dgemv, 
        c_cblas.CblasColMajor, c_cblas.CblasNoTrans, n_meas, n_coef, 1.0, arr, 
        n_meas, x_vec, 1, 0.0, b_vec, 1);  
        
    assert np.all(np.isclose(b_vec_gt, b_vec))
    
    

if __name__ == "__main__":
    test_cblas_dgemm()
    test_cblas_dgemv()
