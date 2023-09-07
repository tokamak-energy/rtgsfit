import ctypes
import numpy as np
import os 
import inspect
import matplotlib.pyplot as plt
from findiff import FinDiff
import c_gradient
import pytest


class TestGradient:
    
#    @classmethod
    def setup_class(
            self, 
            func = lambda x, y: np.sin(x*np.pi)*np.sin(y*np.pi), 
            n_row = 11, 
            n_col = 11,
            thresh = 1e-10,
            show = False
            ):
        
        self.n_row = n_row
        self.d_row = 1. / (n_row - 1)
        self.n_col = n_col
        self.d_col = 1. / (n_col - 1)
        self.x = np.linspace(0., 1., n_col)
        self.y = np.linspace(0., 1., n_row)
        self.xx = np.meshgrid(self.x, self.y)[0]
        self.yy = np.meshgrid(self.x, self.y)[1]
        self.data = np.ascontiguousarray(func(self.xx, self.yy))
        self.p_data = self.data.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        self.show = show
        self.thresh = thresh
    
        
    def test_gradient_row(self):
    
        grad_row = FinDiff((0, self.d_row, 1), acc=2)
        truth = grad_row(self.data)
        
        estimate = np.ascontiguousarray(np.zeros((self.n_row, self.n_col)))
        p_estimate = estimate.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        
        c_gradient.gradient_row(self.p_data, np.int32(self.n_row), 
                np.int32(self.n_col), self.d_row, p_estimate)
                
        estimate = np.ctypeslib.as_array(p_estimate, shape=(self.n_row, self.n_col))    
            
        for est, tru in zip(estimate.flatten(), truth.flatten()):
            assert(np.abs(est-tru) < self.thresh)
        
        self.plot(estimate, truth)
    
        
    def test_gradient_col(self):
    
        grad_col = FinDiff((1, self.d_col, 1), acc=2)
        truth = grad_col(self.data)
        
        estimate = np.ascontiguousarray(np.zeros((self.n_row, self.n_col)))
        p_estimate = estimate.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        
        c_gradient.gradient_col(self.p_data, np.int32(self.n_row), 
                np.int32(self.n_col), self.d_col, p_estimate)
                
        estimate = np.ctypeslib.as_array(p_estimate, shape=(self.n_row, self.n_col))    
            
        for est, tru in zip(estimate.flatten(), truth.flatten()):
            assert(np.abs(est-tru) < self.thresh)
        
        self.plot(estimate, truth)        
    
        
    def test_hessian_row_row(self):
    
        hess_row_row = FinDiff((0, self.d_row, 2), acc=2)
        truth = hess_row_row(self.data)
        
        estimate = np.ascontiguousarray(np.zeros((self.n_row, self.n_col)))
        p_estimate = estimate.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        
        c_gradient.hessian_row_row(self.p_data, np.int32(self.n_row), 
                np.int32(self.n_col), self.d_row, p_estimate)
                
        estimate = np.ctypeslib.as_array(p_estimate, shape=(self.n_row, self.n_col))    
            
        for est, tru in zip(estimate.flatten(), truth.flatten()):
            assert(np.abs(est-tru) < self.thresh)
        
        self.plot(estimate, truth)   


    def test_hessian_col_col(self):
    
        hess_col_col= FinDiff((1, self.d_row, 2), acc=2)
        truth = hess_col_col(self.data)
        
        estimate = np.ascontiguousarray(np.zeros((self.n_row, self.n_col)))
        p_estimate = estimate.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        
        c_gradient.hessian_col_col(self.p_data, np.int32(self.n_row), 
                np.int32(self.n_col), self.d_col, p_estimate)
                
        estimate = np.ctypeslib.as_array(p_estimate, shape=(self.n_row, self.n_col))    
            
        for est, tru in zip(estimate.flatten(), truth.flatten()):
            assert(np.abs(est-tru) < self.thresh)
        
        self.plot(estimate, truth)   
    
        
    def test_hessian_row_col(self):
    
        hess_row_col = FinDiff((0, self.d_row, 1), (1, self.d_col, 1), acc=2)
        truth = hess_row_col(self.data)
        
        estimate = np.ascontiguousarray(np.zeros((self.n_row, self.n_col)))
        p_estimate = estimate.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

        estimate_tmp = np.ascontiguousarray(np.zeros((self.n_row, self.n_col)))
        p_estimate_tmp = estimate_tmp.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                
        c_gradient.gradient_row(self.p_data, np.int32(self.n_row), 
                np.int32(self.n_col), self.d_row, p_estimate_tmp)

        c_gradient.gradient_col(p_estimate_tmp, np.int32(self.n_row), 
                np.int32(self.n_col), self.d_col, p_estimate)
                                
        estimate = np.ctypeslib.as_array(p_estimate, shape=(self.n_row, self.n_col))    
            
        for est, tru in zip(estimate.flatten(), truth.flatten()):
            assert(np.abs(est-tru) < self.thresh)
        
        self.plot(estimate, truth) 
    

    def test_gradient_bound(self):
    
        grad_row = FinDiff((0, self.d_row, 1), acc=2)
        truth_row = grad_row(self.data)
        
        grad_col = FinDiff((1, self.d_col, 1), acc=2)
        truth_col = grad_col(self.data)        
        
        truth = np.concatenate([truth_col[:, 0]*self.d_row, truth_row[0, :]*self.d_col, 
                -truth_col[:, -1]*self.d_row, -truth_row[-1, :]*self.d_col])
        
        estimate = np.ascontiguousarray(np.zeros((self.n_row, self.n_col)))
        p_estimate = estimate.ctypes.data_as(ctypes.POINTER(ctypes.c_double))         
               
        c_gradient.gradient_bound(self.p_data, np.int32(self.n_row), 
                np.int32(self.n_col), self.d_row, self.d_col, p_estimate)
                                
        estimate = np.ctypeslib.as_array(p_estimate, shape=(2*(self.n_row + 
                self.n_col),))    
        
        if self.show:
            plt.plot(estimate)
            plt.plot(truth)
            plt.show()                             
        
        for est, tru in zip(estimate.flatten(), truth.flatten()):
            assert(np.abs(est-tru) < self.thresh)
        
                                
    def plot(self, estimate, truth):
        
        if self.show:
            fig, ax = plt.subplots(1,3)
            ax[0].imshow(truth)
            ax[1].imshow(estimate)
            ax[2].plot(truth, '+')
            ax[2].set_prop_cycle(None)
            ax[2].plot(estimate, 'x')
            plt.show()

 

