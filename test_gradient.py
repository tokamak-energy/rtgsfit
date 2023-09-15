import ctypes
import numpy as np
import os 
import inspect
import matplotlib.pyplot as plt
from findiff import FinDiff
import c_gradient
import pytest
from scipy.io import loadmat


class TestGradient:
    
#    @classmethod
    def setup_class(
            self, 
            func = lambda x, y: np.sin(x*np.pi)*np.sin(y*np.pi), 
            datafile = '12001000_RUN01_for_python.mat',
            thresh = 1e-10,
            show = False
            ):
        
        data = loadmat(datafile, squeeze_me=True)
        
        self.n_z = data['n_z']
        self.dz = data['dz']
        self.n_r = data['n_r']
        self.dr = data['dr']        
        x = np.linspace(0, 1, self.n_r)
        y = np.linspace(0, 1, self.n_z)
        xx, yy = np.meshgrid(x, y)
        self.data = np.ascontiguousarray(func(xx, yy))
        self.p_data = self.data.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        self.show = show
        self.thresh = thresh
    
        
    def test_gradient_z(self):
    
        gradz = FinDiff((0, self.dz, 1), acc=2)
        truth = gradz(self.data)
        
        estimate = np.ascontiguousarray(np.zeros((self.n_z, self.n_r)))
        p_estimate = estimate.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        
        c_gradient.gradient_z(self.p_data, p_estimate)
                
        estimate = np.ctypeslib.as_array(p_estimate, shape=(self.n_z, self.n_r))    
            
        for est, tru in zip(estimate.flatten(), truth.flatten()):
            assert(np.abs(est-tru) < self.thresh)
        
        self.plot(estimate, truth)
    
        
    def test_gradient_r(self):
    
        gradr = FinDiff((1, self.dr, 1), acc=2)
        truth = gradr(self.data)
        
        estimate = np.ascontiguousarray(np.zeros((self.n_z, self.n_r)))
        p_estimate = estimate.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        
        c_gradient.gradient_r(self.p_data,  p_estimate)
                
        estimate = np.ctypeslib.as_array(p_estimate, shape=(self.n_z, self.n_r))    
            
        for est, tru in zip(estimate.flatten(), truth.flatten()):
            assert(np.abs(est-tru) < self.thresh)
        
        self.plot(estimate, truth)        
    
        
    def test_hessian_zz(self):
    
        hess_row_row = FinDiff((0, self.dz, 2), acc=2)
        truth = hess_row_row(self.data)
        
        estimate = np.ascontiguousarray(np.zeros((self.n_z, self.n_r)))
        p_estimate = estimate.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        
        c_gradient.hessian_zz(self.p_data, p_estimate)
                
        estimate = np.ctypeslib.as_array(p_estimate, shape=(self.n_z, self.n_r))    
            
        for est, tru in zip(estimate.flatten(), truth.flatten()):
            assert(np.abs(est-tru) < self.thresh)
        
        self.plot(estimate, truth)   


    def test_hessian_rr(self):
    
        hess_col_col= FinDiff((1, self.dr, 2), acc=2)
        truth = hess_col_col(self.data)
        
        estimate = np.ascontiguousarray(np.zeros((self.n_z, self.n_r)))
        p_estimate = estimate.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        
        c_gradient.hessian_rr(self.p_data, p_estimate)
                
        estimate = np.ctypeslib.as_array(p_estimate, shape=(self.n_z, self.n_r))    
            
        for est, tru in zip(estimate.flatten(), truth.flatten()):
            assert(np.abs(est-tru) < self.thresh)
        
        self.plot(estimate, truth)   
    
        
    def test_hessian_zr(self):
    
        hess_row_col = FinDiff((0, self.dz, 1), (1, self.dr, 1), acc=2)
        truth = hess_row_col(self.data)
        
        estimate = np.ascontiguousarray(np.zeros((self.n_z, self.n_r)))
        p_estimate = estimate.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

        estimate_tmp = np.ascontiguousarray(np.zeros((self.n_z, self.n_r)))
        p_estimate_tmp = estimate_tmp.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                
        c_gradient.gradient_z(self.p_data, p_estimate_tmp)

        c_gradient.gradient_r(p_estimate_tmp, p_estimate)
                                
        estimate = np.ctypeslib.as_array(p_estimate, shape=(self.n_z, self.n_r))    
            
        for est, tru in zip(estimate.flatten(), truth.flatten()):
            assert(np.abs(est-tru) < self.thresh)
        
        self.plot(estimate, truth) 
    

    def test_gradient_bound(self):
    
        gradz = FinDiff((0, self.dz, 1), acc=2)
        truth_row = gradz(self.data)
        
        gradr = FinDiff((1, self.dr, 1), acc=2)
        truth_col = gradr(self.data)        
        
        truth = np.concatenate([truth_col[:, 0]*self.dz, -truth_row[-1, 1:-1]*self.dr, 
                -truth_col[::-1, -1]*self.dz, truth_row[0, -2:0:-1]*self.dr])
        print(truth.shape)
        estimate = np.ascontiguousarray(np.zeros(2*(self.n_z+self.n_r)-4))
        p_estimate = estimate.ctypes.data_as(ctypes.POINTER(ctypes.c_double))         
               
        c_gradient.gradient_bound(self.p_data, p_estimate)
                                
        estimate = np.ctypeslib.as_array(p_estimate, shape=(2*(self.n_z+self.n_r)-4, ))    
        
        if self.show:
            plt.plot(estimate, '+')
            plt.plot(truth, 'x')
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

 

