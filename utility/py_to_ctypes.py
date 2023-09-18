import ctypes
import numpy as np


def py2c(var):
    
    if isinstance(var, int) or isinstance(var, np.int64):
        c_var = ctypes.c_int(np.int32(var))
    elif isinstance(var, float):
        c_var = ctypes.c_double(var)
    elif isinstance(var, np.ndarray):
        if var.dtype == np.float64:
            c_var = var.copy().astype(ctypes.c_double)
            c_var = c_var.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        elif var.dtype == np.int64:
            c_var = var.copy().astype(ctypes.c_int)
            c_var = c_var.ctypes.data_as(ctypes.POINTER(ctypes.c_int))      
        else:
            raise Exception('numpy array type not handled')
    else:
        breakpoint()
        raise Exception('data type not handled')  
        
    return c_var


def c2py(c_var, var):
    
    if isinstance(c_var, ctypes.c_int):
        var = c_var.value
    elif isinstance(c_var, ctypes.c_double):
        var = c_var.value
    else:
        var = np.ctypeslib.as_array(c_var, var.shape)
    return var
    
    
def run_c_func(func, *args):
    
    c_args = [py2c(xx) for xx in args]
    
    func(*c_args)
    
    py_args = [c2py(xx, yy) for xx, yy in zip(c_args, args)]
    
    return py_args
