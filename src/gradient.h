#ifndef GRADIENT_H_   /* Include guard */
#define GRADIENT_H_

void gradient_z(double* arr, double* grad_z);   

void gradient_r(double* arr, double* grad_r);
        
void hessian_zz(double* arr, double* hess_zz);
        
void hessian_rr(double* arr, double* hess_rr);
        
void hessian_rz(double* arr, double* hess_rz);
        
void gradient_bound(double* arr, double* grad_bound);

#endif 
