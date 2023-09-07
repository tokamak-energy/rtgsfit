#ifndef GRADIENT_H_   /* Include guard */
#define GRADIENT_H_

void gradient_row(double* arr, int n_row, int n_col, double d_row, 
        double* grad_row);   

void gradient_col(double* arr, int n_row, int n_col, double d_col,
        double* grad_col);
        
void hessian_row_row(double* arr, int n_row, int n_col, double d_row,
        double* hess_row_row);
        
void hessian_col_col(double* arr, int n_row, int n_col, double d_col,
        double* hess_col_col);
        
void hessian_row_col(double* arr, int n_row, int n_col, double d_row,
        double d_col, double* hess_row_col);
        
void gradient_bound(double* arr, int n_row, int n_col, double d_row, 
        double d_col, double* grad_bound);

#endif 
