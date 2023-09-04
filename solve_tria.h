#ifndef SOLVE_TRIA_H_   /* Include guard */
#define SOLVE_TRIA_H_

void permute(double* arr_orig, int n_row, int* idx_final, double* arr_swap);

void back_sub_upper_short(int n_row, int n_rep, double* upper, double* b_vec, 
        double* out);

void back_sub_lower_short(int n_row, int n_rep, double* lower, double* b_vec,
        double* out);
        
void solve_tria(int n_row, int n_rep, double* lower, double* upper, 
        double* b_vec, int* idx_final, double* out);

#endif 
