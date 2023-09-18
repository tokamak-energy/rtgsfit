#ifndef SOLVE_TRIA_H_   /* Include guard */
#define SOLVE_TRIA_H_

void permute(double* arr_orig, double* arr_swap);

void back_sub_upper_short(double* b_vec, double* out);

void back_sub_lower_short(double* b_vec, double* out);
        
void solve_tria(double* b_vec, double* out);

#endif 
