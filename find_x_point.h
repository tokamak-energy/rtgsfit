#ifndef FIND_X_POINT_H_
#define FIND_X_POINT_H_


        
void find_null_in_gradient(double *psi, 
        double *opt_r, double *opt_z, double *opt_psi, int *i_opt,
        double *xpt_r, double *xpt_z, double *xpt_psi, int *i_xpt);
        
                
void find_lcfs_rz(double *psi, double psi_lcfs, 
        double *r_lcfs, double *z_lcfs, int *n_lcfs);
        
                
void inside_lcfs(double r_opt, double z_opt, double *r_lcfs, 
        double *z_lcfs, int n_lcfs, int *mask);           
        

#endif

/*int max_idx(int n_arr, double* arr);*/


/*void find_opt_xpt(double* psi, double* grad_z, double* psi_lcfs,*/
/*        double* psi_axis, double* opt_chosen_r, double* opt_chosen_z);*/

