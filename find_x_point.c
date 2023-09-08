#include "find_x_point.h"
#include <stdio.h>

#define NP 4

void find_zero_on_edge(
        double* grad_patch,
        double thresh, 
        double* cross_row, 
        double* cross_col, 
        int* count
        )
{
    
    int ii;
    int nn_c[NP] = {1, 1, 0, 0};
    int nn_r[NP] = {0, 1, 0, 1};
    int nn_ind[NP] = {1, 3, 0, 2};
    int rowNum[NP] = {0, 0, 1, 1};
    int colNum[NP] = {0, 1, 0, 1};
    int diff_col, diff_row;
    double off;
    
    *count = 0;
    
    for (ii=0; ii<NP; ii++)
        {
        if (fabs(grad_patch[ii]) < thresh)
        {
            cross_col[*count] = colNum[ii];
            cross_row[*count] = rowNum[ii];
            *count = *count + 1;
        }
        else if (grad_patch[ii] * grad_patch[nn_ind[ii]] < 0)
        {
            diff_col = -colNum[ii] + nn_c[ii];
            diff_row = -rowNum[ii] + nn_r[ii];
            off = fabs(grad_patch[ii]) / 
                  (fabs(grad_patch[ii]) + fabs(grad_patch[nn_ind[ii]]));
            cross_col[*count] = colNum[ii] + diff_col * off;
            cross_row[*count] = rowNum[ii] + diff_row * off;
            *count = *count + 1;
        }
    }
}

void find_null_in_gradient(
        double dr,
        double dz,
        int n_row,
        int n_col,
        double* r_arr,
        double* z_arr,
        double* psi,
        double* grad_r,
        double* grad_z,
        double* hess_rr,
        double* hess_zz,
        double* hess_rz,
        double* opt_r,
        double* opt_z,
        double* opt_psi,
        int* i_opt,
        double* xpt_r,
        double* xpt_z,
        double* xpt_psi,
        int* i_xpt       
        )
{
    
    int idx, i_row, i_col;
    double hess_det;
    double psi_up, psi_down, psi_at_null, dist_to_null_r, dist_to_null_z;
    double abs_dist_null_r, abs_dist_null_z;
    double inv_dr = 1.0/dr;
    double inv_dz = 1.0/dz;
    
    for (i_row=0; i_row<n_row; i_row++) 
    {
        for (i_col=0; i_col<n_col; i_col++)
        {
            
            idx = i_row*n_col + i_col;
            hess_det = hess_rr[idx]*hess_zz[idx] - hess_rz[idx]*hess_rz[idx];
            
            dist_to_null_r = (grad_z[idx]*hess_rz[idx] - hess_zz[idx]*grad_r[idx])/hess_det;
            dist_to_null_z = (grad_r[idx]*hess_rz[idx] - hess_rr[idx]*grad_z[idx])/hess_det;
            
            abs_dist_null_r = fabs(dist_to_null_r);
            abs_dist_null_z = fabs(dist_to_null_z);
            
            if (abs_dist_null_r < 0.5*dr && abs_dist_null_z < 0.5*dz)
            {
                if ((hess_det > 0.0 && hess_rr[idx] < 0.0) || (hess_det < 0.0))
                {
                    // get psi at point from interpolation 
                    if (dist_to_null_r > 0.0)
                    {
                        if (dist_to_null_z > 0.0)
                        {
                            psi_up = inv_dr*((1 - abs_dist_null_r)*psi[idx + n_col] + abs_dist_null_r*psi[idx + n_col + 1]);
                            psi_down = inv_dr*((1 - abs_dist_null_r)*psi[idx] + abs_dist_null_r*psi[idx + 1]);
                            psi_at_null = inv_dz*((1 - abs_dist_null_z)*psi_down + abs_dist_null_z*psi_up);
                        }
                        else
                        {
                            psi_up = inv_dr*((1 - abs_dist_null_r)*psi[idx] + abs_dist_null_r*psi[idx + 1]);
                            psi_up = inv_dr*((1 - abs_dist_null_r)*psi[idx - n_col] + abs_dist_null_r*psi[idx - n_col + 1]);
                            psi_at_null = inv_dz*((1 - abs_dist_null_z)*psi_up + abs_dist_null_z*psi_down);
                        }
                    }
                    else
                    {
                        if (dist_to_null_z > 0.0)
                        {
                            psi_up = inv_dr*((1 - abs_dist_null_r)*psi[idx + n_col] + abs_dist_null_r*psi[idx + n_col - 1]);
                            psi_down = inv_dr*((1 - abs_dist_null_r)*psi[idx] + abs_dist_null_r*psi[idx - 1]);
                            psi_at_null = inv_dz*((1 - abs_dist_null_z)*psi_down + abs_dist_null_z*psi_up);
                        
                        }
                        else
                        {
                            psi_up = inv_dr*((1 - abs_dist_null_r)*psi[idx + n_col] + abs_dist_null_r*psi[idx + n_col - 1]);
                            psi_down = inv_dr*((1 - abs_dist_null_r)*psi[idx] + abs_dist_null_r*psi[idx - 1]);
                            psi_at_null = inv_dz*((1 - abs_dist_null_z)*psi_up + abs_dist_null_z*psi_down);                   
                        }
                    }
                    
                    if (hess_det > 0.0 && hess_rr[idx] < 0.0) 
                    {
                        opt_r[*i_opt] = r_arr[i_col] + dist_to_null_r;
                        opt_z[*i_opt] = z_arr[i_row] + dist_to_null_z;    
                        opt_psi[*i_opt] = psi_at_null;
                        *i_opt += 1;
                    }
                    else if (hess_det < 0.0)
                    {
                        xpt_r[*i_xpt] = r_arr[i_col] + dist_to_null_r;
                        xpt_z[*i_xpt] = z_arr[i_row] + dist_to_null_z;    
                        xpt_psi[*i_xpt] = psi_at_null;
                        *i_xpt += 1;      
                    }              
                    
                } 
            }
        }
    }
}
                
void find_lcfs(
        double dr,
        double dz,
        int n_row,
        int n_col,
        double* r_arr,
        double* z_arr,
        double* psi,
        double psi_bound,
        double thresh, 
        double* r_lcfs,
        double* z_lcfs,
        double* n_lcfs
        )
{
    
    int i_row, i_col, idx, count, i_count;
    double patch[NP], cross_row[n_row*n_col], cross_col[n_row*n_col];
    
    n_lcfs = 0;
    
    for (idx=0; idx<(n_row*n_col); idx++)
    {
        psi[idx] -= psi_bound;
    }
    
    for (i_row=0; i_row<n_row-1; i_row++)
    {
        for (i_col=0; i_col<n_col-1; i_col++)
        {
            count = 0;
            idx = i_row*n_col + i_col;
            patch[0] = psi[idx];
            patch[1] = psi[idx+1];
            patch[2] = psi[idx+n_col];
            patch[3] = psi[idx+n_col+1];
            find_zero_on_edge(patch, thresh, cross_row, cross_col, &count);
            n_lcfs += count;
            
            for (i_count=0; i_count<count; i_count++)
            {
                r_lcfs[i_count] = r_arr[i_col] + cross_col[i_count]*dr;
                z_lcfs[i_count] = z_arr[i_row] + cross_row[i_count]*dz;
            }
        }
    }
} 
       
       
void inside_lcfs(
        double dr,
        double dz,
        int n_row,
        int n_col,
        double r_opt,
        double z_opt,
        double thresh,
        double* r_arr,
        double* z_arr,
        double* psi,
        double psi_bound,
        double* r_lcfs,
        double* z_lcfs,
        double* n_lcfs, 
        int* idx,
        int* n_idx
        )        
{
    
    int i_row, i_col, i_lcfs, col_start, col_end, row_start, row_end;
    int flag = 0;
    double z_nearest;
    double r_start, r_end, r_tmp;
    double z_start, z_end, z_tmp;
    *n_idx = 0;
    z_nearest = round((z_opt - z_arr[0])/dz)*dz + z_arr[0];
    
    for (i_lcfs=0; i_lcfs<*n_lcfs; i_lcfs++)
    {
        if (fabs(z_lcfs[i_lcfs] - z_nearest) < thresh)
        {
            if (flag == 0)
            {
                r_tmp = r_lcfs[i_lcfs];
                flag += 1;
            }
            else if (r_tmp < r_lcfs[i_lcfs])
            {
                r_start = r_tmp;
                r_end = r_lcfs[i_lcfs];
            }
            else
            {
                r_start = r_lcfs[i_lcfs];
                r_end = r_tmp;
            }
        }
    }
    
    col_start = (int) ceil((r_start - r_arr[0] ) / dr);
    col_end = (int) floor((r_end - r_arr[0])/dr);
    
    for (i_col=col_start; i_col<=col_end; i_col++)
    {
        flag = 0;
        for (i_lcfs=0; i_lcfs<*n_lcfs; i_lcfs++)
        {
            if (fabs(r_lcfs[i_lcfs] - r_arr[i_col]) < thresh)
            {
                if (flag == 0)
                {
                    z_tmp = z_lcfs[i_lcfs];
                    flag += 1;
                }
                else if (z_tmp < z_lcfs[i_lcfs])
                {
                    z_start = z_tmp;
                    z_end = z_lcfs[i_lcfs];
                }
                else
                {
                    z_start = z_lcfs[i_lcfs];
                    z_end = z_tmp;
                }
            }
        }
        row_start = (int) ceil((z_start - z_arr[0] ) / dz);
        row_end = (int) floor((z_end - z_arr[0])/dz);
        
        for (i_row=row_start; i_row<=row_end; i_row++)
        {
            idx[*n_idx] = i_row*n_col + i_col;
            *n_idx += 1;
        }
    }
}  
            
   
   
   
