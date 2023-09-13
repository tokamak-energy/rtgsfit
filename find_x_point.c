#include "find_x_point.h"
#include <stdio.h>

#define MAX_N 32
#define MAX_LCFS_N 2048

int max_idx(n_arr, arr)
{   
    int i_arr;
    int i_max = 0;
    double arr_max = arr[0]; 
    
    for (i_arr=1; i_arr<n_arr; i_arr++)
    {
        if (arr_max < arr[i_arr])
        {
            arr_max = arr[i_arr];
            i_max = i_arr;
        }
    }
    return i_max
} 

void find_opt_xpt(
        double dr,
        double dz,
        int n_row,
        int n_col,
        double* r_vec,
        double* z_vec,
        double* psi,
        double* grad_r,
        double* grad_z,
        double* hess_rr,
        double* hess_zz,
        double* hess_rz,
        double* psi_bound,
        double* psi_axis,
        double* opt_chosen_r,
        double* opt_chosen_z,       
        )
{
    
    double opt_r[MAX_N], opt_z[MAX_N], opt_psi[MAX_N];
    double xpt_r[MAX_N], xpt_z[MAX_N], xpt_psi[MAX_N];
    int n_opt, n_xpt, i_xpt_max, i_opt_max;
    double xpt_psi_max;
    double frac = 0.999;
    
    find_null_in_gradient(dr, dz, n_row, n_col, r_vec, z_vec, psi, grad_r, 
            grad_z, hess_rr, hess_zz, hess_rz, opt_r, opt_z, opt_psi, &n_opt, 
            xpt_r, xpt_z, xpt_psi, &n_xpt);
   
    // extract largest opt/xpt psi values
    i_xpt_max = max_idx(n_xpt, xpt_psi);
    xpt_psi_max = xpt_psi[i_xpt_max];
    
    i_opt_max = max_idx(n_opt, opt_psi);
    
    *psi_axis = opt_psi[i_opt_max];
    *opt_chosen_r = opt_r[i_opt_max];
    *opt_chosen_z = opt_z[i_opt_max];    
    *psi_bound = frac * xpt_psi_max + (1-frac)*psi_axis;
}


void find_lcfs_mask(
        double d_col,
        double d_row,
        int n_row,
        int n_col,
        double* r_vec,
        double* z_vec,
        double* psi,
        double psi_bound,
        int* mask
        )
{
    double thresh = 1e-10;
    double r_lcfs[MAX_LCFS_N], z_lcfs[MAX_LCFS_N];
    int n_lcfs;
    
    // extract LCFS
    find_lcfs(d_col, d_row, n_row, n_col, r_vec, z_vec, psi, psi_bound, thresh, 
            r_lcfs, z_lcfs, &n_lcfs);         
            
    // extract inside of LCFS
    inside_lcfs(dr, dz, n_row, n_col, opt_r, opt_z, thresh, r_vec, z_vec, 
            r_lcfs, z_lcfs, &n_lcfs, mask);
            
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
    double rel_dist_null_r, rel_dist_null_z;

    for (i_row=0; i_row<(n_row-1); i_row++) 
    {
        for (i_col=0; i_col<(n_col-1); i_col++)
        {
            
            idx = i_row*n_col + i_col;
            hess_det = hess_rr[idx]*hess_zz[idx] - hess_rz[idx]*hess_rz[idx];
            
            dist_to_null_r = (grad_z[idx]*hess_rz[idx] - hess_zz[idx]*grad_r[idx])/hess_det;
            dist_to_null_z = (grad_r[idx]*hess_rz[idx] - hess_rr[idx]*grad_z[idx])/hess_det;
            
            if (dist_to_null_r >= 0.0 && dist_to_null_r < dr && 
                    dist_to_null_z >= 0.0 && dist_to_null_z < dz)
            {
                if ((hess_det > 0.0 && hess_rr[idx] < 0.0) || (hess_det < 0.0))
                {
                    rel_dist_null_r = dist_to_null_r / dr;
                    rel_dist_null_z = dist_to_null_z / dz;     
                       
                    // get psi at point from interpolation 
                    if (dist_to_null_r > 0.0)
                    {
                        if (dist_to_null_z > 0.0)
                        {
                            psi_up = (1 - rel_dist_null_r)*psi[idx + n_col] + rel_dist_null_r*psi[idx + n_col + 1];
                            psi_down = (1 - rel_dist_null_r)*psi[idx] + rel_dist_null_r*psi[idx + 1];
                            psi_at_null = (1 - rel_dist_null_z)*psi_down + rel_dist_null_z*psi_up;
                        }
                        else
                        {
                            psi_up = (1 - rel_dist_null_r)*psi[idx] + rel_dist_null_r*psi[idx + 1];
                            psi_up = (1 - rel_dist_null_r)*psi[idx - n_col] + rel_dist_null_r*psi[idx - n_col + 1];
                            psi_at_null = (1 - rel_dist_null_z)*psi_up + rel_dist_null_z*psi_down;
                        }
                    }
                    else
                    {
                        if (dist_to_null_z > 0.0)
                        {
                            psi_up = (1 - rel_dist_null_r)*psi[idx + n_col] + rel_dist_null_r*psi[idx + n_col - 1];
                            psi_down = (1 - rel_dist_null_r)*psi[idx] + rel_dist_null_r*psi[idx - 1];
                            psi_at_null = (1 - rel_dist_null_z)*psi_down + rel_dist_null_z*psi_up;
                        
                        }
                        else
                        {
                            psi_up = (1 - rel_dist_null_r)*psi[idx + n_col] + rel_dist_null_r*psi[idx + n_col - 1];
                            psi_down = (1 - rel_dist_null_r)*psi[idx] + rel_dist_null_r*psi[idx - 1];
                            psi_at_null = (1 - rel_dist_null_z)*psi_up + rel_dist_null_z*psi_down;                   
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
        int* n_lcfs
        )
{
    
    int i_row, i_col, idx;
    double off;
    
    *n_lcfs = 0;
    
    for (idx=0; idx<(n_row*n_col); idx++)
    {
        psi[idx] -= psi_bound;
    }
    
    for (i_row=0; i_row<n_row-1; i_row++)
    {
        for (i_col=0; i_col<n_col-1; i_col++)
        {
            idx = i_row*n_col + i_col;
            
            if (fabs(psi[idx]) < thresh)
            {

                r_lcfs[*n_lcfs] = r_arr[i_col];
                z_lcfs[*n_lcfs] = z_arr[i_row];
                *n_lcfs = *n_lcfs + 1;
            }
            else 
            {
                if (psi[idx] * psi[idx+1] < 0)
                {
                    off = fabs(psi[idx]) / (fabs(psi[idx]) + fabs(psi[idx+1]));
                    r_lcfs[*n_lcfs] = r_arr[i_col] + off*dr;
                    z_lcfs[*n_lcfs ] = z_arr[i_row];
                    *n_lcfs = *n_lcfs + 1;
                }
                if (psi[idx] * psi[idx+n_col] < 0)
                {
                    off = fabs(psi[idx]) / (fabs(psi[idx]) + fabs(psi[idx+n_col]));
                    r_lcfs[*n_lcfs] = r_arr[i_col];
                    z_lcfs[*n_lcfs ] = z_arr[i_row] + off*dz;
                    *n_lcfs = *n_lcfs + 1;
                }
            }            
         
        }
    }
    
    for (i_row=0; i_row<n_row-1; i_row++)
    {
        idx = i_row*n_col + n_col - 1;
        
        if (fabs(psi[idx]) < thresh)
        {

            r_lcfs[*n_lcfs] = r_arr[i_col];
            z_lcfs[*n_lcfs] = z_arr[i_row];
            *n_lcfs = *n_lcfs + 1;
        }
        else if (psi[idx] * psi[idx+n_col] < 0)
        {
            off = fabs(psi[idx]) / (fabs(psi[idx]) + fabs(psi[idx+n_col]));
            r_lcfs[*n_lcfs] = r_arr[i_col];
            z_lcfs[*n_lcfs ] = z_arr[i_row] + off*dz;
            *n_lcfs = *n_lcfs + 1;
        }            
     
    }
    
    for (i_col=0; i_col<n_col-1; i_col++)
    {
        idx = (n_row-1)*n_col + i_col;
        
        if (fabs(psi[idx]) < thresh)
        {

            r_lcfs[*n_lcfs] = r_arr[i_col];
            z_lcfs[*n_lcfs] = z_arr[i_row];
            *n_lcfs = *n_lcfs + 1;
        }
        else if (psi[idx] * psi[idx+1] < 0)
        {
            off = fabs(psi[idx]) / (fabs(psi[idx]) + fabs(psi[idx+1]));
            r_lcfs[*n_lcfs] = r_arr[i_col] + off*dr;
            z_lcfs[*n_lcfs ] = z_arr[i_row];
            *n_lcfs = *n_lcfs + 1;
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
        double* r_lcfs,
        double* z_lcfs,
        int* n_lcfs, 
        int* mask
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
            mask[i_row*n_col + i_col] = 1;
        }
    }
}  
            

/*void inside_lcfs(*/
/*        double dr,*/
/*        double dz,*/
/*        int n_row,*/
/*        int n_col,*/
/*        double r_opt,*/
/*        double z_opt,*/
/*        double thresh,*/
/*        double* r_arr,*/
/*        double* z_arr,*/
/*        double* r_lcfs,*/
/*        double* z_lcfs,*/
/*        int* n_lcfs, */
/*        int* idx,*/
/*        int* n_idx*/
/*        )        */
/*{*/
/*    */
/*    int i_row, i_col, i_lcfs, col_start, col_end, row_start, row_end;*/
/*    int flag = 0;*/
/*    double z_nearest;*/
/*    double r_start, r_end, r_tmp;*/
/*    double z_start, z_end, z_tmp;*/
/*    *n_idx = 0;*/
/*    z_nearest = round((z_opt - z_arr[0])/dz)*dz + z_arr[0];*/
/*    */
/*    for (i_lcfs=0; i_lcfs<*n_lcfs; i_lcfs++)*/
/*    {*/
/*        if (fabs(z_lcfs[i_lcfs] - z_nearest) < thresh)*/
/*        {*/
/*            if (flag == 0)*/
/*            {*/
/*                r_tmp = r_lcfs[i_lcfs];*/
/*                flag += 1;*/
/*            }*/
/*            else if (r_tmp < r_lcfs[i_lcfs])*/
/*            {*/
/*                r_start = r_tmp;*/
/*                r_end = r_lcfs[i_lcfs];*/
/*            }*/
/*            else*/
/*            {*/
/*                r_start = r_lcfs[i_lcfs];*/
/*                r_end = r_tmp;*/
/*            }*/
/*        }*/
/*    }*/
/*    */
/*    col_start = (int) ceil((r_start - r_arr[0] ) / dr);*/
/*    col_end = (int) floor((r_end - r_arr[0])/dr);*/
/*    */
/*    for (i_col=col_start; i_col<=col_end; i_col++)*/
/*    {*/
/*        flag = 0;*/
/*        for (i_lcfs=0; i_lcfs<*n_lcfs; i_lcfs++)*/
/*        {*/
/*            if (fabs(r_lcfs[i_lcfs] - r_arr[i_col]) < thresh)*/
/*            {*/
/*                if (flag == 0)*/
/*                {*/
/*                    z_tmp = z_lcfs[i_lcfs];*/
/*                    flag += 1;*/
/*                }*/
/*                else if (z_tmp < z_lcfs[i_lcfs])*/
/*                {*/
/*                    z_start = z_tmp;*/
/*                    z_end = z_lcfs[i_lcfs];*/
/*                }*/
/*                else*/
/*                {*/
/*                    z_start = z_lcfs[i_lcfs];*/
/*                    z_end = z_tmp;*/
/*                }*/
/*            }*/
/*        }*/
/*        row_start = (int) ceil((z_start - z_arr[0] ) / dz);*/
/*        row_end = (int) floor((z_end - z_arr[0])/dz);*/
/*        */
/*        for (i_row=row_start; i_row<=row_end; i_row++)*/
/*        {*/
/*            idx[*n_idx] = i_row*n_col + i_col;*/
/*            *n_idx += 1;*/
/*        }*/
/*    }*/
/*}  */
/*          */

/*#define NP 4*/

/*void find_zero_on_edge(*/
/*        double* grad_patch,*/
/*        double thresh, */
/*        double* cross_row, */
/*        double* cross_col, */
/*        int* count*/
/*        )*/
/*{*/
/*    */
/*    int ii;*/
/*    int nn_c[NP] = {1, 1, 0, 0};*/
/*    int nn_r[NP] = {0, 1, 0, 1};*/
/*    int nn_ind[NP] = {1, 3, 0, 2};*/
/*    int rowNum[NP] = {0, 0, 1, 1};*/
/*    int colNum[NP] = {0, 1, 0, 1};*/
/*    int diff_col, diff_row;*/
/*    double off;*/
/*    */
/*    *count = 0;*/
/*    */
/*    for (ii=0; ii<NP; ii++)*/
/*        {*/
/*        if (fabs(grad_patch[ii]) < thresh)*/
/*        {*/
/*            cross_col[*count] = colNum[ii];*/
/*            cross_row[*count] = rowNum[ii];*/
/*            *count = *count + 1;*/
/*        }*/
/*        else if (grad_patch[ii] * grad_patch[nn_ind[ii]] < 0)*/
/*        {*/
/*            diff_col = -colNum[ii] + nn_c[ii];*/
/*            diff_row = -rowNum[ii] + nn_r[ii];*/
/*            off = fabs(grad_patch[ii]) / */
/*                  (fabs(grad_patch[ii]) + fabs(grad_patch[nn_ind[ii]]));*/
/*            cross_col[*count] = colNum[ii] + diff_col * off;*/
/*            cross_row[*count] = rowNum[ii] + diff_row * off;*/
/*            *count = *count + 1;*/
/*        }*/
/*    }*/
/*}*/
   
/*void find_lcfs_old(*/
/*        double dr,*/
/*        double dz,*/
/*        int n_row,*/
/*        int n_col,*/
/*        double* r_arr,*/
/*        double* z_arr,*/
/*        double* psi,*/
/*        double psi_bound,*/
/*        double thresh, */
/*        double* r_lcfs,*/
/*        double* z_lcfs,*/
/*        int* n_lcfs*/
/*        )*/
/*{*/
/*    */
/*    int i_row, i_col, idx, count, i_count;*/
/*    double patch[NP], cross_row[n_row*n_col], cross_col[n_row*n_col];*/
/*    */
/*    *n_lcfs = 0;*/
/*    */
/*    for (idx=0; idx<(n_row*n_col); idx++)*/
/*    {*/
/*        psi[idx] -= psi_bound;*/
/*    }*/
/*    */
/*    for (i_row=0; i_row<n_row-1; i_row++)*/
/*    {*/
/*        for (i_col=0; i_col<n_col-1; i_col++)*/
/*        {*/
/*            count = 0;*/
/*            idx = i_row*n_col + i_col;*/
/*            patch[0] = psi[idx];*/
/*            patch[1] = psi[idx+1];*/
/*            patch[2] = psi[idx+n_col];*/
/*            patch[3] = psi[idx+n_col+1];*/
/*            find_zero_on_edge(patch, thresh, cross_row, cross_col, &count);*/
/*            */
/*            for (i_count=0; i_count<count; i_count++)*/
/*            {*/
/*                r_lcfs[*n_lcfs + i_count] = r_arr[i_col] + cross_col[i_count]*dr;*/
/*                z_lcfs[*n_lcfs + i_count] = z_arr[i_row] + cross_row[i_count]*dz;*/
/*            }*/
/*            *n_lcfs += count;*/
/*        }*/
/*    }*/
/*}    */
   
