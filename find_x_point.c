#include "find_x_point.h"
#include <stdio.h>
#include <math.h>
#include "constants.h"
#include <float.h>
#include "gradient.h"

void find_null_in_gradient(
        double* flux,
        double* opt_r,
        double* opt_z,
        double* opt_flux,
        int* opt_n,
        double* xpt_r,
        double* xpt_z,
        double* xpt_flux,
        int* xpt_n       
        )
{
    
    int idx, i_row, i_col;
    double hess_det;
    double flux_up, flux_down, flux_at_null, dist_to_null_r, dist_to_null_z;
    double rel_dist_null_r, rel_dist_null_z;
    double grad_z[N_GRID], grad_r[N_GRID], hess_zz[N_GRID], hess_rr[N_GRID], hess_rz[N_GRID];

    // find hessian & gradients 
    gradient_z(flux, grad_z);   
    gradient_r(flux, grad_r);
    hessian_zz(flux, hess_zz);
    hessian_rr(flux, hess_rr);
    gradient_z(flux, hess_rz);
    gradient_r(flux, hess_rz);
    
    for (i_row=0; i_row<(N_Z-1); i_row++) 
    {
        for (i_col=0; i_col<(N_R-1); i_col++)
        {
            
            idx = i_row*N_R + i_col;
            hess_det = hess_rr[idx]*hess_zz[idx] - hess_rz[idx]*hess_rz[idx];
            
            dist_to_null_r = (grad_z[idx]*hess_rz[idx] - hess_zz[idx]*grad_r[idx])/hess_det;
            dist_to_null_z = (grad_r[idx]*hess_rz[idx] - hess_rr[idx]*grad_z[idx])/hess_det;
            
            if (dist_to_null_r >= 0.0 && dist_to_null_r < DR && 
                    dist_to_null_z >= 0.0 && dist_to_null_z < DZ &&
                    ((hess_det > 0.0 && hess_rr[idx] < 0.0) || (hess_det < 0.0)))
            {
                rel_dist_null_r = dist_to_null_r / DR;
                rel_dist_null_z = dist_to_null_z / DZ;     
                   
                // get flux at point from interpolation 
                flux_up = (1 - rel_dist_null_r)*flux[idx + N_R] + rel_dist_null_r*flux[idx + N_R + 1];
                flux_down = (1 - rel_dist_null_r)*flux[idx] + rel_dist_null_r*flux[idx + 1];
                flux_at_null = (1 - rel_dist_null_z)*flux_down + rel_dist_null_z*flux_up;
                         
                if (hess_det > 0.0 && hess_rr[idx] < 0.0) 
                {
                    opt_r[*opt_n] = R_VEC[i_col] + dist_to_null_r;
                    opt_z[*opt_n] = Z_VEC[i_row] + dist_to_null_z;    
                    opt_flux[*opt_n] = flux_at_null;
                    *opt_n += 1;
                }
                else if (hess_det < 0.0)
                {
                    xpt_r[*xpt_n] = R_VEC[i_col] + dist_to_null_r;
                    xpt_z[*xpt_n] = Z_VEC[i_row] + dist_to_null_z;    
                    xpt_flux[*xpt_n] = flux_at_null;
                    *xpt_n += 1;      
                }
            }
        }
    }
}
                  
       
void find_lcfs_rz(
        double* flux,
        double flux_lcfs,
        double* lcfs_r,
        double* lcfs_z,
        int* lcfs_n
        )
{
    
    int i_row, i_col, idx;
    double off;
    
    *lcfs_n = 0;
    
    for (idx=0; idx<(N_Z*N_R); idx++)
    {
        flux[idx] -= flux_lcfs;
    }
    
    for (i_row=0; i_row<N_Z-1; i_row++)
    {
        for (i_col=0; i_col<N_R-1; i_col++)
        {
            idx = i_row*N_R + i_col;
            
            if (fabs(flux[idx]) < DBL_MIN)
            {

                lcfs_r[*lcfs_n] = R_VEC[i_col];
                lcfs_z[*lcfs_n] = Z_VEC[i_row];
                *lcfs_n = *lcfs_n + 1;
            }
            else 
            {
                if (flux[idx] * flux[idx+1] < 0)
                {
                    off = fabs(flux[idx]) / (fabs(flux[idx]) + fabs(flux[idx+1]));
                    lcfs_r[*lcfs_n] = R_VEC[i_col] + off*DR;
                    lcfs_z[*lcfs_n ] = Z_VEC[i_row];
                    *lcfs_n = *lcfs_n + 1;
                }
                if (flux[idx] * flux[idx+N_R] < 0)
                {
                    off = fabs(flux[idx]) / (fabs(flux[idx]) + fabs(flux[idx+N_R]));
                    lcfs_r[*lcfs_n] = R_VEC[i_col];
                    lcfs_z[*lcfs_n ] = Z_VEC[i_row] + off*DZ;
                    *lcfs_n = *lcfs_n + 1;
                }
            }            
         
        }
    }
    
    for (i_row=0; i_row<N_Z-1; i_row++)
    {
        idx = i_row*N_R + N_R - 1;
        
        if (fabs(flux[idx]) < DBL_MIN)
        {

            lcfs_r[*lcfs_n] = R_VEC[i_col];
            lcfs_z[*lcfs_n] = Z_VEC[i_row];
            *lcfs_n = *lcfs_n + 1;
        }
        else if (flux[idx] * flux[idx+N_R] < 0)
        {
            off = fabs(flux[idx]) / (fabs(flux[idx]) + fabs(flux[idx+N_R]));
            lcfs_r[*lcfs_n] = R_VEC[i_col];
            lcfs_z[*lcfs_n ] = Z_VEC[i_row] + off*DZ;
            *lcfs_n = *lcfs_n + 1;
        }            
     
    }
    
    for (i_col=0; i_col<N_R-1; i_col++)
    {
        idx = (N_Z-1)*N_R + i_col;
        
        if (fabs(flux[idx]) < DBL_MIN)
        {

            lcfs_r[*lcfs_n] = R_VEC[i_col];
            lcfs_z[*lcfs_n] = Z_VEC[i_row];
            *lcfs_n = *lcfs_n + 1;
        }
        else if (flux[idx] * flux[idx+1] < 0)
        {
            off = fabs(flux[idx]) / (fabs(flux[idx]) + fabs(flux[idx+1]));
            lcfs_r[*lcfs_n] = R_VEC[i_col] + off*DR;
            lcfs_z[*lcfs_n ] = Z_VEC[i_row];
            *lcfs_n = *lcfs_n + 1;
        }
    }
} 
       
              
void inside_lcfs(
        double r_opt,
        double z_opt,
        double* lcfs_r,
        double* lcfs_z,
        int lcfs_n, 
        int* mask
        )        
{
    
    int i_row, i_col, i_lcfs, col_start, col_end, row_start, row_end;
    int flag = 0;
    double z_nearest;
    double r_start, r_end, r_tmp;
    double z_start, z_end, z_tmp;

    z_nearest = round((z_opt - Z_VEC[0])/DZ)*DZ + Z_VEC[0];
    
    for (i_lcfs=0; i_lcfs<lcfs_n; i_lcfs++)
    {
        if (fabs(lcfs_z[i_lcfs] - z_nearest) < DBL_MIN)
        {
            if (flag == 0)
            {
                r_tmp = lcfs_r[i_lcfs];
                flag += 1;
            }
            else if (r_tmp < lcfs_r[i_lcfs])
            {
                r_start = r_tmp;
                r_end = lcfs_r[i_lcfs];
            }
            else
            {
                r_start = lcfs_r[i_lcfs];
                r_end = r_tmp;
            }
        }
    }
    
    col_start = (int) ceil((r_start - R_VEC[0] ) / DR);
    col_end = (int) floor((r_end - R_VEC[0])/DR);
    
    for (i_col=col_start; i_col<=col_end; i_col++)
    {
        flag = 0;
        for (i_lcfs=0; i_lcfs<lcfs_n; i_lcfs++)
        {
            if (fabs(lcfs_r[i_lcfs] - R_VEC[i_col]) < DBL_MIN)
            {
                if (flag == 0)
                {
                    z_tmp = lcfs_z[i_lcfs];
                    flag += 1;
                }
                else if (z_tmp < lcfs_z[i_lcfs])
                {
                    z_start = z_tmp;
                    z_end = lcfs_z[i_lcfs];
                }
                else
                {
                    z_start = lcfs_z[i_lcfs];
                    z_end = z_tmp;
                }
            }
        }
        row_start = (int) ceil((z_start - Z_VEC[0] ) / DZ);
        row_end = (int) floor((z_end - Z_VEC[0])/DZ);
        
        for (i_row=row_start; i_row<=row_end; i_row++)
        {
            mask[i_row*N_R + i_col] = 1;
        }
    }
}  
            



/*void find_opt_xpt(      */
/*        double* flux,*/
/*        double* grad_z,*/
/*        double* flux_lcfs,*/
/*        double* flux_axis,*/
/*        double* opt_chosen_r,*/
/*        double* opt_chosen_z,       */
/*        )*/
/*{*/
/*    */
/*    double opt_r[N_MAX_XPT], opt_z[N_MAX_XPT], opt_flux[N_MAX_XPT];*/
/*    double xpt_r[N_MAX_XPT], xpt_z[N_MAX_XPT], xpt_flux[N_MAX_XPT];*/
/*    int n_opt, n_xpt, i_xpt_max, i_opt_max;*/
/*    double xpt_flux_max;*/
/*    double frac = 0.99;*/
/*    */
/*    find_null_in_gradient(flux, grad_z, opt_r, opt_z, opt_flux, &n_opt, xpt_r, */
/*            xpt_z, xpt_flux, &n_xpt);*/
/*   */
/*    // extract largest opt/xpt flux values*/
/*    i_xpt_max = max_idx(n_xpt, xpt_flux);*/
/*    xpt_flux_max = xpt_flux[i_xpt_max];*/
/*    */
/*    i_opt_max = max_idx(n_opt, opt_flux);    */
/*    *flux_axis = opt_flux[i_opt_max];*/
/*    *opt_chosen_r = opt_r[i_opt_max];*/
/*    *opt_chosen_z = opt_z[i_opt_max];  */
/*      */
/*    *flux_lcfs = frac * xpt_flux_max + (1-frac)*flux_axis;*/
/*}*/
/*void inside_lcfs(*/
/*        double DR,*/
/*        double DZ,*/
/*        int N_Z,*/
/*        int N_R,*/
/*        double r_opt,*/
/*        double z_opt,*/
/*        double DBL_MIN,*/
/*        double* R_VEC,*/
/*        double* Z_VEC,*/
/*        double* lcfs_r,*/
/*        double* lcfs_z,*/
/*        int* lcfs_n, */
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
/*    z_nearest = round((z_opt - Z_VEC[0])/DZ)*DZ + Z_VEC[0];*/
/*    */
/*    for (i_lcfs=0; i_lcfs<*lcfs_n; i_lcfs++)*/
/*    {*/
/*        if (fabs(lcfs_z[i_lcfs] - z_nearest) < DBL_MIN)*/
/*        {*/
/*            if (flag == 0)*/
/*            {*/
/*                r_tmp = lcfs_r[i_lcfs];*/
/*                flag += 1;*/
/*            }*/
/*            else if (r_tmp < lcfs_r[i_lcfs])*/
/*            {*/
/*                r_start = r_tmp;*/
/*                r_end = lcfs_r[i_lcfs];*/
/*            }*/
/*            else*/
/*            {*/
/*                r_start = lcfs_r[i_lcfs];*/
/*                r_end = r_tmp;*/
/*            }*/
/*        }*/
/*    }*/
/*    */
/*    col_start = (int) ceil((r_start - R_VEC[0] ) / DR);*/
/*    col_end = (int) floor((r_end - R_VEC[0])/DR);*/
/*    */
/*    for (i_col=col_start; i_col<=col_end; i_col++)*/
/*    {*/
/*        flag = 0;*/
/*        for (i_lcfs=0; i_lcfs<*lcfs_n; i_lcfs++)*/
/*        {*/
/*            if (fabs(lcfs_r[i_lcfs] - R_VEC[i_col]) < DBL_MIN)*/
/*            {*/
/*                if (flag == 0)*/
/*                {*/
/*                    z_tmp = lcfs_z[i_lcfs];*/
/*                    flag += 1;*/
/*                }*/
/*                else if (z_tmp < lcfs_z[i_lcfs])*/
/*                {*/
/*                    z_start = z_tmp;*/
/*                    z_end = lcfs_z[i_lcfs];*/
/*                }*/
/*                else*/
/*                {*/
/*                    z_start = lcfs_z[i_lcfs];*/
/*                    z_end = z_tmp;*/
/*                }*/
/*            }*/
/*        }*/
/*        row_start = (int) ceil((z_start - Z_VEC[0] ) / DZ);*/
/*        row_end = (int) floor((z_end - Z_VEC[0])/DZ);*/
/*        */
/*        for (i_row=row_start; i_row<=row_end; i_row++)*/
/*        {*/
/*            idx[*n_idx] = i_row*N_R + i_col;*/
/*            *n_idx += 1;*/
/*        }*/
/*    }*/
/*}  */
/*          */

/*#define NP 4*/

/*void find_zero_on_edge(*/
/*        double* grad_patch,*/
/*        double DBL_MIN, */
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
/*        if (fabs(grad_patch[ii]) < DBL_MIN)*/
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
/*        double DR,*/
/*        double DZ,*/
/*        int N_Z,*/
/*        int N_R,*/
/*        double* R_VEC,*/
/*        double* Z_VEC,*/
/*        double* flux,*/
/*        double flux_lcfs,*/
/*        double DBL_MIN, */
/*        double* lcfs_r,*/
/*        double* lcfs_z,*/
/*        int* lcfs_n*/
/*        )*/
/*{*/
/*    */
/*    int i_row, i_col, idx, count, i_count;*/
/*    double patch[NP], cross_row[N_Z*N_R], cross_col[N_Z*N_R];*/
/*    */
/*    *lcfs_n = 0;*/
/*    */
/*    for (idx=0; idx<(N_Z*N_R); idx++)*/
/*    {*/
/*        flux[idx] -= flux_lcfs;*/
/*    }*/
/*    */
/*    for (i_row=0; i_row<N_Z-1; i_row++)*/
/*    {*/
/*        for (i_col=0; i_col<N_R-1; i_col++)*/
/*        {*/
/*            count = 0;*/
/*            idx = i_row*N_R + i_col;*/
/*            patch[0] = flux[idx];*/
/*            patch[1] = flux[idx+1];*/
/*            patch[2] = flux[idx+N_R];*/
/*            patch[3] = flux[idx+N_R+1];*/
/*            find_zero_on_edge(patch, DBL_MIN, cross_row, cross_col, &count);*/
/*            */
/*            for (i_count=0; i_count<count; i_count++)*/
/*            {*/
/*                lcfs_r[*lcfs_n + i_count] = R_VEC[i_col] + cross_col[i_count]*DR;*/
/*                lcfs_z[*lcfs_n + i_count] = Z_VEC[i_row] + cross_row[i_count]*DZ;*/
/*            }*/
/*            *lcfs_n += count;*/
/*        }*/
/*    }*/
/*}    */
   
