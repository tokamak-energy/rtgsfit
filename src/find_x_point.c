#include "find_x_point.h"
#include <stdio.h>
#include <math.h>
#include "constants.h"
#include <float.h>
#include "gradient.h"


#define NP 4

void find_zero_on_edge(
        double* grad_patch,
        double* cross_row, 
        double* cross_col, 
        int* count
        )
{

    int ii;
    const int nn_c[NP] = {1, 1, 0, 0};
    const int nn_r[NP] = {0, 1, 0, 1};
    const int nn_ind[NP] = {1, 3, 0, 2};
    const int rowNum[NP] = {0, 0, 1, 1};
    const int colNum[NP] = {0, 1, 0, 1};
    int diff_col, diff_row;
    double off;

    *count = 0;

    for (ii=0; ii<NP; ii++)
    {
        if (fabs(grad_patch[ii]) < THRESH)
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

void n_comb(
        int n_count,
        int* n_comb,
        int* idx_a,
        int* idx_b
        )
{
    
    switch (n_count)
    {
        case 1:
            *n_comb = 1;
            
            idx_a[0] = 0;
            idx_b[0] = 0;
            
            break;
        case 2:
            *n_comb = 1;
            
            idx_a[0] = 0;
            idx_b[0] = 1;
            
            break;            
        case 3:
            *n_comb = 3;
            
            idx_a[0] = 0;
            idx_b[0] = 1;
            
            idx_a[1] = 1;
            idx_b[1] = 2;
            
            idx_a[2] = 0;
            idx_b[2] = 2;            
            break;  
        case 4:
            *n_comb = 6;
            
            idx_a[0] = 0;
            idx_b[0] = 1;
            
            idx_a[1] = 1;
            idx_b[1] = 2;
            
            idx_a[2] = 2;
            idx_b[2] = 3;
            
            idx_a[3] = 0;
            idx_b[3] = 2;
            
            idx_a[4] = 0;
            idx_b[4] = 3;
            
            idx_a[5] = 1;
            idx_b[5] = 3;                                   
            break;   
    }
}            


double min(
        double aa,
        double bb
        )
{
    if (aa < bb)
    {
        return aa;
    }
    else
    {
        return bb;
    }
}


double max(
        double aa,
        double bb
        )
{
    if (aa > bb)
    {
        return aa;
    }
    else
    {
        return bb;
    }
}


void find_xpoint_cross(
        double a_start_row,
        double a_end_row,
        double a_start_col,
        double a_end_col,
        double b_start_row,
        double b_end_row,
        double b_start_col,
        double b_end_col,
        double* col_off, 
        double* row_off,
        int* flag
        )
{   
    double a_diff_row = a_end_row - a_start_row;
    double a_diff_col = a_end_col - a_start_col;
    double b_diff_row = b_end_row - b_start_row;
    double b_diff_col = b_end_col - b_start_col;
    double det = -a_diff_col*b_diff_row + b_diff_col*a_diff_row;   
    double ab_start_col = a_start_col - b_start_col;
    double ab_start_row = (a_start_row - b_start_row);
    double lambda_a = (b_diff_row * ab_start_col - b_diff_col * ab_start_row) / det;
    double lambda_b = (a_diff_row * ab_start_col  - a_diff_col * ab_start_row) / det;
    
/*    printf("a_start_row %lf, a_end_row %lf, a_start_col %lf, a_end_col %lf, b_start_row %lf, \*/
/*        b_end_row %lf, b_start_col%lf , b_end_col %lf\n", a_start_row, a_end_row, a_start_col, a_end_col, b_start_row,*/
/*        b_end_row, b_start_col, b_end_col);*/
/*    printf("lambda_a %lf, lambda_b %lf\n", lambda_a, lambda_b);*/
    
    *flag = 0;
    if ((lambda_a >= 0) && (lambda_a <= 1) && (lambda_b >= 0) && (lambda_b <= 1))
    {
        *col_off = a_start_col + a_diff_col * lambda_a;
        *row_off =  a_start_row + a_diff_row * lambda_a;
        *flag = 1;
    }
}

/*    double matrix_inv[] = { b_diff_row / det, -b_diff_col / det, */
/*                          -a_diff_row / det, a_diff_col / det };*/
/*   */
/*    double b_vec[] = {a_start_col - b_start_col, a_row_start - b_row_start};*/
/*    double matrix = {a_diff_col, b_diff_col, a_diff_row, b_diff_row};*/
/*    double matrix[] = {a_diff_row, 0.0, -1.0, 0.0,*/
/*                 0.0, b_diff_row, -1.0, 0.0,*/
/*                 a_diff_col, 0.0, 0.0, -1.0,*/
/*                 0.0, b_diff_col, 0.0, -1.0};*/
/*                 */
/*    double vector[] = {*/
/*                 */
                 
                 

void find_null_in_gradient_march(
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
    
    int idx, i_row, i_col, ii, jj;
    double gr_patch[NP], gz_patch[NP], grs_row[NP], gzs_row[NP], grs_col[NP], gzs_col[NP];
    double grad_z[N_GRID], grad_r[N_GRID], hess_zz[N_GRID], hess_rr[N_GRID], hess_rz[N_GRID];
    double hess_det_at_null;
    double hess_rr_at_null;
    double hess_zz_at_null;
    double hess_rz_at_null;    
    int count_r, count_z;
    int is_stnry;
    double col_off, row_off;
    int n_comb_r, n_comb_z;
    int idx_r_start[NP], idx_r_end[NP], idx_z_start[NP], idx_z_end[NP];
    double gr_row_start, gr_row_end, gr_col_start, gr_col_end;
    double gz_row_start, gz_row_end, gz_col_start, gz_col_end;   
     
    *opt_n = 0;
    *xpt_n = 0;
    
    // find hessian & gradients 
    gradient_z(flux, grad_z);   
    gradient_r(flux, grad_r);
    hessian_zz(flux, hess_zz);
    hessian_rr(flux, hess_rr);
    gradient_r(grad_z, hess_rz);
    
    for (i_row = 0; i_row < N_Z-1; i_row++)
    {
        for (i_col = 0; i_col < N_R-1; i_col++)
        {
            idx = i_row*N_R + i_col;
            
            gr_patch[0] = grad_r[idx];
            gr_patch[1] = grad_r[idx+1]; 
            gr_patch[2] = grad_r[idx+N_R];
            gr_patch[3] = grad_r[idx+N_R+1];

            gz_patch[0] = grad_z[idx];
            gz_patch[1] = grad_z[idx+1]; 
            gz_patch[2] = grad_z[idx+N_R];
            gz_patch[3] = grad_z[idx+N_R+1];

            find_zero_on_edge(gr_patch, grs_row, grs_col, &count_r);
            find_zero_on_edge(gz_patch, gzs_row, gzs_col, &count_z);

            if ((count_r > 0) && (count_z > 0))
            {

                n_comb(count_r, &n_comb_r, idx_r_start, idx_r_end);            
                n_comb(count_z, &n_comb_z, idx_z_start, idx_z_end);            
                                                
                for (ii = 0; ii < n_comb_r; ii++)
                {
                    gr_row_start = grs_row[idx_r_start[ii]];
                    gr_row_end = grs_row[idx_r_end[ii]];                
            
                    gr_col_start = grs_col[idx_r_start[ii]];
                    gr_col_end = grs_col[idx_r_end[ii]]; 
                     
                    for (jj = 0; jj < n_comb_z; jj++)
                    {
                    
                        gz_row_start = gzs_row[idx_z_start[ii]];
                        gz_row_end = gzs_row[idx_z_end[ii]];                
                
                        gz_col_start = gzs_col[idx_z_start[ii]];
                        gz_col_end = gzs_col[idx_z_end[ii]];                  
                                           
                        if ((min(gz_row_start, gz_row_end) <= max(gr_row_start, gr_row_end)) &&
                                (min(gr_row_start, gr_row_end) <= max(gz_row_start, gz_row_end)) &&
                                (min(gz_col_start, gz_col_end) <= max(gr_col_start, gr_col_end)) &&
                                (min(gr_col_start, gr_col_end) <= max(gz_col_start, gz_col_end)))
                        {
                        
                            find_xpoint_cross(gz_row_start, gz_row_end, gz_col_start,
                                    gz_col_end, gr_row_start, gr_row_end, gr_col_start,
                                    gr_col_end, &col_off, &row_off, &is_stnry);
                                    
                            if (is_stnry)
                            {

                                hess_rr_at_null = lin_intrp_2(hess_rr, idx, col_off, row_off);
                                hess_zz_at_null = lin_intrp_2(hess_zz, idx, col_off, row_off);
                                hess_rz_at_null = lin_intrp_2(hess_rz, idx, col_off, row_off);
                            
                                hess_det_at_null = hess_rr_at_null*hess_zz_at_null - \
                                        hess_rz_at_null*hess_rz_at_null;                               
                                
                                if (hess_det_at_null > 0.0 && hess_rr_at_null < 0.0 && MASK_LIM[idx]) 
                                {
                                    opt_r[*opt_n] = R_VEC[i_col] + DR*col_off;
                                    opt_z[*opt_n] = Z_VEC[i_row] + DZ*row_off;    
                                    opt_flux[*opt_n] = lin_intrp_2(flux, idx, col_off, row_off);
                                    *opt_n += 1;
                                }
                                else if (hess_det_at_null < 0.0)
                                {
                                    xpt_r[*xpt_n] = R_VEC[i_col] + DR*col_off;
                                    xpt_z[*xpt_n] = Z_VEC[i_row] + DZ*row_off;    
                                    xpt_flux[*xpt_n] = lin_intrp_2(flux, idx, col_off, row_off);
                                    *xpt_n += 1;     
                                }
                            }
                        }
                    }
                }
            }
        }
    }                   
}



double lin_intrp_2(
        double* psi,
        int idx, 
        double frac_dist_null_r,
        double frac_dist_null_z
        )
{
    double psi_up, psi_down, psi_at_null;

    psi_up = (1.0 - frac_dist_null_r)*psi[idx + N_R] + frac_dist_null_r*psi[idx + N_R + 1];
    psi_down = (1.0 - frac_dist_null_r)*psi[idx] + frac_dist_null_r*psi[idx + 1];
    psi_at_null = (1.0 - frac_dist_null_z)*psi_down + frac_dist_null_z*psi_up;
    
    return psi_at_null;
}

    
double lin_intrp(
        double* psi,
        int idx, 
        double dist_to_null_r, 
        double dist_to_null_z,
        double abs_dist_null_r,
        double abs_dist_null_z
        )
{
    double psi_up, psi_down, psi_at_null;

    abs_dist_null_r = abs_dist_null_r / DR;
    abs_dist_null_z = abs_dist_null_z / DZ;

    if (dist_to_null_r > 0.0)
    {
        if (dist_to_null_z > 0.0)
        {
            psi_up = (1 - abs_dist_null_r)*psi[idx + N_R] + abs_dist_null_r*psi[idx + N_R + 1];
            psi_down = (1 - abs_dist_null_r)*psi[idx] + abs_dist_null_r*psi[idx + 1];
            psi_at_null = (1 - abs_dist_null_z)*psi_down + abs_dist_null_z*psi_up;
        }
        else
        {
            psi_up = (1 - abs_dist_null_r)*psi[idx] + abs_dist_null_r*psi[idx + 1];
            psi_down = (1 - abs_dist_null_r)*psi[idx - N_R] + abs_dist_null_r*psi[idx - N_R + 1];
            psi_at_null = (1 - abs_dist_null_z)*psi_up + abs_dist_null_z*psi_down;
        }
    }
    else
    {
        if (dist_to_null_z > 0.0)
        {
            psi_up = (1 - abs_dist_null_r)*psi[idx + N_R] + abs_dist_null_r*psi[idx + N_R - 1];
            psi_down = (1 - abs_dist_null_r)*psi[idx] + abs_dist_null_r*psi[idx - 1];
            psi_at_null = (1 - abs_dist_null_z)*psi_down + abs_dist_null_z*psi_up;
        
        }
        else
        {
            psi_up = (1 - abs_dist_null_r)*psi[idx] + abs_dist_null_r*psi[idx - 1];
            psi_down = (1 - abs_dist_null_r)*psi[idx - N_R] + abs_dist_null_r*psi[idx - N_R - 1];
            psi_at_null = (1 - abs_dist_null_z)*psi_up + abs_dist_null_z*psi_down;                
        }
    }
    return psi_at_null;
}


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
    double dist_to_null_r, dist_to_null_z, abs_dist_null_r, abs_dist_null_z;
    double grad_z[N_GRID], grad_r[N_GRID], hess_zz[N_GRID], hess_rr[N_GRID], hess_rz[N_GRID];
    double hess_det_at_null;
    double hess_rr_at_null;
    double hess_zz_at_null;
    double hess_rz_at_null;    
    
    *opt_n = 0;
    *xpt_n = 0;
    
    // find hessian & gradients 
    gradient_z(flux, grad_z);   
    gradient_r(flux, grad_r);
    hessian_zz(flux, hess_zz);
    hessian_rr(flux, hess_rr);
    gradient_r(grad_z, hess_rz);
    
    for (i_row=1; i_row<(N_Z-1); i_row++) 
    {
        for (i_col=1; i_col<(N_R-1); i_col++)
        {
            
            idx = i_row*N_R + i_col;
            hess_det = hess_zz[idx]*hess_rr[idx] - hess_rz[idx]*hess_rz[idx];
            dist_to_null_r = (grad_z[idx]*hess_rz[idx] - hess_zz[idx]*grad_r[idx])/hess_det;
            dist_to_null_z = (grad_r[idx]*hess_rz[idx] - hess_rr[idx]*grad_z[idx])/hess_det;
            
            abs_dist_null_r = fabs(dist_to_null_r);
            abs_dist_null_z = fabs(dist_to_null_z);
            
            if (abs_dist_null_r < 0.5*DR && abs_dist_null_z < 0.5*DZ && MASK_LIM[idx])
            {   
            
                hess_rr_at_null = lin_intrp(hess_rr, idx, dist_to_null_r, 
                            dist_to_null_z, abs_dist_null_r, abs_dist_null_z);
                hess_zz_at_null = lin_intrp(hess_zz, idx, dist_to_null_r, 
                            dist_to_null_z, abs_dist_null_r, abs_dist_null_z);
                hess_rz_at_null = lin_intrp(hess_rz, idx, dist_to_null_r, 
                            dist_to_null_z, abs_dist_null_r, abs_dist_null_z);  
                            
                hess_det_at_null = hess_rr_at_null*hess_zz_at_null - \
                        hess_rz_at_null*hess_rz_at_null;
                                                                                       
                if (hess_det_at_null > 0.0 && hess_rr_at_null < 0.0) 
                {
                    opt_r[*opt_n] = R_VEC[i_col] + dist_to_null_r;
                    opt_z[*opt_n] = Z_VEC[i_row] + dist_to_null_z;    
                    opt_flux[*opt_n] = lin_intrp(flux, idx, dist_to_null_r, 
                            dist_to_null_z,  abs_dist_null_r, abs_dist_null_z);
                    *opt_n += 1;
                }
                else if (hess_det_at_null < 0.0)
                {
                    xpt_r[*xpt_n] = R_VEC[i_col] + dist_to_null_r;
                    xpt_z[*xpt_n] = Z_VEC[i_row] + dist_to_null_z;    
                    xpt_flux[*xpt_n] = lin_intrp(flux, idx, dist_to_null_r, 
                            dist_to_null_z, abs_dist_null_r, abs_dist_null_z);
                    *xpt_n += 1;     
                }
            }
        }
    }
}


       
void find_lcfs_rz(
        double* flux_orig,
        double flux_lcfs,
        double* lcfs_r,
        double* lcfs_z,
        int* lcfs_n
        )
{
    
    int i_row, i_col, idx;
    double off;
    double flux[N_GRID];
    
    *lcfs_n = 0;
    
    for (idx=0; idx<(N_Z*N_R); idx++)
    {
        flux[idx] = flux_orig[idx] - flux_lcfs;
    }
    
    for (i_row=0; i_row<N_Z-1; i_row++)
    {
        for (i_col=0; i_col<N_R-1; i_col++)
        {
            idx = i_row*N_R + i_col;
            
            if ((fabs(flux[idx]) < THRESH)  && MASK_LIM[idx])
            {

                lcfs_r[*lcfs_n] = R_VEC[i_col];
                lcfs_z[*lcfs_n] = Z_VEC[i_row];
                *lcfs_n = *lcfs_n + 1;
            }
            else 
            {
                if ((flux[idx] * flux[idx+1] < 0)  && (MASK_LIM[idx] || MASK_LIM[idx+1]))
                {
                    off = fabs(flux[idx]) / (fabs(flux[idx]) + fabs(flux[idx+1]));
                    if (((off <= 0.5) && MASK_LIM[idx]) || ((off >= 0.5) && MASK_LIM[idx+1]))
                    {
                        lcfs_r[*lcfs_n] = R_VEC[i_col] + off*DR;
                        lcfs_z[*lcfs_n ] = Z_VEC[i_row];
                        *lcfs_n = *lcfs_n + 1;
                    }
                }
                if ((flux[idx] * flux[idx+N_R] < 0) && (MASK_LIM[idx] || MASK_LIM[idx+N_R]))
                {
                    off = fabs(flux[idx]) / (fabs(flux[idx]) + fabs(flux[idx+N_R]));
                    if (((off <= 0.5) && MASK_LIM[idx]) || ((off >= 0.5) && MASK_LIM[idx+N_R]))
                    {
                        lcfs_r[*lcfs_n] = R_VEC[i_col];
                        lcfs_z[*lcfs_n ] = Z_VEC[i_row] + off*DZ;
                        *lcfs_n = *lcfs_n + 1;
                    }
                }
            }            
         
        }
    }
    
    for (i_row=0; i_row<N_Z-1; i_row++)
    {
        idx = i_row*N_R + N_R - 1;
        
        if ((fabs(flux[idx]) < THRESH) && MASK_LIM[idx])
        {

            lcfs_r[*lcfs_n] = R_VEC[i_col];
            lcfs_z[*lcfs_n] = Z_VEC[i_row];
            *lcfs_n = *lcfs_n + 1;
        }
        else if ((flux[idx] * flux[idx+N_R] < 0) && (MASK_LIM[idx] || MASK_LIM[idx+N_R]))
        {
            off = fabs(flux[idx]) / (fabs(flux[idx]) + fabs(flux[idx+N_R]));
            if (((off <= 0.5) && MASK_LIM[idx]) || ((off >= 0.5) && MASK_LIM[idx+N_R]))
            { 
                lcfs_r[*lcfs_n] = R_VEC[i_col];
                lcfs_z[*lcfs_n ] = Z_VEC[i_row] + off*DZ;
                *lcfs_n = *lcfs_n + 1;
            }
        }            
     
    }
    
    for (i_col=0; i_col<N_R-1; i_col++)
    {
        idx = (N_Z-1)*N_R + i_col;
        
        if ((fabs(flux[idx]) < THRESH) && MASK_LIM[idx])
        {

            lcfs_r[*lcfs_n] = R_VEC[i_col];
            lcfs_z[*lcfs_n] = Z_VEC[i_row];
            *lcfs_n = *lcfs_n + 1;
        }
        else if ((flux[idx] * flux[idx+1] < 0) && (MASK_LIM[idx] || MASK_LIM[idx+1]))
        {
            off = fabs(flux[idx]) / (fabs(flux[idx]) + fabs(flux[idx+1]));
            if (((off <= 0.5) && MASK_LIM[idx]) || ((off >= 0.5) && MASK_LIM[idx+1]))
            {
                lcfs_r[*lcfs_n] = R_VEC[i_col] + off*DR;
                lcfs_z[*lcfs_n ] = Z_VEC[i_row];
                *lcfs_n = *lcfs_n + 1;
            }
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
    
    int i_row, i_col, i_lcfs, col_start, col_end, row_start, row_end, i_grid, i_count;
    int count = 0;
    double z_nearest;
    double r_start, r_end;
    double r_tmp[N_XPT_MAX];
    double z_start, z_end;
    double z_tmp[N_XPT_MAX];

    for (i_grid=0; i_grid<N_GRID; i_grid++)
    {
        mask[i_grid] = 0 ;
    }

    z_nearest = round((z_opt - Z_VEC[0])/DZ)*DZ + Z_VEC[0];
    
    for (i_lcfs=0; i_lcfs<lcfs_n; i_lcfs++)
    {
        if (fabs(lcfs_z[i_lcfs] - z_nearest) < THRESH)
        {
            r_tmp[count] = lcfs_r[i_lcfs];
            count += 1;
        }
    }
    r_start = -DBL_MAX;
    r_end = DBL_MAX;
    
    for (i_count=0; i_count< count; i_count++)
    {
        if (r_tmp[i_count] < r_opt)
        {
            if (r_tmp[i_count] > r_start)
            {
                r_start = r_tmp[i_count];
            }
        }
        else 
        {
            if (r_tmp[i_count] < r_end)
            {
                r_end = r_tmp[i_count];
            }
        }
    }

    col_start = (int) ceil((r_start - R_VEC[0] ) / DR);
    col_end = (int) floor((r_end - R_VEC[0])/DR);

    for (i_col=col_start; i_col<=col_end; i_col++)
    {
        count = 0;
        for (i_lcfs=0; i_lcfs<lcfs_n; i_lcfs++)
        {
            if (fabs(lcfs_r[i_lcfs] - R_VEC[i_col]) < THRESH)
            {
                z_tmp[count] = lcfs_z[i_lcfs];
                count += 1;
            }
        }

        z_start = -DBL_MAX;
        z_end = DBL_MAX;

        for (i_count=0; i_count< count; i_count++)
        {
            if (z_tmp[i_count] < z_opt)
            {
                if (z_tmp[i_count] > z_start)
                {
                    z_start = z_tmp[i_count];
                }
            }
            else 
            {
                if (z_tmp[i_count] < z_end)
                {
                    z_end = z_tmp[i_count];
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
