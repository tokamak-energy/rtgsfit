#ifndef __REPLAY_RTGSFIT__H__
#define __REPLAY_RTGSFIT__H__

#define status_ok(status) (((status) & 1) == 1)
#include <stdio.h>
#include <string.h>
#include <mdslib.h> /* to read/write from/to MDSplus */

//#include "ST40PCS.h"
#include "utils.h"
#include "mds_tools.h"
#include "constants.h"
#include "rtgsfit.h"

char *pcs1 = "acq2106_032";
char *pcs2 = "acq2106_033";
char filename[] = "replay_rtgsfit";
char *decm_sig = "acq2106_032.decimate";
char *fn_sensors_idx = "../data/sensor_index.txt";
char *fn_pf_coil_idx = "../data/pfcoil_index.txt";
char *fn_meas = "../data/meas.txt";
char *fn_coil_curr = "../data/coil_curr.txt";
char *fn_flux_norm = "../data/flux_norm.txt";
char *fn_mask = "../data/mask.txt";
char *fn_psi_total = "../data/psi_total.txt";
float decimate;
int dec;
FILE *p_fid;
uint32_t sensors_num, pf_coils_num;
int signal_name_len = 100;
unsigned int num_pcs1_channels = 160;
unsigned int num_pcs2_channels = 128;
unsigned int num_of_out_channels = 336;
size_t num_of_in_channels = 384;
int dtype_float = DTYPE_FLOAT;
int null = 0;
float *pcs1_data; /* array of floats pcs1 data */
float *pcs1_time; /* array of floats pcs1 time */
float *pcs1_dec_data; /* decimated pcs1 data */
float *pcs1_dec_time; /* decimated pcs1 time axis */
float *pcs_data; /* array of floats pcs1 data */
size_t size_pcs1_data; /* total size of pcs1 data including pcs1 time axis
                          in bytes*/
size_t pcs1_dec_len; /* length of decimated pcs1 data */
size_t size_pcs1_DI;
int len_pcs1_sign; /* length of one pcs1 signal */
int len_pcs1_DI;
size_t pcs1_data_len; /* length of pcs2 data in floats including time axis */
float *pcs2_data; /* array of floats pcs2 data */
float *pcs2_time; /* array of floats used for timebase */
float *pcs2_dec_data; /* decimated pcs2 data */
float *pcs2_dec_time; /* decimated pcs2 time axis */
int len_pcs2_sign; /* length of one pcs2 signal */
int len_pcs2_DI; /* length of 1 DI signal */
size_t size_pcs2_data; /* total size of pcs2 data including pcs2 time axis
                       in bytes */
size_t size_pcs2_DI;
size_t pcs2_data_len; /* length of pcs2 data in floats including time axis */
size_t pcs2_dec_len; /* length of decimated pcs2 data */

uint32_t pulseNo; /* pulse number to replay */
float t_start; /* start time [s] from plasma breakdown */
float rtgsfit_step; /* step time [s] from t_start till t_end at each */
  /* rtgsfit_step rtgsfit will be calculated */
float t_end; /* end time [s] from plasma breakdown */
size_t N_iter; /* number of iterations; equal the length of time axis */
size_t dec_N_iter; /* number of decimated iterations */

#endif /* __REPLAY_RTGSFIT__H__ */
