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

int pulseNo; /* pulse number to replay */
size_t N_iter; /* number of iterations; equal the length of time axis */
size_t dec_N_iter; /* number of decimated iterations */

#endif /* __REPLAY_RTGSFIT__H__ */
