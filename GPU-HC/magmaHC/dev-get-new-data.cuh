#ifndef dev_get_new_data_cuh_
#define dev_get_new_data_cuh_
// ==============================================================================================================
// Device function for calculating vectors in Runge-Kutta algorithm and compute solution norm
//
// Modifications
//    Chiang-Heng Chien  21-05-20:   Originally created
//    Chiang-Heng Chien  23-12-28:   Add macros and replace reduction sum by shuffle operation in a warp level
//
//> (c) LEMS, Brown University
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// ==============================================================================================================
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>

#include <cuda_runtime.h>

#include "magma_v2.h"
#include "magma_lapack.h"
#include "magma_internal.h"
#undef max
#undef min
#include "magma_templates.h"
#include "sync.cuh"
#undef max
#undef min
#include "shuffle.cuh"
#undef max
#undef min
#include "batched_kernel_param.h"

#include "definitions.hpp"

__device__ __inline__ void
create_x_for_k2(
  const int tx, float &t, float delta_t, float one_half_delta_t, magmaFloatComplex *s_sols,
  magmaFloatComplex *s_track, magmaFloatComplex *sB, magmaFloatComplex gc )
{
  s_sols[tx] += sB[tx] * delta_t * gc * 1.0/6.0;     // -- s = s + (\Delta t) (k1/6) --
  sB[tx] *= one_half_delta_t * gc;                   // -- k1 * (\Delta t)/2 --
  s_track[tx] += sB[tx];                        // -- x = x + k1 * (\Delta t)/2 --
  t += one_half_delta_t;                        // -- t = t + (\Delta t)/2 --
}

__device__ __inline__ void
create_x_for_k3(
  const int tx, float delta_t, float one_half_delta_t, magmaFloatComplex *s_sols, \
  magmaFloatComplex *s_track, magmaFloatComplex *s_track_pred_init, magmaFloatComplex *sB, \
  magmaFloatComplex gc05 )
{
  s_sols[tx] += sB[tx] * delta_t * gc05 * 1.0/3.0;     // -- s = s + (\Delta t) (k1/6 + k2/3) --
  s_track[tx] = s_track_pred_init[tx];          // -- copy the initial prior prediction solution --
  sB[tx] *= one_half_delta_t * gc05;                   // -- k2 * (\Delta t)/2 --
  s_track[tx] += sB[tx];                        // -- x = x + k2 * (\Delta t)/2 --
}

__device__ __inline__ void
create_x_for_k4(
  const int tx, float &t, float delta_t, float one_half_delta_t, magmaFloatComplex *s_sols,
  magmaFloatComplex *s_track, magmaFloatComplex *s_track_pred_init, magmaFloatComplex *sB, \
  magmaFloatComplex gc05 )
{
  s_sols[tx] += sB[tx] * delta_t * gc05 * 1.0/3.0;     // -- s = s + (\Delta t) (k1/6 + k2/3 + k3/3) --
  s_track[tx] = s_track_pred_init[tx];          // -- copy the initial prior prediction solution --
  sB[tx] *= delta_t * gc05;                            // -- k3 * (\Delta t) --
  s_track[tx] += sB[tx];                        // -- x = x + k3 * (\Delta t) --
  t += one_half_delta_t;                        // -- now t becomes t = t + (\Delta t) --
}


// __device__ __inline__ void
// compute_norm2(
//   const int tx, magmaFloatComplex *sB, magmaFloatComplex *s_track,
//   float *s_sqrt_sols, float *s_sqrt_corr, float* s_norm)
// {
//   s_sqrt_sols[tx] = MAGMA_C_REAL(sB[tx])*MAGMA_C_REAL(sB[tx]) + MAGMA_C_IMAG(sB[tx])*MAGMA_C_IMAG(sB[tx]);
//   s_sqrt_corr[tx] = MAGMA_C_REAL(s_track[tx])*MAGMA_C_REAL(s_track[tx]) + MAGMA_C_IMAG(s_track[tx])*MAGMA_C_IMAG(s_track[tx]);
//   magmablas_syncwarp();
//   magma_sum_reduce< NUM_OF_VARS >( tx, s_sqrt_sols );
//   magma_sum_reduce< NUM_OF_VARS >( tx, s_sqrt_corr );

//   if ( tx == 0 ) {
//       s_norm[0] = s_sqrt_sols[0];
//       s_norm[1] = s_sqrt_corr[0];
//   }
// }

#endif
