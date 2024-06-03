#ifndef KERNEL_HC_solver_5PT_REL_POS_FULL_FORM_QUAT_CU
#define KERNEL_HC_solver_5PT_REL_POS_FULL_FORM_QUAT_CU
// =============================================================================================================
// GPU homotopy continuation solver for 5-point relative pose, full form, Quaternion representation
//
// Modifications
//    Chiang-Heng Chien  24-05-25:   Build on top of 5pt_rel_pos_alg_form_quat
//
//> (c) LEMS, Brown University
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// =============================================================================================================
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>

#include <cuda.h>
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

#include "../definitions.hpp"
#include "magmaHC-kernels.hpp"

//> device functions
#include "../gpu-idx-evals/dev-eval-indxing-5pt_rel_pos_full_form_quat.cuh"
#include "../dev-cgesv-batched-small.cuh"
#include "../dev-get-new-data.cuh"

template< int Num_Of_Vars, int Num_of_Coeffs_from_Params, int Max_Order_of_t, \
          int dHdx_Max_Terms, int dHdx_Max_Parts, int dHdx_Entry_Offset, int dHdx_Row_Offset, \
          int dHdt_Max_Terms, int dHdt_Max_Parts, int dHdt_Row_Offset, \
          unsigned Full_Parallel_Offset, \
          unsigned Partial_Parallel_Thread_Offset, \
          unsigned Partial_Parallel_Index_Offset, \
          unsigned Max_Order_of_t_Plus_One, \
          unsigned Partial_Parallel_Index_Offset_Hx, \
          unsigned Partial_Parallel_Index_Offset_Ht >
__global__ void
HC_solver_5pt_rel_pos_full_form_quat(
  int HC_max_steps, int HC_max_correction_steps, int HC_delta_t_incremental_steps,
  magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
  magma_int_t* d_Hx_indices, magma_int_t* d_Ht_indices,
  magmaFloatComplex_ptr d_phc_coeffs_Hx, magmaFloatComplex_ptr d_phc_coeffs_Ht,
  bool* d_is_GPU_HC_Sol_Converge, bool* d_is_GPU_HC_Sol_Infinity,
  magmaFloatComplex* d_Debug_Purpose
)
{
  extern __shared__ magmaFloatComplex zdata[];
  const int tx = threadIdx.x;
  const int batchid = blockIdx.x ;

  magmaFloatComplex* d_startSols   = d_startSols_array[batchid];
  magmaFloatComplex* d_track       = d_Track_array[batchid];
  const int* __restrict__ d_Hx_idx = d_Hx_indices;
  const int* __restrict__ d_Ht_idx = d_Ht_indices;
  const magmaFloatComplex* __restrict__ d_const_phc_coeffs_Hx = d_phc_coeffs_Hx;
  // const magmaFloatComplex* __restrict__ d_const_phc_coeffs_Ht = d_phc_coeffs_Ht;

  //> registers declarations
  magmaFloatComplex r_cgesvA[Num_Of_Vars] = {MAGMA_C_ZERO};
  magmaFloatComplex r_cgesvB = MAGMA_C_ZERO;
  int linfo = 0, rowid = tx;
  float t0 = 0.0, t_step = 0.0, delta_t = 0.05;
  bool end_zone = 0;

  //> shared memory declarations
  magmaFloatComplex *s_sols               = (magmaFloatComplex*)(zdata);
  magmaFloatComplex *s_track              = s_sols                   + (Num_Of_Vars+1);
  magmaFloatComplex *s_track_last_success = s_track                  + (Num_Of_Vars+1);
  magmaFloatComplex *sB                   = s_track_last_success     + (Num_Of_Vars+1);
  magmaFloatComplex *sx                   = sB                       + Num_Of_Vars;
  magmaFloatComplex *s_phc_coeffs_Hx      = sx                       + Num_Of_Vars;
  magmaFloatComplex *s_phc_coeffs_Ht      = s_phc_coeffs_Hx          + (Num_of_Coeffs_from_Params+1);
  float* dsx                              = (float*)(s_phc_coeffs_Ht + (Num_of_Coeffs_from_Params+1));
  int* sipiv                              = (int*)(dsx               + Num_Of_Vars);
#if USE_LOOPY_RUNGE_KUTTA
  float* s_delta_t_scale                  = (float*)(sipiv           + (Num_Of_Vars+1));
  int* s_RK_Coeffs                        = (int*)(s_delta_t_scale   + 1);
#endif

  s_sols[tx] = d_startSols[tx];
  s_track[tx] = d_track[tx];
  s_track_last_success[tx] = s_track[tx];
  if (tx == 0) {
    s_sols[Num_Of_Vars]               = MAGMA_C_MAKE(1.0, 0.0);
    s_track[Num_Of_Vars]              = MAGMA_C_MAKE(1.0, 0.0);
    s_track_last_success[Num_Of_Vars] = MAGMA_C_MAKE(1.0, 0.0);
    sipiv[Num_Of_Vars]                = 0;
  }
  magmablas_syncwarp();

  float one_half_delta_t;   //> 1/2 \Delta t
  float r_sqrt_sols;
  float r_sqrt_corr;
  bool r_isSuccessful;
  bool r_isInfFail;

  //> This is special for this problem as the polynomial degree = 1, rendering dH/dt independent of t
  s_phc_coeffs_Ht[tx] = d_phc_coeffs_Ht[tx];
  if (tx < 5) {
    s_phc_coeffs_Ht[tx + Num_Of_Vars] = d_phc_coeffs_Ht[tx + Num_Of_Vars];
  }

  //#pragma unroll
  volatile int hc_max_steps = HC_max_steps;
  for (int step = 0; step <= hc_max_steps; step++) {
    if (t0 < 1.0 && (1.0-t0 > 0.0000001)) {

      // ===================================================================
      //> Decide delta t at end zone
      // ===================================================================
      if (!end_zone && fabs(1 - t0) <= (0.0500001)) {
        end_zone = true;
      }

      if (end_zone) {
        if (delta_t > fabs(1 - t0))
          delta_t = fabs(1 - t0);
      }
      else if (delta_t > fabs(1 - 0.05 - t0)) {
        delta_t = fabs(1 - 0.05 - t0);
      }

      t_step = t0;
      one_half_delta_t = 0.5 * delta_t;

      // ===================================================================
      //> Runge-Kutta Predictor
      // ===================================================================
      unsigned char scales[3] = {1, 0, 1};
      if (tx == 0) {
        s_delta_t_scale[0] = 0.0;
        s_RK_Coeffs[0] = 1;
      }
      magmablas_syncwarp();

      //> For simplicity, let's stay with no gamma-trick mode
      #pragma no unroll
      for (int rk_step = 0; rk_step < 4; rk_step++ ) {

        //> Evaluate parameter homotopy
        eval_parameter_homotopy< Num_Of_Vars, Max_Order_of_t, Full_Parallel_Offset, Partial_Parallel_Thread_Offset, Partial_Parallel_Index_Offset, \
                                Max_Order_of_t_Plus_One, Partial_Parallel_Index_Offset_Hx, Partial_Parallel_Index_Offset_Ht > \
                                ( tx, t0, s_phc_coeffs_Hx, d_const_phc_coeffs_Hx );

        //> Evaluate dH/dx and dH/dt
        eval_Jacobian_Hx< Num_Of_Vars, dHdx_Max_Terms, dHdx_Max_Parts, dHdx_Entry_Offset, dHdx_Row_Offset >( tx, s_track, r_cgesvA, d_Hx_idx, s_phc_coeffs_Hx );
        eval_Jacobian_Ht< dHdt_Max_Terms, dHdt_Max_Parts, dHdt_Row_Offset >( tx, s_track, r_cgesvB, d_Ht_idx, s_phc_coeffs_Ht );

        //> linear system solver: solve for k1, k2, k3, or k4
        cgesv_batched_small_device< Num_Of_Vars >( tx, r_cgesvA, sipiv, r_cgesvB, sB, sx, dsx, rowid, linfo );
        magmablas_syncwarp();

        if (rk_step < 3) {

          s_sols[tx] += sB[tx] * delta_t * (s_RK_Coeffs[0] * 1.0/6.0);
          s_track[tx] = (s_RK_Coeffs[0] > 1) ? s_track_last_success[tx] : s_track[tx];

          if (tx == 0) {
            s_delta_t_scale[0] += scales[rk_step] * one_half_delta_t;
            s_RK_Coeffs[0] = s_RK_Coeffs[0] << scales[rk_step];           //> Shift one bit
          }
          magmablas_syncwarp();

          sB[tx] *= s_delta_t_scale[0];
          s_track[tx] += sB[tx];
          t0 += scales[rk_step] * one_half_delta_t;
        }
        magmablas_syncwarp();
      }
      //> Make prediction
      s_sols[tx] += sB[tx] * delta_t * 1.0/6.0;
      s_track[tx] = s_sols[tx];
      magmablas_syncwarp();

      // ===================================================================
      //> Gauss-Newton Corrector
      // ===================================================================
      //#pragma unroll
      volatile int hc_max_xorrection_steps = HC_max_correction_steps;
      #pragma no unroll
      for(int i = 0; i < hc_max_xorrection_steps; i++) {

        eval_Jacobian_Hx< Num_Of_Vars, dHdx_Max_Terms, dHdx_Max_Parts, dHdx_Entry_Offset, dHdx_Row_Offset >( tx, s_track, r_cgesvA, d_Hx_idx, s_phc_coeffs_Hx );
        eval_Homotopy< dHdt_Max_Terms, dHdt_Max_Parts, dHdt_Row_Offset >( tx, s_track, r_cgesvB, d_Ht_idx, s_phc_coeffs_Hx );

        //> G-N corrector first solve
        cgesv_batched_small_device< Num_Of_Vars >( tx, r_cgesvA, sipiv, r_cgesvB, sB, sx, dsx, rowid, linfo );
        magmablas_syncwarp();

        //> correct the sols
        s_track[tx] -= sB[tx];
        magmablas_syncwarp();

        r_sqrt_sols = MAGMA_C_REAL(sB[tx])*MAGMA_C_REAL(sB[tx]) + MAGMA_C_IMAG(sB[tx])*MAGMA_C_IMAG(sB[tx]);
        r_sqrt_corr = MAGMA_C_REAL(s_track[tx])*MAGMA_C_REAL(s_track[tx]) + MAGMA_C_IMAG(s_track[tx])*MAGMA_C_IMAG(s_track[tx]);
        magmablas_syncwarp();

        for (int offset = WARP_SIZE/2; offset > 0; offset /= 2 ) {
            r_sqrt_sols += __shfl_down_sync(__activemask(), r_sqrt_sols, offset);
            r_sqrt_corr += __shfl_down_sync(__activemask(), r_sqrt_corr, offset);
        }

        if ( tx == 0 ) {
            r_isSuccessful = r_sqrt_sols < 0.000001 * r_sqrt_corr;
            r_isInfFail = (r_sqrt_corr > 1e14) ? (true) : (false);
        }
        //> Broadcast the values of r_isSuccessful and r_isInfFail from thread 0 to all the rest of the threads
        r_isSuccessful = __shfl_sync(__activemask(), r_isSuccessful, 0);
        r_isInfFail = __shfl_sync(__activemask(), r_isInfFail, 0);

        if (r_isInfFail) break;
        if (r_isSuccessful) break;
      }

      if ( r_isInfFail ) break;

      // ===================================================================
      //> Decide Track Changes
      // ===================================================================
      if (!r_isSuccessful) {
        delta_t *= 0.5;
        //> should be the last successful tracked sols
        s_track[tx] = s_track_last_success[tx];
        s_sols[tx] = s_track_last_success[tx];
        if (tx == 0) sipiv[Num_Of_Vars] = 0;
        magmablas_syncwarp();
        t0 = t_step;
      }
      else {
        if (tx == 0) sipiv[Num_Of_Vars]++;
        s_track_last_success[tx] = s_track[tx];
        s_sols[tx] = s_track[tx];
        magmablas_syncwarp();
        if (sipiv[Num_Of_Vars] >= HC_delta_t_incremental_steps) {
          if (tx == 0) sipiv[Num_Of_Vars] = 0;
          delta_t *= 2;
        }
        magmablas_syncwarp();
      }
    }
    else {
      break;
    }
  }

  d_track[tx] = s_track[tx];
  if (tx == 0) {
    d_is_GPU_HC_Sol_Converge[ batchid ] = (t0 >= 1.0 || (1.0-t0 <= 0.0000001)) ? (1) : (0);
    d_is_GPU_HC_Sol_Infinity[ batchid ] = (r_isInfFail) ? (1) : (0);
  }

#if GPU_DEBUG
  d_Debug_Purpose[ batchid ] = (t0 >= 1.0 || (1.0-t0 <= 0.0000001)) ? MAGMA_C_MAKE(1.0, 0.0) : MAGMA_C_MAKE(t0, delta_t);
#endif
}

real_Double_t
kernel_HC_Solver_5pt_rel_pos_full_form_quat(
  magma_queue_t         my_queue,
  int                   HC_max_steps, 
  int                   HC_max_correction_steps, 
  int                   HC_delta_t_incremental_steps,
  magmaFloatComplex**   d_startSols_array, 
  magmaFloatComplex**   d_Track_array,
  magma_int_t*          d_Hx_idx_array,
  magma_int_t*          d_Ht_idx_array, 
  magmaFloatComplex_ptr d_phc_coeffs_Hx, 
  magmaFloatComplex_ptr d_phc_coeffs_Ht,
  bool*                 d_is_GPU_HC_Sol_Converge,
  bool*                 d_is_GPU_HC_Sol_Infinity,
  magmaFloatComplex*    d_Debug_Purpose
)
{
  //> Hard-coded for each problem
  const int num_of_vars               = 16;
  const int num_of_tracks             = 40;
  const int num_of_coeffs_from_params = 20;
  const int dHdx_Max_Terms            = 7;
  const int dHdx_Max_Parts            = 4;
  const int dHdt_Max_Terms            = 9;
  const int dHdt_Max_Parts            = 5;
  const int max_order_of_t            = 1;

  const int dHdx_Entry_Offset = dHdx_Max_Terms * dHdx_Max_Parts;
  const int dHdx_Row_Offset   = num_of_vars * dHdx_Entry_Offset;
  const int dHdt_Row_Offset   = dHdt_Max_Terms * dHdt_Max_Parts;

  real_Double_t gpu_time;
  dim3 threads(num_of_vars, 1, 1);
  dim3 grid(num_of_tracks, 1, 1);
  cudaError_t e = cudaErrorInvalidValue;

  //> Constant values for evaluating the Jacobians, passed as template
  const unsigned Full_Parallel_Offset                 = (num_of_coeffs_from_params+1)/(num_of_vars);
  const unsigned Partial_Parallel_Thread_Offset       = (num_of_coeffs_from_params+1) - (num_of_vars)*(Full_Parallel_Offset);
  const unsigned Partial_Parallel_Index_Offset        = (num_of_vars)*(Full_Parallel_Offset);
  const unsigned Max_Order_of_t_Plus_One              = max_order_of_t + 1;
  const unsigned Partial_Parallel_Index_Offset_for_Hx = (num_of_vars-1)*(Max_Order_of_t_Plus_One) + (max_order_of_t) + (Full_Parallel_Offset-1)*(Max_Order_of_t_Plus_One)*(num_of_vars) + 1;
  const unsigned Partial_Parallel_Index_Offset_for_Ht = (num_of_vars-1)*(max_order_of_t) + (max_order_of_t-1) + (Full_Parallel_Offset-1)*(max_order_of_t)*(num_of_vars) + 1;

  //> declare shared memory
  magma_int_t shmem  = 0;
  shmem += (num_of_vars+1)               * sizeof(magmaFloatComplex);   // startSols
  shmem += (num_of_vars+1)               * sizeof(magmaFloatComplex);   // track
  shmem += (num_of_vars+1)               * sizeof(magmaFloatComplex);   // track_pred_init
  shmem += (num_of_coeffs_from_params+1) * sizeof(magmaFloatComplex);   // s_phc_coeffs_Hx
  shmem += (num_of_coeffs_from_params+1) * sizeof(magmaFloatComplex);   // s_phc_coeffs_Ht
  shmem += num_of_vars                   * sizeof(magmaFloatComplex);   // sB
  shmem += num_of_vars                   * sizeof(magmaFloatComplex);   // sx
  shmem += num_of_vars                   * sizeof(float);               // dsx
  shmem += num_of_vars                   * sizeof(int);                 // pivot
  shmem += 1 * sizeof(int);                                             // predictor_success counter
  shmem += 1 * sizeof(int);                                             // Loopy Runge-Kutta coefficients
  shmem += 1 * sizeof(float);                                           // Loopy Runge-Kutta delta t

  //> Get max. dynamic shared memory on the GPU
  int nthreads_max, shmem_max = 0;
  cudacheck( cudaDeviceGetAttribute(&nthreads_max, cudaDevAttrMaxThreadsPerBlock, 0) );
#if CUDA_VERSION >= 9000
  cudacheck( cudaDeviceGetAttribute (&shmem_max, cudaDevAttrMaxSharedMemoryPerBlockOptin, 0) );
  if (shmem <= shmem_max) {
    cudacheck( cudaFuncSetAttribute(HC_solver_5pt_rel_pos_full_form_quat \
                                    < num_of_vars, num_of_coeffs_from_params, max_order_of_t, \
                                      dHdx_Max_Terms, dHdx_Max_Parts, dHdx_Entry_Offset, dHdx_Row_Offset, \
                                      dHdt_Max_Terms, dHdt_Max_Parts, dHdt_Row_Offset, \
                                      Full_Parallel_Offset, \
                                      Partial_Parallel_Thread_Offset, \
                                      Partial_Parallel_Index_Offset, \
                                      Max_Order_of_t_Plus_One, \
                                      Partial_Parallel_Index_Offset_for_Hx, \
                                      Partial_Parallel_Index_Offset_for_Ht>, \
                                    cudaFuncAttributeMaxDynamicSharedMemorySize, shmem) );
  }
#else
  cudacheck( cudaDeviceGetAttribute (&shmem_max, cudaDevAttrMaxSharedMemoryPerBlock, 0) );
#endif

  //> Message of overuse shared memory
  if ( shmem > shmem_max ) printf("Error: kernel %s requires too many threads or too much shared memory\n", __func__);

  void *kernel_args[] = { &HC_max_steps, &HC_max_correction_steps, &HC_delta_t_incremental_steps, \
                          &d_startSols_array, &d_Track_array, \
                          &d_Hx_idx_array, &d_Ht_idx_array, \
                          &d_phc_coeffs_Hx, &d_phc_coeffs_Ht, \
                          &d_is_GPU_HC_Sol_Converge, &d_is_GPU_HC_Sol_Infinity, \
                          &d_Debug_Purpose };

  gpu_time = magma_sync_wtime( my_queue );

  // cudacheck( cudaEventRecord(start) );
  e = cudaLaunchKernel((void*)HC_solver_5pt_rel_pos_full_form_quat \
                        < num_of_vars, num_of_coeffs_from_params, max_order_of_t, \
                          dHdx_Max_Terms, dHdx_Max_Parts, dHdx_Entry_Offset, dHdx_Row_Offset, \
                          dHdt_Max_Terms, dHdt_Max_Parts, dHdt_Row_Offset, \
                          Full_Parallel_Offset, \
                          Partial_Parallel_Thread_Offset, \
                          Partial_Parallel_Index_Offset, \
                          Max_Order_of_t_Plus_One, \
                          Partial_Parallel_Index_Offset_for_Hx, \
                          Partial_Parallel_Index_Offset_for_Ht>, \
                        grid, threads, kernel_args, shmem, my_queue->cuda_stream());

  gpu_time = magma_sync_wtime( my_queue ) - gpu_time;
  if( e != cudaSuccess ) printf("cudaLaunchKernel of HC_solver_5pt_rel_pos_full_form_quat is not successful!\n");

  return gpu_time;
}

#endif
