#ifndef kernel_HC_Solver_5pt_rel_pos_geo_form_quat_cu
#define kernel_HC_Solver_5pt_rel_pos_geo_form_quat_cu
// =======================================================================================
// GPU homotopy continuation solver for 5-point relative pose problem (Geometric Form)
//
// Modifications
//    Chiang-Heng Chien  22-10-31:   Initially created
//    Chiang-Heng Chien  24-01-04:   Add macro definitions for computing coefficients from parameter homotopy
//
//> (c) LEMS, Brown University
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// =======================================================================================
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>

#include <cuda.h>
#include <cuda_runtime.h>

//> MAGMA
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

#include "magmaHC-kernels.hpp"

#include "../definitions.hpp"
#include "../gpu-idx-evals/dev-eval-indxing-5pt_rel_pos_geo_form_quat.cuh"
#include "../dev-cgesv-batched-small.cuh"
#include "../dev-get-new-data.cuh"

//namespace GPU_Device {

  template< unsigned Full_Parallel_Offset, \
            unsigned Partial_Parallel_Thread_Offset, \
            unsigned Partial_Parallel_Index_Offset, \
            unsigned Max_Order_of_t_Plus_One, \
            unsigned Partial_Parallel_Index_Offset_Hx, \
            unsigned Partial_Parallel_Index_Offset_Ht >
  __global__ void
  homotopy_continuation_solver_5pt_rel_pos_geo_form_quat(
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

    magmaFloatComplex* d_startSols = d_startSols_array[batchid];
    magmaFloatComplex* d_track     = d_Track_array[batchid];
    const magmaFloatComplex* __restrict__ d_const_phc_coeffs_Hx = d_phc_coeffs_Hx;
    const magmaFloatComplex* __restrict__ d_const_phc_coeffs_Ht = d_phc_coeffs_Ht;
    const int* __restrict__ d_Hx_idx = d_Hx_indices;
    const int* __restrict__ d_Ht_idx = d_Ht_indices;

    //> Registers
    magmaFloatComplex r_cgesvA[NUM_OF_VARS] = {MAGMA_C_ZERO};
    magmaFloatComplex r_cgesvB = MAGMA_C_ZERO;
    int linfo = 0, rowid = tx;
    float t0 = 0.0, t_step = 0.0, delta_t = 0.05;
    bool inf_failed = 0;
    bool end_zone = 0;

    //> Shared Memory
    magmaFloatComplex *s_sols               = (magmaFloatComplex*)(zdata);
    magmaFloatComplex *s_track              = s_sols + (NUM_OF_VARS+1);
    magmaFloatComplex *s_track_last_success = s_track + (NUM_OF_VARS+1);
    magmaFloatComplex *sB                   = s_track_last_success + (NUM_OF_VARS+1);
    magmaFloatComplex *sx                   = sB + NUM_OF_VARS;
    magmaFloatComplex *s_phc_coeffs_Hx      = sx + NUM_OF_VARS;
    magmaFloatComplex *s_phc_coeffs_Ht      = s_phc_coeffs_Hx + (NUM_OF_COEFFS_FROM_PARAMS+1);
    float* dsx                              = (float*)(s_phc_coeffs_Ht + (NUM_OF_COEFFS_FROM_PARAMS+1));
    int* sipiv                              = (int*)(dsx + NUM_OF_VARS);
    float *s_sqrt_sols                      = (float*)(sipiv + NUM_OF_VARS);
    float *s_sqrt_corr                      = s_sqrt_sols + NUM_OF_VARS;
    float *s_norm                           = s_sqrt_corr + NUM_OF_VARS;
    bool s_isSuccessful                     = (bool)(s_norm + 2);
    int s_pred_success_count                = (int)(s_isSuccessful + 1);

    s_sols[tx] = d_startSols[tx];
    s_track[tx] = d_track[tx];
    s_track_last_success[tx] = s_track[tx];
    s_sqrt_sols[tx] = 0;
    s_sqrt_corr[tx] = 0;
    s_isSuccessful = 0;
    s_pred_success_count = 0;
    if (tx == 0) {
      s_sols[NUM_OF_VARS] = MAGMA_C_MAKE(1.0, 0.0);
      s_track[NUM_OF_VARS] = MAGMA_C_MAKE(1.0, 0.0);
      s_track_last_success[NUM_OF_VARS] = MAGMA_C_MAKE(1.0, 0.0);
    }
    __syncthreads();

    float one_half_delta_t;   // -- 1/2 \Delta t --

    #pragma unroll
    for (int step = 0; step <= HC_MAX_STEPS; step++) {
      if (t0 < 1.0 && (1.0-t0 > 0.0000001)) {

        // ===================================================================
        // -- Decide delta t at end zone --
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
        // -- Runge-Kutta Predictor --
        // ===================================================================
        // -- get HxHt for k1 --
        eval_parameter_homotopy<Full_Parallel_Offset, Partial_Parallel_Thread_Offset, Partial_Parallel_Index_Offset, \
                                Max_Order_of_t_Plus_One, Partial_Parallel_Index_Offset_Hx, Partial_Parallel_Index_Offset_Ht> \
                                ( tx, t0, s_phc_coeffs_Hx, s_phc_coeffs_Ht, d_const_phc_coeffs_Hx, d_const_phc_coeffs_Ht );
        eval_Jacobian_Hx< HX_MAXIMAL_TERMS*HX_MAXIMAL_PARTS, NUM_OF_VARS*HX_MAXIMAL_TERMS*HX_MAXIMAL_PARTS>( tx, s_track, r_cgesvA, d_Hx_idx, s_phc_coeffs_Hx );
        eval_Jacobian_Ht< HT_MAXIMAL_TERMS*HT_MAXIMAL_PARTS >( tx, s_track, r_cgesvB, d_Ht_idx, s_phc_coeffs_Ht );

        // -- solve k1 --
        cgesv_batched_small_device< NUM_OF_VARS >( tx, r_cgesvA, sipiv, r_cgesvB, sB, sx, dsx, rowid, linfo );
        magmablas_syncwarp();

        // -- compute x for the creation of HxHt for k2 --
        create_x_for_k2( tx, t0, delta_t, one_half_delta_t, s_sols, s_track, sB );
        magmablas_syncwarp();

        // -- get HxHt for k2 --
        eval_parameter_homotopy<Full_Parallel_Offset, Partial_Parallel_Thread_Offset, Partial_Parallel_Index_Offset, \
                                Max_Order_of_t_Plus_One, Partial_Parallel_Index_Offset_Hx, Partial_Parallel_Index_Offset_Ht> \
                                ( tx, t0, s_phc_coeffs_Hx, s_phc_coeffs_Ht, d_const_phc_coeffs_Hx, d_const_phc_coeffs_Ht );
        eval_Jacobian_Hx< HX_MAXIMAL_TERMS*HX_MAXIMAL_PARTS, NUM_OF_VARS*HX_MAXIMAL_TERMS*HX_MAXIMAL_PARTS>( tx, s_track, r_cgesvA, d_Hx_idx, s_phc_coeffs_Hx );
        eval_Jacobian_Ht< HT_MAXIMAL_TERMS*HT_MAXIMAL_PARTS >( tx, s_track, r_cgesvB, d_Ht_idx, s_phc_coeffs_Ht );

        // -- solve k2 --
        cgesv_batched_small_device< NUM_OF_VARS >( tx, r_cgesvA, sipiv, r_cgesvB, sB, sx, dsx, rowid, linfo );
        magmablas_syncwarp();

        // -- compute x for the generation of HxHt for k3 --
        create_x_for_k3( tx, delta_t, one_half_delta_t, s_sols, s_track, s_track_last_success, sB );
        magmablas_syncwarp();

        // -- get HxHt for k3 --
        eval_Jacobian_Hx< HX_MAXIMAL_TERMS*HX_MAXIMAL_PARTS, NUM_OF_VARS*HX_MAXIMAL_TERMS*HX_MAXIMAL_PARTS>( tx, s_track, r_cgesvA, d_Hx_idx, s_phc_coeffs_Hx );
        eval_Jacobian_Ht< HT_MAXIMAL_TERMS*HT_MAXIMAL_PARTS >( tx, s_track, r_cgesvB, d_Ht_idx, s_phc_coeffs_Ht );

        // -- solve k3 --
        cgesv_batched_small_device< NUM_OF_VARS >( tx, r_cgesvA, sipiv, r_cgesvB, sB, sx, dsx, rowid, linfo );
        magmablas_syncwarp();

        // -- compute x for the generation of HxHt for k4 --
        create_x_for_k4( tx, t0, delta_t, one_half_delta_t, s_sols, s_track, s_track_last_success, sB );
        magmablas_syncwarp();

        // -- get HxHt for k4 --
        eval_parameter_homotopy<Full_Parallel_Offset, Partial_Parallel_Thread_Offset, Partial_Parallel_Index_Offset, \
                                Max_Order_of_t_Plus_One, Partial_Parallel_Index_Offset_Hx, Partial_Parallel_Index_Offset_Ht> \
                                ( tx, t0, s_phc_coeffs_Hx, s_phc_coeffs_Ht, d_const_phc_coeffs_Hx, d_const_phc_coeffs_Ht );
        eval_Jacobian_Hx< HX_MAXIMAL_TERMS*HX_MAXIMAL_PARTS, NUM_OF_VARS*HX_MAXIMAL_TERMS*HX_MAXIMAL_PARTS>( tx, s_track, r_cgesvA, d_Hx_idx, s_phc_coeffs_Hx );
        eval_Jacobian_Ht< HT_MAXIMAL_TERMS*HT_MAXIMAL_PARTS >( tx, s_track, r_cgesvB, d_Ht_idx, s_phc_coeffs_Ht );

        // -- solve k4 --
        cgesv_batched_small_device< NUM_OF_VARS >( tx, r_cgesvA, sipiv, r_cgesvB, sB, sx, dsx, rowid, linfo );
        magmablas_syncwarp();

        // -- make prediction --
        s_sols[tx] += sB[tx] * delta_t * 1.0/6.0;
        s_track[tx] = s_sols[tx];
        __syncthreads();

        // ===================================================================
        // -- Gauss-Newton Corrector --
        // ===================================================================
        //#pragma unroll
        for(int i = 0; i < HC_MAX_CORRECTION_STEPS; i++) {

          eval_Jacobian_Hx< HX_MAXIMAL_TERMS*HX_MAXIMAL_PARTS, NUM_OF_VARS*HX_MAXIMAL_TERMS*HX_MAXIMAL_PARTS>( tx, s_track, r_cgesvA, d_Hx_idx, s_phc_coeffs_Hx );
          eval_Homotopy< HT_MAXIMAL_TERMS*HT_MAXIMAL_PARTS >( tx, s_track, r_cgesvB, d_Ht_idx, s_phc_coeffs_Hx );

          // -- G-NUM_OF_VARS corrector first solve --
          cgesv_batched_small_device< NUM_OF_VARS >( tx, r_cgesvA, sipiv, r_cgesvB, sB, sx, dsx, rowid, linfo );
          magmablas_syncwarp();

          // -- correct the sols --
          s_track[tx] -= sB[tx];
          __syncthreads();

          // -- compute the norms; norm[0] is norm(sB), norm[1] is norm(sol) --
          compute_norm2( tx, sB, s_track, s_sqrt_sols, s_sqrt_corr, s_norm );
          __syncthreads();

          s_isSuccessful = s_norm[0] < 0.000001 * s_norm[1];
          __syncthreads();

          if (s_isSuccessful)
	           break;
        }

        //> stop if the values of the solution is too large
        if ((s_norm[1] > 1e14) && (t0 < 1.0) && (1.0-t0 > 0.001)) {
          inf_failed = 1;
          break;
        }

        // ===================================================================
        // -- Decide Track Changes --
        // ===================================================================
        if (!s_isSuccessful) {
          s_pred_success_count = 0;
          delta_t *= 0.5;
          // -- should be the last successful tracked sols --
          s_track[tx] = s_track_last_success[tx];
          s_sols[tx] = s_track_last_success[tx];
          __syncthreads();
          t0 = t_step;
        }
        else {
          s_track_last_success[tx] = s_track[tx];
          s_sols[tx] = s_track[tx];
          __syncthreads();
          s_pred_success_count++;
          if (s_pred_success_count >= HC_NUM_OF_STEPS_TO_INCREASE_DELTA_T) {
            s_pred_success_count = 0;
            delta_t *= 2;
          }
        }
      }
      else {
        break;
      }
    }

    d_track[tx] = s_track[tx];
    if (tx == 0) {
      d_is_GPU_HC_Sol_Converge[ batchid ] = (t0 >= 1.0 || (1.0-t0 <= 0.0000001)) ? (1) : (0);
      d_is_GPU_HC_Sol_Infinity[ batchid ] = (inf_failed) ? (1) : (0);
    }

#if GPU_DEBUG
    d_Debug_Purpose[ batchid ] = (t0 >= 1.0 || (1.0-t0 <= 0.0000001)) ? MAGMA_C_MAKE(1.0, 0.0) : MAGMA_C_MAKE(t0, delta_t);
#endif
  }

  real_Double_t
  kernel_HC_Solver_5pt_rel_pos_geo_form_quat(
    magma_queue_t my_queue, \
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array, \
    magma_int_t* d_Hx_idx_array,           magma_int_t* d_Ht_idx_array, \
    magmaFloatComplex_ptr d_phc_coeffs_Hx, magmaFloatComplex_ptr d_phc_coeffs_Ht, \
    bool* d_is_GPU_HC_Sol_Converge,        bool* d_is_GPU_HC_Sol_Infinity, \
    magmaFloatComplex* d_Debug_Purpose
  )
  {
    real_Double_t gpu_time;
    dim3 threads(NUM_OF_VARS, 1, 1);
    dim3 grid(NUM_OF_TRACKS, 1, 1);
    cudaError_t e = cudaErrorInvalidValue;

    //> Constant values for evaluating the Jacobians, passed as template
    const unsigned Full_Parallel_Offset                 = (NUM_OF_COEFFS_FROM_PARAMS+1)/(NUM_OF_VARS);
    const unsigned Partial_Parallel_Thread_Offset       = (NUM_OF_COEFFS_FROM_PARAMS+1) - (NUM_OF_VARS)*(Full_Parallel_Offset);
    const unsigned Partial_Parallel_Index_Offset        = (NUM_OF_VARS)*(Full_Parallel_Offset);
    const unsigned Max_Order_of_t_Plus_One              = MAX_ORDER_OF_T + 1;
    const unsigned Partial_Parallel_Index_Offset_for_Hx = (NUM_OF_VARS-1)*(Max_Order_of_t_Plus_One) + (MAX_ORDER_OF_T) + (Full_Parallel_Offset-1)*(Max_Order_of_t_Plus_One)*(NUM_OF_VARS) + 1;
    const unsigned Partial_Parallel_Index_Offset_for_Ht = (NUM_OF_VARS-1)*(MAX_ORDER_OF_T) + (MAX_ORDER_OF_T-1) + (Full_Parallel_Offset-1)*(MAX_ORDER_OF_T)*(NUM_OF_VARS) + 1;

    magma_int_t shmem  = 0;
    shmem += (NUM_OF_VARS+1) * sizeof(magmaFloatComplex);       // startSols
    shmem += (NUM_OF_VARS+1) * sizeof(magmaFloatComplex);       // track
    shmem += (NUM_OF_VARS+1) * sizeof(magmaFloatComplex);       // track_pred_init

    shmem += (NUM_OF_COEFFS_FROM_PARAMS+1) * sizeof(magmaFloatComplex);  // -- s_phc_coeffs_Hx --
    shmem += (NUM_OF_COEFFS_FROM_PARAMS+1) * sizeof(magmaFloatComplex);  // -- s_phc_coeffs_Ht --

    shmem += NUM_OF_VARS * sizeof(magmaFloatComplex); // sB
    shmem += NUM_OF_VARS * sizeof(magmaFloatComplex); // sx
    shmem += NUM_OF_VARS * sizeof(float);            // dsx
    shmem += NUM_OF_VARS * sizeof(int);               // pivot
    shmem += NUM_OF_VARS * sizeof(float);             // s_sqrt for sol norm-2 in G-NUM_OF_VARS corrector
    shmem += NUM_OF_VARS * sizeof(float);             // s_sqrt for corr norm-2 in G-NUM_OF_VARS corrector
    shmem += 2 * sizeof(float);             // s_norm for norm-2 in G-NUM_OF_VARS corrector
    shmem += 1 * sizeof(bool);              // is_successful 
    shmem += 1 * sizeof(int);               // predictor_success counter

    void *kernel_args[] = { &d_startSols_array, &d_Track_array, \
                            &d_Hx_idx_array, &d_Ht_idx_array, \
                            &d_phc_coeffs_Hx, &d_phc_coeffs_Ht, \
                            &d_is_GPU_HC_Sol_Converge, &d_is_GPU_HC_Sol_Infinity, \
                            &d_Debug_Purpose };

    gpu_time = magma_sync_wtime( my_queue );
    
    e = cudaLaunchKernel((void*)homotopy_continuation_solver_5pt_rel_pos_geo_form_quat \
                          <Full_Parallel_Offset, \
                           Partial_Parallel_Thread_Offset, \
                           Partial_Parallel_Index_Offset, \
                           Max_Order_of_t_Plus_One, \
                           Partial_Parallel_Index_Offset_for_Hx, \
                           Partial_Parallel_Index_Offset_for_Ht>, \
                          grid, threads, kernel_args, shmem, my_queue->cuda_stream());

    gpu_time = magma_sync_wtime( my_queue ) - gpu_time;
    if( e != cudaSuccess ) printf("cudaLaunchKernel of homotopy_continuation_solver_5pt_rel_pos_geo_form_quat is not successful!\n");

    return gpu_time;
  }

//}

#endif
