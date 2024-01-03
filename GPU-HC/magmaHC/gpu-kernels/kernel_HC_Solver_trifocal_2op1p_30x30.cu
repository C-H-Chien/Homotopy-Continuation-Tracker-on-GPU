#ifndef kernel_HC_Solver_trifocal_2op1p_30x30_cu
#define kernel_HC_Solver_trifocal_2op1p_30x30_cu
// ===========================================================================================
// GPU homotopy continuation solver for the trifocal 2op1p 30x30 problem
// Version 2: Direct evaluation of parameter homotopy. The coefficient 
//            part of each polynomial is not expanded to an uni-variable 
//            polynomial. Rather, depending on the order of t, the parameter 
//            homotopy formulation is explicitly hard-coded such that we do not 
//            need to compute coefficients from parameters first then use 
//            index-based to evaluate the Jacobians. The required amount of data 
//            to be stored in a kernel in this method is reduced which expects 
//            to speedup over the first version.
//
// Major Modifications
//    Chiang-Heng Chien  22-10-03:   Edited from the first version 
//                                   (kernel_HC_Solver_trifocal_2op1p_30.cu)
//    Chiang-Heng Chien  23-12-28:   Add macros and circular arc homotopy for gamma-trick
//
// ============================================================================================
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>

#include <cuda.h>
#include <cuda_runtime.h>

//> MAGMA
#include "flops.h"
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

//> device functions
#include "../gpu-idx-evals/dev-eval-indxing-trifocal_2op1p_30x30.cuh"
#include "../dev-cgesv-batched-small.cuh"
#include "../dev-get-new-data.cuh"

//namespace GPU_Device {

  __global__ void
  HC_solver_trifocal_2op1p_30x30(
    magmaFloatComplex** d_startSols_array,
    magmaFloatComplex** d_Track_array,
    magmaFloatComplex*  d_startParams,
    magmaFloatComplex*  d_targetParams,
    magmaFloatComplex*  d_diffParams,
    const magma_int_t* __restrict__ d_Hx_indices,
    const magma_int_t* __restrict__ d_Ht_indices,
    bool* d_is_GPU_HC_Sol_Converge,
    bool* d_is_GPU_HC_Sol_Infinity,
    magmaFloatComplex* d_Debug_Purpose
  )
  {
    extern __shared__ magmaFloatComplex zdata[];
    const int tx = threadIdx.x;
    const int batchid = blockIdx.x ;

    //> define pointers to the arrays
    magmaFloatComplex* d_startSols    = d_startSols_array[batchid];
    magmaFloatComplex* d_track        = d_Track_array[batchid];

    //> declarations of registers
    magmaFloatComplex r_cgesvA[NUM_OF_VARS] = {MAGMA_C_ZERO};
    magmaFloatComplex r_cgesvB              = MAGMA_C_ZERO;
    int linfo = 0, rowid = tx;
    float t0 = 0.0, t_step = 0.0, delta_t = 0.01;
    bool end_zone = 0;
    int hc_step = 0;
    bool inf_failed = 0;

    //> declarations of shared memories
    magmaFloatComplex *s_startParams        = (magmaFloatComplex*)(zdata);
    magmaFloatComplex *s_targetParams       = s_startParams + (NUM_OF_PARAMS + 1);
    magmaFloatComplex *s_diffParams         = s_targetParams + (NUM_OF_PARAMS + 1);
    magmaFloatComplex *s_param_homotopy     = s_diffParams + (NUM_OF_PARAMS + 1);
    magmaFloatComplex *s_sols               = s_param_homotopy + (NUM_OF_PARAMS + 1);
    magmaFloatComplex *s_track              = s_sols + (NUM_OF_VARS+1);
    magmaFloatComplex *s_track_last_success = s_track + (NUM_OF_VARS+1);
    magmaFloatComplex *sB                   = s_track_last_success + (NUM_OF_VARS+1);
    magmaFloatComplex *sx                   = sB + NUM_OF_VARS;
    float *dsx                              = (float*)(sx + NUM_OF_VARS);
    float *s_sqrt_sols                      = dsx + NUM_OF_VARS;
    float *s_sqrt_corr                      = s_sqrt_sols + NUM_OF_VARS;
    float *s_norm                           = s_sqrt_corr + NUM_OF_VARS;
    int   *sipiv                            = (int*)(s_norm + 2);
    bool   s_isSuccessful                   = (bool)(sipiv + NUM_OF_VARS);
    int    s_pred_success_count             = (int)(s_isSuccessful + 1);

    //> read data from global memory to shared memories or do initializations
    s_sols[tx]               = d_startSols[tx];
    s_track[tx]              = d_track[tx];
    s_track_last_success[tx] = s_track[tx];
    s_sqrt_sols[tx]          = 0;
    s_sqrt_corr[tx]          = 0;
    s_isSuccessful           = 0;
    s_pred_success_count     = 0;

    //> start and target parameters
    s_startParams[tx]  = d_startParams[tx];
    s_targetParams[tx] = d_targetParams[tx];
    s_diffParams[tx]   = d_diffParams[tx];

    if (tx == 0) {
      //> the rest of the start and target parameters
      #pragma unroll
      for(int i = NUM_OF_VARS; i <= NUM_OF_PARAMS; i++) {
        s_startParams[i]  = d_startParams[i];
        s_targetParams[i] = d_targetParams[i];
        s_diffParams[i]   = d_diffParams[i];
      }
      s_sols[NUM_OF_VARS]               = MAGMA_C_MAKE(1.0, 0.0);
      s_track[NUM_OF_VARS]              = MAGMA_C_MAKE(1.0, 0.0);
      s_track_last_success[NUM_OF_VARS] = MAGMA_C_MAKE(1.0, 0.0);
      s_param_homotopy[NUM_OF_PARAMS]   = MAGMA_C_ONE;
    }
    __syncthreads();

    //> 1/2 \Delta t
    float one_half_delta_t;

    //#pragma unroll
    for (int step = 0; step <= HC_MAX_STEPS; step++) {
      if (t0 < 1.0 && (1.0-t0 > 0.0000001)) {

        // ===================================================================
        // Decide delta t at end zone
        // ===================================================================
        if (!end_zone && fabs(1 - t0) <= (0.0500001)) {
          end_zone = true;
        }

        if (end_zone) {
          if (delta_t > fabs(1 - t0))
            delta_t = fabs(1 - t0);
        }
        else if (delta_t > fabs(0.95 - t0)) {
          delta_t = fabs(0.95 - t0);
        }

        t_step = t0;
        one_half_delta_t = 0.5 * delta_t;
        // ===================================================================
        // Prediction: 4-th order Runge-Kutta method
        // ===================================================================
        //> get HxHt for k1
        compute_param_homotopy( tx, t0, s_param_homotopy, s_startParams, s_targetParams );
        eval_Jacobian_Hx< HX_MAXIMAL_TERMS*HX_MAXIMAL_PARTS >( tx, t0, r_cgesvA, s_track, s_startParams, s_targetParams, s_param_homotopy, d_Hx_indices );
        eval_Jacobian_Ht( tx, t0, r_cgesvB, s_track, s_startParams, s_targetParams, s_param_homotopy, d_Ht_indices, s_diffParams );

        //> solve k1
        cgesv_batched_small_device< NUM_OF_VARS >( tx, r_cgesvA, sipiv, r_cgesvB, sB, sx, dsx, rowid, linfo );
        magmablas_syncwarp();

        //> compute x for the creation of HxHt for k2
        create_x_for_k2( tx, t0, delta_t, one_half_delta_t, s_sols, s_track, sB );
        magmablas_syncwarp();

        //> get HxHt for k2
        compute_param_homotopy( tx, t0, s_param_homotopy, s_startParams, s_targetParams );
        eval_Jacobian_Hx< HX_MAXIMAL_TERMS*HX_MAXIMAL_PARTS >( tx, t0, r_cgesvA, s_track, s_startParams, s_targetParams, s_param_homotopy, d_Hx_indices );
        eval_Jacobian_Ht( tx, t0, r_cgesvB, s_track, s_startParams, s_targetParams, s_param_homotopy, d_Ht_indices, s_diffParams );

        //> solve k2
        cgesv_batched_small_device< NUM_OF_VARS >( tx, r_cgesvA, sipiv, r_cgesvB, sB, sx, dsx, rowid, linfo );
        magmablas_syncwarp();

        //> compute x for the generation of HxHt for k3
        create_x_for_k3( tx, delta_t, one_half_delta_t, s_sols, s_track, s_track_last_success, sB );
        magmablas_syncwarp();

        //> get HxHt for k3
        //compute_param_homotopy<NUM_OF_VARS>( tx, t0, s_param_homotopy, s_start_params, s_target_params );
        eval_Jacobian_Hx< HX_MAXIMAL_TERMS*HX_MAXIMAL_PARTS >( tx, t0, r_cgesvA, s_track, s_startParams, s_targetParams, s_param_homotopy, d_Hx_indices );
        eval_Jacobian_Ht( tx, t0, r_cgesvB, s_track, s_startParams, s_targetParams, s_param_homotopy, d_Ht_indices, s_diffParams );

        //> solve k3
        cgesv_batched_small_device< NUM_OF_VARS >( tx, r_cgesvA, sipiv, r_cgesvB, sB, sx, dsx, rowid, linfo );
        magmablas_syncwarp();

        //> compute x for the generation of HxHt for k4
        create_x_for_k4( tx, t0, delta_t, one_half_delta_t, s_sols, s_track, s_track_last_success, sB );
        magmablas_syncwarp();

        //> get HxHt for k4
        compute_param_homotopy( tx, t0, s_param_homotopy, s_startParams, s_targetParams );
        eval_Jacobian_Hx< HX_MAXIMAL_TERMS*HX_MAXIMAL_PARTS >( tx, t0, r_cgesvA, s_track, s_startParams, s_targetParams, s_param_homotopy, d_Hx_indices );
        eval_Jacobian_Ht( tx, t0, r_cgesvB, s_track, s_startParams, s_targetParams, s_param_homotopy, d_Ht_indices, s_diffParams );

        //> solve k4
        cgesv_batched_small_device< NUM_OF_VARS >( tx, r_cgesvA, sipiv, r_cgesvB, sB, sx, dsx, rowid, linfo );
        magmablas_syncwarp();

        //> make prediction
        s_sols[tx] += sB[tx] * delta_t * 1.0/6.0;
        s_track[tx] = s_sols[tx];
        __syncthreads();

        // ===================================================================
        //> Gauss-Newton Corrector
        // ===================================================================
        //#pragma unroll
        for(int i = 0; i < HC_MAX_CORRECTION_STEPS; i++) {

          //> evaluate the Jacobian Hx
          eval_Jacobian_Hx<HX_MAXIMAL_TERMS*HX_MAXIMAL_PARTS>( tx, t0, r_cgesvA, s_track, s_startParams, s_targetParams, s_param_homotopy, d_Hx_indices );

          //> evaluate the parameter homotopy
          eval_Parameter_Homotopy( tx, t0, r_cgesvB, s_track, s_startParams, s_targetParams, s_param_homotopy, d_Ht_indices );

          //> G-NUM_OF_VARS corrector first solve
          cgesv_batched_small_device< NUM_OF_VARS >( tx, r_cgesvA, sipiv, r_cgesvB, sB, sx, dsx, rowid, linfo );
          magmablas_syncwarp();

          //> correct the sols
          s_track[tx] -= sB[tx];
          __syncthreads();

          //> compute the norms; norm[0] is norm(sB), norm[1] is norm(sol)
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
        // Decide Track Changes
        // ===================================================================
        if (!s_isSuccessful) {
          s_pred_success_count = 0;
          delta_t *= 0.5;
          //> should be the last successful tracked sols
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
        hc_step++;
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

/*
my_queue, d_startSols_array,  d_Track_array, \
                d_startParams, d_targetParams, d_diffParams, \
                d_Hx_idx, d_Ht_idx, \
                d_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Infinity, \
                d_Debug_Purpose
*/

  real_Double_t
  kernel_HC_Solver_trifocal_2op1p_30x30(
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, 
    magmaFloatComplex** d_Track_array,
    magmaFloatComplex*  d_startParams,
    magmaFloatComplex*  d_targetParams,
    magmaFloatComplex*  d_diffParams,
    magma_int_t*        d_Hx_indx, 
    magma_int_t*        d_Ht_indx,
    bool*               d_is_GPU_HC_Sol_Converge,
    bool*               d_is_GPU_HC_Sol_Infinity,
    magmaFloatComplex*  d_Debug_Purpose
  )
  {
    real_Double_t gpu_time;
    dim3 threads(NUM_OF_VARS, 1, 1);
    dim3 grid(NUM_OF_TRACKS, 1, 1);
    cudaError_t e = cudaErrorInvalidValue;

    //> declare the amount of shared memory for the use of the kernel
    magma_int_t shmem  = 0;
    shmem += (NUM_OF_PARAMS+1) * sizeof(magmaFloatComplex);           //> start parameters
    shmem += (NUM_OF_PARAMS+1) * sizeof(magmaFloatComplex);           //> target parameters
    shmem += (NUM_OF_PARAMS+1) * sizeof(magmaFloatComplex);           //> difference of start and target parameters
    shmem += (NUM_OF_PARAMS+1) * sizeof(magmaFloatComplex);           //> parameter homotopy used when t is changed
    shmem += (NUM_OF_VARS+1)   * sizeof(magmaFloatComplex);           //> start solutions
    shmem += (NUM_OF_VARS+1)   * sizeof(magmaFloatComplex);           //> intermediate solutions
    shmem += (NUM_OF_VARS+1)   * sizeof(magmaFloatComplex);           //> last successful intermediate solutions
    shmem += (NUM_OF_VARS)     * sizeof(magmaFloatComplex);           //> linear system solution
    shmem += (NUM_OF_VARS)     * sizeof(magmaFloatComplex);           //> intermediate varaible for cgesv
    shmem += (NUM_OF_VARS)     * sizeof(float);                       //> intermediate varaible for cgesv
    shmem += (NUM_OF_VARS)     * sizeof(int);                         //> squared solution
    shmem += (NUM_OF_VARS)     * sizeof(float);                       //> squared correction solution
    shmem += (NUM_OF_VARS)     * sizeof(float);                       //> solution norm
    shmem += (2)               * sizeof(float);                       //> pivot for cgesv
    shmem += (1)               * sizeof(bool);                        //> is_successful 
    shmem += (1)               * sizeof(int);                         //> predictor successes counter

    //> declare kernel arguments  
    void *kernel_args[] = {&d_startSols_array, &d_Track_array,
                           &d_startParams, &d_targetParams, &d_diffParams,
                           &d_Hx_indx, &d_Ht_indx,
                           &d_is_GPU_HC_Sol_Converge,
                           &d_is_GPU_HC_Sol_Infinity,
                           &d_Debug_Purpose
                          };

    gpu_time = magma_sync_wtime( my_queue );

    //> launch the GPU kernel
    e = cudaLaunchKernel((void*)HC_solver_trifocal_2op1p_30x30, \
                          grid, threads, kernel_args, shmem, my_queue->cuda_stream());

    gpu_time = magma_sync_wtime( my_queue ) - gpu_time;
    if( e != cudaSuccess ) printf("cudaLaunchKernel of HC_solver_trifocal_2op1p_30x30 is not successful!\n");

    return gpu_time;
  }

//}

#endif
