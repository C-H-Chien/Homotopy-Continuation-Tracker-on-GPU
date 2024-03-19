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
  unsigned ACTIVE_MASK = __activemask();

  //> define pointers to the arrays
  magmaFloatComplex* d_startSols    = d_startSols_array[batchid];
  magmaFloatComplex* d_track        = d_Track_array[batchid];

  //> declarations of registers
  magmaFloatComplex r_cgesvA[NUM_OF_VARS] = {MAGMA_C_ZERO};
  magmaFloatComplex r_cgesvB              = MAGMA_C_ZERO;
  int linfo = 0, rowid = tx;
  float t0 = 0.0, t_step = 0.0, delta_t = 0.01;
  bool end_zone = 0;
#if APPLY_GAMMA_TRICK
  magmaFloatComplex gammified_t0;
  magmaFloatComplex gc;
#endif

  if (tx == 0 && batchid == 0) {
    printf("Inside the kernel, test ...\n");
  }

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
  float* dsx                              = (float*)(sx + NUM_OF_VARS);
  int* sipiv                              = (int*)(dsx + NUM_OF_VARS);
#if USE_LOOPY_RUNGE_KUTTA
  float* s_delta_t_scale                  = (float*)(sipiv + (NUM_OF_VARS+1));
  int* s_RK_Coeffs                        = (int*)(s_delta_t_scale + 1);
#endif

  //> read data from global memory to shared memories or do initializations
  s_sols[tx]               = d_startSols[tx];
  s_track[tx]              = d_track[tx];
  s_track_last_success[tx] = s_track[tx];
  
  //> start and target parameters
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>> REUSE MEMORY? >>>>>>>>>>>>>>>>>>>>>>>
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
    sipiv[NUM_OF_VARS]                = 0;
  }
  __syncthreads();

  //> 1/2 \Delta t
  float one_half_delta_t;
  float r_sqrt_sols;
  float r_sqrt_corr;
  bool r_isSuccessful;
  bool r_isInfFail;

#if USE_LOOPY_RUNGE_KUTTA
  bool scales[3];
  scales[0] = 1;
  scales[1] = 0;
  scales[2] = 1;
#endif

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
#if USE_LOOPY_RUNGE_KUTTA
      if (tx == 0) {
        s_delta_t_scale[0] = 0.0;
        s_RK_Coeffs[0] = 1;
      }
      __syncthreads();

      //> For simplicity, let's stay with no gamma-trick mode
      for (int rk_step = 0; rk_step < 4; rk_step++ ) {

        //> Evaluate parameter homotopy
        compute_param_homotopy< float >( tx, t0, s_param_homotopy, s_startParams, s_targetParams );

        //> Evaluate dH/dx and dH/dt
        eval_Jacobian_Hx< HX_MAXIMAL_TERMS*HX_MAXIMAL_PARTS >( tx, r_cgesvA, s_track, s_startParams, s_targetParams, s_param_homotopy, d_Hx_indices );
        eval_Jacobian_Ht( tx, r_cgesvB, s_track, s_startParams, s_targetParams, s_param_homotopy, d_Ht_indices, s_diffParams );

        //> linear system solver: solve for k1, k2, k3, or k4
        cgesv_batched_small_device< NUM_OF_VARS >( tx, r_cgesvA, sipiv, r_cgesvB, sB, sx, dsx, rowid, linfo );
        magmablas_syncwarp();

        if (rk_step < 3) {

          s_sols[tx] += sB[tx] * delta_t * (s_RK_Coeffs[0] * 1.0/6.0);
          s_track[tx] = (s_RK_Coeffs[0] > 1) ? s_track_last_success[tx] : s_track[tx];
          
          if (tx == 0) {
            s_delta_t_scale[0] += scales[rk_step] * one_half_delta_t;
            s_RK_Coeffs[0] = s_RK_Coeffs[0] << scales[rk_step];           //> Shift one bit
          }
          __syncthreads();

          sB[tx] *= s_delta_t_scale[0];
          s_track[tx] += sB[tx];
          t0 += scales[rk_step] * one_half_delta_t;
        }
        magmablas_syncwarp();
      }
      //> Make prediction
      s_sols[tx] += sB[tx] * delta_t * 1.0/6.0;
      s_track[tx] = s_sols[tx];
      __syncthreads();
#else
      //> get HxHt for k1
#if APPLY_GAMMA_TRICK
      gammified_t0 = GAMMA * t0 / (MAGMA_C_ONE + GAMMA_MINUS_ONE * t0);                                      //> t0
      compute_param_homotopy< magmaFloatComplex >( tx, gammified_t0, s_param_homotopy, s_startParams, s_targetParams );
#else
      compute_param_homotopy< float >( tx, t0, s_param_homotopy, s_startParams, s_targetParams );
#endif
      eval_Jacobian_Hx< HX_MAXIMAL_TERMS*HX_MAXIMAL_PARTS >( tx, r_cgesvA, s_track, s_startParams, s_targetParams, s_param_homotopy, d_Hx_indices );
      eval_Jacobian_Ht( tx, r_cgesvB, s_track, s_startParams, s_targetParams, s_param_homotopy, d_Ht_indices, s_diffParams );
      
      //> solve k1
      cgesv_batched_small_device< NUM_OF_VARS >( tx, r_cgesvA, sipiv, r_cgesvB, sB, sx, dsx, rowid, linfo );
      magmablas_syncwarp();

      //> compute x for the creation of HxHt for k2 and get HxHt for k2
#if APPLY_GAMMA_TRICK
      gc = GAMMA / ((GAMMA_MINUS_ONE * t0 + 1.0) * (GAMMA_MINUS_ONE * t0 + 1.0));
      create_x_for_k2( tx, t0, delta_t, one_half_delta_t, s_sols, s_track, sB, gc );
      gammified_t0 = GAMMA * t0 / (MAGMA_C_ONE + GAMMA_MINUS_ONE * t0); //> After create_x_for_k2, this t0 is actually (t0 + one_half_delta_t)
      magmablas_syncwarp();
      compute_param_homotopy< magmaFloatComplex >( tx, gammified_t0, s_param_homotopy, s_startParams, s_targetParams );
#else
      create_x_for_k2( tx, t0, delta_t, one_half_delta_t, s_sols, s_track, sB, MAGMA_C_ONE );
      magmablas_syncwarp();
      compute_param_homotopy< float >( tx, t0, s_param_homotopy, s_startParams, s_targetParams );
#endif
      eval_Jacobian_Hx< HX_MAXIMAL_TERMS*HX_MAXIMAL_PARTS >( tx, r_cgesvA, s_track, s_startParams, s_targetParams, s_param_homotopy, d_Hx_indices );
      eval_Jacobian_Ht( tx, r_cgesvB, s_track, s_startParams, s_targetParams, s_param_homotopy, d_Ht_indices, s_diffParams );

      //> solve k2
      cgesv_batched_small_device< NUM_OF_VARS >( tx, r_cgesvA, sipiv, r_cgesvB, sB, sx, dsx, rowid, linfo );
      magmablas_syncwarp();

      //> compute x for the generation of HxHt for k3 and get HxHt for k3
#if APPLY_GAMMA_TRICK
      //gc = GAMMA / (((GAMMA - 1.0) * (t0 + one_half_delta_t) + 1.0) * ((GAMMA - 1.0) * (t0 + one_half_delta_t) + 1.0));
      gc = GAMMA / ((GAMMA_MINUS_ONE * t0 + 1.0) * (GAMMA_MINUS_ONE * t0 + 1.0));
      create_x_for_k3( tx, delta_t, one_half_delta_t, s_sols, s_track, s_track_last_success, sB, gc );
      magmablas_syncwarp();
#else
      create_x_for_k3( tx, delta_t, one_half_delta_t, s_sols, s_track, s_track_last_success, sB, MAGMA_C_ONE );
      magmablas_syncwarp();
#endif
      eval_Jacobian_Hx< HX_MAXIMAL_TERMS*HX_MAXIMAL_PARTS >( tx, r_cgesvA, s_track, s_startParams, s_targetParams, s_param_homotopy, d_Hx_indices );
      eval_Jacobian_Ht( tx, r_cgesvB, s_track, s_startParams, s_targetParams, s_param_homotopy, d_Ht_indices, s_diffParams );

      //> solve k3
      cgesv_batched_small_device< NUM_OF_VARS >( tx, r_cgesvA, sipiv, r_cgesvB, sB, sx, dsx, rowid, linfo );
      magmablas_syncwarp();

      //> compute x for the generation of HxHt for k4 and get HxHt for k4
#if APPLY_GAMMA_TRICK
      create_x_for_k4( tx, t0, delta_t, one_half_delta_t, s_sols, s_track, s_track_last_success, sB, gc );
      gammified_t0 = GAMMA * t0 / (MAGMA_C_ONE + GAMMA_MINUS_ONE * t0); //> After create_x_for_k4, this t0 is actually (t0 + delta_t)
      magmablas_syncwarp();
      compute_param_homotopy< magmaFloatComplex >( tx, gammified_t0, s_param_homotopy, s_startParams, s_targetParams );
#else
      create_x_for_k4( tx, t0, delta_t, one_half_delta_t, s_sols, s_track, s_track_last_success, sB, MAGMA_C_ONE );
      magmablas_syncwarp();
      compute_param_homotopy< float >( tx, t0, s_param_homotopy, s_startParams, s_targetParams );
#endif
      eval_Jacobian_Hx< HX_MAXIMAL_TERMS*HX_MAXIMAL_PARTS >( tx, r_cgesvA, s_track, s_startParams, s_targetParams, s_param_homotopy, d_Hx_indices );
      eval_Jacobian_Ht( tx, r_cgesvB, s_track, s_startParams, s_targetParams, s_param_homotopy, d_Ht_indices, s_diffParams );

      //> solve k4
      cgesv_batched_small_device< NUM_OF_VARS >( tx, r_cgesvA, sipiv, r_cgesvB, sB, sx, dsx, rowid, linfo );
      magmablas_syncwarp();

      //> make prediction
#if APPLY_GAMMA_TRICK
      gc = GAMMA / ((GAMMA_MINUS_ONE * t0 + 1.0) * (GAMMA_MINUS_ONE * t0 + 1.0));
      s_sols[tx] += sB[tx] * delta_t * gc * 1.0/6.0;
      s_track[tx] = s_sols[tx];
      __syncthreads();
#else
      s_sols[tx] += sB[tx] * delta_t * 1.0/6.0;
      s_track[tx] = s_sols[tx];
      __syncthreads();
#endif

#endif  //> USE_LOOPY_RUNGE_KUTTA

      // ===================================================================
      //> Gauss-Newton Corrector
      // ===================================================================
      //#pragma unroll
      for(int i = 0; i < HC_MAX_CORRECTION_STEPS; i++) {

        //> evaluate the Jacobian Hx and the parameter homotopy
        eval_Jacobian_Hx<HX_MAXIMAL_TERMS*HX_MAXIMAL_PARTS>( tx, r_cgesvA, s_track, s_startParams, s_targetParams, s_param_homotopy, d_Hx_indices );
        eval_Parameter_Homotopy( tx, r_cgesvB, s_track, s_startParams, s_targetParams, s_param_homotopy, d_Ht_indices );

        //> G-NUM_OF_VARS corrector first solve
        cgesv_batched_small_device< NUM_OF_VARS >( tx, r_cgesvA, sipiv, r_cgesvB, sB, sx, dsx, rowid, linfo );
        magmablas_syncwarp();

        //> correct the sols
        s_track[tx] -= sB[tx];
        __syncthreads();

        r_sqrt_sols = MAGMA_C_REAL(sB[tx])*MAGMA_C_REAL(sB[tx]) + MAGMA_C_IMAG(sB[tx])*MAGMA_C_IMAG(sB[tx]);
        r_sqrt_corr = MAGMA_C_REAL(s_track[tx])*MAGMA_C_REAL(s_track[tx]) + MAGMA_C_IMAG(s_track[tx])*MAGMA_C_IMAG(s_track[tx]);
        __syncthreads();

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
      // Decide Track Changes
      // ===================================================================
      if (!r_isSuccessful) {
        delta_t *= 0.5;
        //> should be the last successful tracked sols
        s_track[tx] = s_track_last_success[tx];
        s_sols[tx] = s_track_last_success[tx];
        if (tx == 0) sipiv[NUM_OF_VARS] = 0;
        __syncthreads();
        t0 = t_step;
      }
      else {
        if (tx == 0) sipiv[NUM_OF_VARS]++;
        s_track_last_success[tx] = s_track[tx];
        s_sols[tx] = s_track[tx];
        __syncthreads();
        if (sipiv[NUM_OF_VARS] >= HC_NUM_OF_STEPS_TO_INCREASE_DELTA_T) {
          if (tx == 0) sipiv[NUM_OF_VARS] = 0;
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
    d_is_GPU_HC_Sol_Infinity[ batchid ] = (r_isInfFail) ? (1) : (0);
  }
#if GPU_DEBUG
  d_Debug_Purpose[ batchid ] = (t0 >= 1.0 || (1.0-t0 <= 0.0000001)) ? MAGMA_C_MAKE(1.0, 0.0) : MAGMA_C_MAKE(t0, delta_t);
#endif
}

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
  shmem += (NUM_OF_VARS+1)   * sizeof(int);                         //> sipiv
  shmem += (1)               * sizeof(int);                         //> predictor successes counter
  shmem += (1)               * sizeof(int);                         //> Loopy Runge-Kutta coefficients
  shmem += (1)               * sizeof(float);                       //> Loopy Runge-Kutta delta t

  //> Get max. dynamic shared memory on the GPU
  int nthreads_max, shmem_max = 0;
  cudacheck( cudaDeviceGetAttribute(&nthreads_max, cudaDevAttrMaxThreadsPerBlock, 0) );
#if CUDA_VERSION >= 9000
  cudacheck( cudaDeviceGetAttribute (&shmem_max, cudaDevAttrMaxSharedMemoryPerBlockOptin, 0) );
  if (shmem <= shmem_max) {
    cudacheck( cudaFuncSetAttribute(HC_solver_trifocal_2op1p_30x30, \
                                    cudaFuncAttributeMaxDynamicSharedMemorySize, shmem) );
  }
#else
  cudacheck( cudaDeviceGetAttribute (&shmem_max, cudaDevAttrMaxSharedMemoryPerBlock, 0) );
#endif

  //> Message of overuse shared memory
  if ( shmem > shmem_max ) printf("Error: kernel %s requires too many threads or too much shared memory\n", __func__);

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

#endif
