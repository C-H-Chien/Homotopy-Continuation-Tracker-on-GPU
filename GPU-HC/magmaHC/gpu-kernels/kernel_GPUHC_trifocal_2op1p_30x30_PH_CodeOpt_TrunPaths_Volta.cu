#ifndef kernel_GPUHC_trifocal_2op1p_30x30_PH_CODEOPT_TRUNPATHS_VOLTA_cu
#define kernel_GPUHC_trifocal_2op1p_30x30_PH_CODEOPT_TRUNPATHS_VOLTA_cu
// ===========================================================================================
// GPU homotopy continuation solver for the trifocal 2op1p 30x30 problem
//
// Major Modifications
//    Chiang-Heng Chien  22-10-03:   Edited from the first version 
//                                   (kernel_HC_Solver_trifocal_2op1p_30.cu)
//    Chiang-Heng Chien  23-12-28:   Add macros
//    Chiang-Heng Chien  24-06-12:   
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
#include "../gpu-idx-evals/dev-eval-indxing-trifocal_2op1p_30x30_inline_LimUnroll.cuh"
#include "../dev-cgesv-batched-small.cuh"
#include "../dev-get-new-data.cuh"

template< int Num_Of_Vars,    int Num_Of_Params, \
          int dHdx_Max_Terms, int dHdx_Max_Parts, int dHdx_Entry_Offset, \
          int dHdt_Max_Terms, int dHdt_Max_Parts, \
          int dHdx_Index_Matrix_Size, int dHdt_Index_Matrix_Size >
__global__ void
kernel_GPUHC_trifocal_pose_PH_CodeOpt_TrunPaths_Volta(
  const int               HC_max_steps, 
  const int               HC_max_correction_steps, 
  const int               HC_delta_t_incremental_steps,
  magmaFloatComplex**     d_startSols_array,
  magmaFloatComplex**     d_Track_array,
  magmaFloatComplex*      d_startParams,
  magmaFloatComplex*      d_targetParams,
  magmaFloatComplex*      d_diffParams,
  const int* __restrict__ d_dHdx_indices,
  const int* __restrict__ d_dHdt_indices,
  bool*                   d_is_GPU_HC_Sol_Converge,
  bool*                   d_is_GPU_HC_Sol_Infinity,
  magmaFloatComplex*      d_Debug_Purpose
)
{
  extern __shared__ magmaFloatComplex zdata[];
  const int tx = threadIdx.x;
  const int batchid = blockIdx.x;
  const int recurrent_batchid = batchid % 312;
  const int ransac_id = batchid / 312;
  unsigned ACTIVE_MASK = __activemask();

  //> define pointers to the arrays
  magmaFloatComplex* d_startSols    = d_startSols_array[recurrent_batchid];
  magmaFloatComplex* d_track        = d_Track_array[batchid];

  //> declarations of registers
  magmaFloatComplex r_cgesvA[Num_Of_Vars] = {MAGMA_C_ZERO};
  magmaFloatComplex r_cgesvB              = MAGMA_C_ZERO;
  int linfo = 0, rowid = tx;
  float t0 = 0.0, t_step = 0.0, delta_t = 0.01;
  bool end_zone = 0;

  //> Dynamic allocated shared memories
  magmaFloatComplex *s_startParams        = (magmaFloatComplex*)(zdata);
  magmaFloatComplex *s_targetParams       = s_startParams        + (Num_Of_Params + 1);
  magmaFloatComplex *s_diffParams         = s_targetParams       + (Num_Of_Params + 1);
  magmaFloatComplex *s_param_homotopy     = s_diffParams         + (Num_Of_Params + 1);
  magmaFloatComplex *s_sols               = s_param_homotopy     + (Num_Of_Params + 1);
  magmaFloatComplex *s_track              = s_sols               + (Num_Of_Vars+1);
  magmaFloatComplex *s_track_last_success = s_track              + (Num_Of_Vars+1);
  magmaFloatComplex *sB                   = s_track_last_success + (Num_Of_Vars+1);
  magmaFloatComplex *sx                   = sB                   + Num_Of_Vars;
  float* dsx                              = (float*)(sx          + Num_Of_Vars);
  int* sipiv                              = (int*)(dsx           + Num_Of_Vars);
  float* s_delta_t_scale                  = (float*)(sipiv       + (Num_Of_Vars+1));
  int* s_RK_Coeffs                        = (int*)(s_delta_t_scale + 1);

  const int* __restrict__ dHdx_indices = d_dHdx_indices;
  const int* __restrict__ dHdt_indices = d_dHdt_indices;

  //> read data from global memory to shared memories or do initializations
  s_sols[tx]               = d_startSols[tx];
  s_track[tx]              = d_track[tx];
  s_track_last_success[tx] = s_track[tx];
  
  //> start and target parameters
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>> REUSE MEMORY >>>>>>>>>>>>>>>>>>>>>>>
  s_startParams[tx]  = d_startParams[tx];
  s_targetParams[tx] = d_targetParams[tx + ransac_id*(Num_Of_Params+1)];
  s_diffParams[tx]   = d_diffParams[tx + ransac_id*(Num_Of_Params+1)];

  if (tx == 0) {
    //> the rest of the start and target parameters
    #pragma unroll
    for(int i = Num_Of_Vars; i <= Num_Of_Params; i++) {
      s_startParams[i]  = d_startParams[i];
      s_targetParams[i] = d_targetParams[i];
      s_diffParams[i]   = d_diffParams[i];
    }
    s_sols[Num_Of_Vars]               = MAGMA_C_MAKE(1.0, 0.0);
    s_track[Num_Of_Vars]              = MAGMA_C_MAKE(1.0, 0.0);
    s_track_last_success[Num_Of_Vars] = MAGMA_C_MAKE(1.0, 0.0);
    s_param_homotopy[Num_Of_Params]   = MAGMA_C_ONE;
    sipiv[Num_Of_Vars]                = 0;
  }
  magmablas_syncwarp();

  //> 1/2 \Delta t
  float one_half_delta_t;
  float r_sqrt_sols;
  float r_sqrt_corr;
  bool r_isSuccessful;
  bool r_isInfFail;

  bool are_Depths_All_Positive = false;
  bool check_depths_sign = true;

  volatile int hc_max_steps = HC_max_steps;
  for (int step = 0; step <= hc_max_steps; step++) {
    if (t0 < 1.0 && (1.0-t0 > 0.0000001)) {

      // ===================================================================
      // Decide delta t at end zone
      // ===================================================================
      if (!end_zone && fabs(1 - t0) <= (0.0500001)) {
        end_zone = true;
      }

      //> Use positive depths to early stop the HC paths      
      if (check_depths_sign) {
        are_Depths_All_Positive = (MAGMA_C_REAL(s_track[0]) > 0) && (MAGMA_C_REAL(s_track[1]) > 0) && (MAGMA_C_REAL(s_track[2]) > 0) && (MAGMA_C_REAL(s_track[3]) > 0) &&
                                  (MAGMA_C_REAL(s_track[4]) > 0) && (MAGMA_C_REAL(s_track[5]) > 0) && (MAGMA_C_REAL(s_track[6]) > 0) && (MAGMA_C_REAL(s_track[7]) > 0);
        if (t0 > 0) check_depths_sign = are_Depths_All_Positive ? false : true;
      }
      if (t0 > 0.95 && check_depths_sign) break;

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
      unsigned char scales[3] = {1, 0, 1};
      if (tx == 0) {
        s_delta_t_scale[0] = 0.0;
        s_RK_Coeffs[0] = 1;
      }
      magmablas_syncwarp();

      //> For simplicity, let's stay with no gamma-trick mode
      for (int rk_step = 0; rk_step < 4; rk_step++ ) {

        //> Evaluate parameter homotopy
        compute_param_homotopy< float, Num_Of_Vars >( tx, t0, s_param_homotopy, s_startParams, s_targetParams );

        //> Evaluate dH/dx and dH/dt
        eval_Jacobian_Hx< Num_Of_Vars, dHdx_Max_Terms, dHdx_Max_Parts, dHdx_Entry_Offset, dHdx_Index_Matrix_Size >( tx, r_cgesvA, s_track, s_startParams, s_targetParams, s_param_homotopy, dHdx_indices );
        eval_Jacobian_Ht< Num_Of_Vars, dHdt_Max_Terms, dHdt_Max_Parts, dHdt_Index_Matrix_Size >( tx, r_cgesvB, s_track, s_startParams, s_targetParams, s_param_homotopy, dHdt_indices, s_diffParams );

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
      volatile int hc_max_correction_steps = HC_max_correction_steps;
      for(int i = 0; i < hc_max_correction_steps; i++) {

        //> evaluate the Jacobian Hx and the parameter homotopy
        eval_Jacobian_Hx< Num_Of_Vars, dHdx_Max_Terms, dHdx_Max_Parts, dHdx_Entry_Offset, dHdx_Index_Matrix_Size >( tx, r_cgesvA, s_track, s_startParams, s_targetParams, s_param_homotopy, dHdx_indices );
        eval_Homotopy< Num_Of_Vars, dHdt_Max_Terms, dHdt_Max_Parts, dHdt_Index_Matrix_Size >( tx, r_cgesvB, s_track, s_startParams, s_targetParams, s_param_homotopy, dHdt_indices );

        //> G-Num_Of_Vars corrector first solve
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
      // Decide Track Changes
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
kernel_GPUHC_trifocal_2op1p_30x30_PH_CodeOpt_TrunPaths_Volta(
  magma_queue_t       my_queue,
  int                 sub_RANSAC_iters,
  int                 HC_max_steps, 
  int                 HC_max_correction_steps, 
  int                 HC_delta_t_incremental_steps,
  magmaFloatComplex** d_startSols_array, 
  magmaFloatComplex** d_Track_array,
  magmaFloatComplex*  d_startParams,
  magmaFloatComplex*  d_targetParams,
  magmaFloatComplex*  d_diffParams,
  int*                d_dHdx_indx, 
  int*                d_dHdt_indx,
  bool*               d_is_GPU_HC_Sol_Converge,
  bool*               d_is_GPU_HC_Sol_Infinity,
  magmaFloatComplex*  d_Debug_Purpose
)
{
  //> Hard-coded for each problem
  const int num_of_params   = 33;
  const int num_of_vars     = 30;
  const int num_of_tracks   = 312;
  const int dHdx_Max_Terms  = 8;
  const int dHdx_Max_Parts  = 5;
  const int dHdt_Max_Terms  = 16;
  const int dHdt_Max_Parts  = 6;
  const int dHdx_Entry_Offset = dHdx_Max_Terms * dHdx_Max_Parts;
  const int dHdx_Index_Matrix_Size = num_of_vars * num_of_vars * dHdx_Max_Terms * dHdx_Max_Parts;
  const int dHdt_Index_Matrix_Size = num_of_vars * dHdt_Max_Terms * dHdt_Max_Parts;

  real_Double_t gpu_time = 0.0;
  dim3 threads(num_of_vars, 1, 1);
  dim3 grid(num_of_tracks*sub_RANSAC_iters, 1, 1);
  cudaError_t e = cudaErrorInvalidValue;

  //> declare the amount of shared memory for the use of the kernel
  magma_int_t shmem  = 0;
  shmem += (num_of_params+1) * sizeof(magmaFloatComplex);           //> start parameters
  shmem += (num_of_params+1) * sizeof(magmaFloatComplex);           //> target parameters
  shmem += (num_of_params+1) * sizeof(magmaFloatComplex);           //> difference of start and target parameters
  shmem += (num_of_params+1) * sizeof(magmaFloatComplex);           //> parameter homotopy used when t is changed
  shmem += (num_of_vars+1)   * sizeof(magmaFloatComplex);           //> start solutions
  shmem += (num_of_vars+1)   * sizeof(magmaFloatComplex);           //> intermediate solutions
  shmem += (num_of_vars+1)   * sizeof(magmaFloatComplex);           //> last successful intermediate solutions
  shmem += (num_of_vars)     * sizeof(magmaFloatComplex);           //> linear system solution
  shmem += (num_of_vars)     * sizeof(magmaFloatComplex);           //> intermediate varaible for cgesv
  shmem += (num_of_vars)     * sizeof(float);                       //> intermediate varaible for cgesv
  shmem += (num_of_vars+1)   * sizeof(int);                         //> sipiv
  shmem += (1)               * sizeof(int);                         //> predictor successes counter
  shmem += (1)               * sizeof(int);                         //> Loopy Runge-Kutta coefficients
  shmem += (1)               * sizeof(float);                       //> Loopy Runge-Kutta delta t

  //> Get max. dynamic shared memory on the GPU
  int nthreads_max, shmem_max = 0;
  cudacheck( cudaDeviceGetAttribute(&nthreads_max, cudaDevAttrMaxThreadsPerBlock, 0) );
#if CUDA_VERSION >= 9000
  cudacheck( cudaDeviceGetAttribute (&shmem_max, cudaDevAttrMaxSharedMemoryPerBlockOptin, 0) );
  if (shmem <= shmem_max) {
    cudacheck( cudaFuncSetAttribute(kernel_GPUHC_trifocal_pose_PH_CodeOpt_TrunPaths_Volta \
                                    <num_of_vars, num_of_params, \
                                     dHdx_Max_Terms, dHdx_Max_Parts, dHdx_Entry_Offset, dHdt_Max_Terms, dHdt_Max_Parts, \
                                     dHdx_Index_Matrix_Size, dHdt_Index_Matrix_Size >, //dHdx_Num_Of_Read_Loops, dHdt_Num_Of_Read_Loops >,
                                    cudaFuncAttributeMaxDynamicSharedMemorySize, shmem) );
  }
#else
  cudacheck( cudaDeviceGetAttribute (&shmem_max, cudaDevAttrMaxSharedMemoryPerBlock, 0) );
#endif

  //> Message of overuse shared memory
  if ( shmem > shmem_max ) printf("Error: kernel %s requires too many threads or too much shared memory\n", __func__);

  //> declare kernel arguments  
  void *kernel_args[] = { &HC_max_steps, 
                          &HC_max_correction_steps, 
                          &HC_delta_t_incremental_steps,
                          &d_startSols_array, &d_Track_array,
                          &d_startParams, &d_targetParams, &d_diffParams,
                          &d_dHdx_indx, &d_dHdt_indx,
                          &d_is_GPU_HC_Sol_Converge,
                          &d_is_GPU_HC_Sol_Infinity,
                          &d_Debug_Purpose
                        };

  // gpu_time = magma_sync_wtime( my_queue );

  //> launch the GPU kernel
  e = cudaLaunchKernel((void*)kernel_GPUHC_trifocal_pose_PH_CodeOpt_TrunPaths_Volta \
                        <num_of_vars, num_of_params, \
                         dHdx_Max_Terms, dHdx_Max_Parts, dHdx_Entry_Offset, dHdt_Max_Terms, dHdt_Max_Parts, \
                         dHdx_Index_Matrix_Size, dHdt_Index_Matrix_Size >, // dHdx_Num_Of_Read_Loops, dHdt_Num_Of_Read_Loops >,
                        grid, threads, kernel_args, shmem, my_queue->cuda_stream());

  // gpu_time = magma_sync_wtime( my_queue ) - gpu_time;
  if( e != cudaSuccess ) printf("cudaLaunchKernel of kernel_GPUHC_trifocal_pose_PH_CodeOpt_TrunPaths_Volta is not successful!\n");

  return gpu_time;
}

#endif
