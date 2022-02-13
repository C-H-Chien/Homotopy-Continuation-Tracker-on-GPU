#ifndef kernel_HC_Solver_5pt_rel_pose_w_depth_recon_cu
#define kernel_HC_Solver_5pt_rel_pose_w_depth_recon_cu
// ============================================================================
// GPU homotopy continuation solver for 5pt_rel_pose_w_depth_recon problem
//
// Modifications
//    Chien  21-12-31:   originally created
//
// ============================================================================
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>

// cuda included
#include <cuda.h>
#include <cuda_runtime.h>

// magma
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

// -- header --
#include "magmaHC-kernels.h"

// -- device function --
#include "../gpu-idx-evals/dev-eval-indxing-5pt_rel_pose_w_depth_recon.cuh"
#include "../dev-cgesv-batched-small.cuh"
#include "../dev-get-new-data.cuh"

namespace magmaHCWrapper {

  template<int N, int phc_coeffsCount, int max_steps, int max_corr_steps, int predSuccessCount, int Hx_max_terms, int Hx_max_parts, int Ht_max_terms, int Ht_max_parts>
  __global__ void
  homotopy_continuation_solver_5pt_rel_pose_w_depth_recon(
    magma_int_t ldda,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array,
    magmaFloatComplex_ptr d_phc_coeffs_Hx, magmaFloatComplex_ptr d_phc_coeffs_Ht
  )
  {
    extern __shared__ magmaFloatComplex zdata[];
    const int tx = threadIdx.x;
    const int batchid = blockIdx.x ;

    magmaFloatComplex* d_startSols = d_startSols_array[batchid];
    magmaFloatComplex* d_track = d_Track_array[batchid];
    magmaFloatComplex* d_cgesvA = d_cgesvA_array[batchid];
    magmaFloatComplex* d_cgesvB = d_cgesvB_array[batchid];
    const magmaFloatComplex* __restrict__ d_const_phc_coeffs_Hx = d_phc_coeffs_Hx;
    const magmaFloatComplex* __restrict__ d_const_phc_coeffs_Ht = d_phc_coeffs_Ht;
    const int* __restrict__ d_Hx_idx = d_Hx_idx_array[0];
    const int* __restrict__ d_Ht_idx = d_Ht_idx_array[0];

    // -- registers declarations --
    magmaFloatComplex r_cgesvA[N] = {MAGMA_C_ZERO};
    magmaFloatComplex r_cgesvB = MAGMA_C_ZERO;
    int linfo = 0, rowid = tx;
    float t0 = 0.0, t_step = 0.0, delta_t = 0.01;
    bool inf_failed = 0;
    bool end_zone = 0;
    int hc_step = 0;

    // -- shared memory declarations --
    magmaFloatComplex *s_sols = (magmaFloatComplex*)(zdata);
    magmaFloatComplex *s_track = s_sols + (N+1);
    magmaFloatComplex *s_track_last_success = s_track + (N+1);
    magmaFloatComplex *sB = s_track_last_success + (N+1);
    magmaFloatComplex *sx = sB + N;
    magmaFloatComplex *s_phc_coeffs_Hx = sx + N;
    magmaFloatComplex *s_phc_coeffs_Ht = s_phc_coeffs_Hx + phc_coeffsCount;
    float* dsx = (float*)(s_phc_coeffs_Ht + phc_coeffsCount);
    int* sipiv = (int*)(dsx + N);
    float *s_sqrt_sols = (float*)(sipiv + N);
    float *s_sqrt_corr = s_sqrt_sols + N;
    float *s_norm = s_sqrt_corr + N;
    bool s_isSuccessful = (bool)(s_norm + 2);
    int s_pred_success_count = (int)(s_isSuccessful + 1);

    // -- initialization: read from gm --
    #pragma unroll
    for(int i = 0; i < N; i++) {
      r_cgesvA[i] = d_cgesvA[ i * ldda + tx ];
    }
    r_cgesvB = d_cgesvB[tx];

    s_sols[tx] = d_startSols[tx];
    s_track[tx] = d_track[tx];
    s_track_last_success[tx] = s_track[tx];
    s_sqrt_sols[tx] = 0;
    s_sqrt_corr[tx] = 0;
    s_isSuccessful = 0;
    s_pred_success_count = 0;
    if (tx == 0) {
      s_sols[N] = MAGMA_C_MAKE(1.0, 0.0);
      s_track[N] = MAGMA_C_MAKE(1.0, 0.0);
      s_track_last_success[N] = MAGMA_C_MAKE(1.0, 0.0);
    }
    __syncthreads();

    float one_half_delta_t;   // -- 1/2 \Delta t --

    //#pragma unroll
    for (int step = 0; step <= max_steps; step++) {
      if (t0 < 1.0 && (1.0-t0 > 0.00001)) {
        t_step = t0;
        one_half_delta_t = 0.5 * delta_t;
        // ===================================================================
        // -- Runge-Kutta Predictor --
        // ===================================================================
        // -- get HxHt for k1 --
        eval_parameter_homotopy<N, N*2>( tx, t0, s_phc_coeffs_Hx, s_phc_coeffs_Ht, d_const_phc_coeffs_Hx, d_const_phc_coeffs_Ht );
        eval_Jacobian_5pt_rel_pose_w_depth_recon< N, Hx_max_terms, Hx_max_parts, Hx_max_terms*Hx_max_parts, N*Hx_max_terms*Hx_max_parts>( tx, s_track, r_cgesvA, d_Hx_idx, s_phc_coeffs_Hx );
        eval_Ht_5pt_rel_pose_w_depth_recon< N, Ht_max_terms, Ht_max_parts, Ht_max_terms*Ht_max_parts >( tx, s_track, r_cgesvB, d_Ht_idx, s_phc_coeffs_Ht );

        // -- solve k1 --
        cgesv_batched_small_device<N>( tx, r_cgesvA, sipiv, r_cgesvB, sB, sx, dsx, rowid, linfo );
        magmablas_syncwarp();

        // -- compute x for the creation of HxHt for k2 --
        create_x_for_k2( tx, t0, delta_t, one_half_delta_t, s_sols, s_track, sB );
        magmablas_syncwarp();

        // -- get HxHt for k2 --
        eval_parameter_homotopy<N, N*2>( tx, t0, s_phc_coeffs_Hx, s_phc_coeffs_Ht, d_const_phc_coeffs_Hx, d_const_phc_coeffs_Ht );
        eval_Jacobian_5pt_rel_pose_w_depth_recon< N, Hx_max_terms, Hx_max_parts, Hx_max_terms*Hx_max_parts, N*Hx_max_terms*Hx_max_parts>( tx, s_track, r_cgesvA, d_Hx_idx, s_phc_coeffs_Hx );
        eval_Ht_5pt_rel_pose_w_depth_recon< N, Ht_max_terms, Ht_max_parts, Ht_max_terms*Ht_max_parts >( tx, s_track, r_cgesvB, d_Ht_idx, s_phc_coeffs_Ht );

        // -- solve k2 --
        cgesv_batched_small_device<N>( tx, r_cgesvA, sipiv, r_cgesvB, sB, sx, dsx, rowid, linfo );
        magmablas_syncwarp();

        // -- compute x for the generation of HxHt for k3 --
        create_x_for_k3( tx, delta_t, one_half_delta_t, s_sols, s_track, s_track_last_success, sB );
        magmablas_syncwarp();

        // -- get HxHt for k3 --
        eval_parameter_homotopy<N, N*2>( tx, t0, s_phc_coeffs_Hx, s_phc_coeffs_Ht, d_const_phc_coeffs_Hx, d_const_phc_coeffs_Ht );
        eval_Jacobian_5pt_rel_pose_w_depth_recon< N, Hx_max_terms, Hx_max_parts, Hx_max_terms*Hx_max_parts, N*Hx_max_terms*Hx_max_parts>( tx, s_track, r_cgesvA, d_Hx_idx, s_phc_coeffs_Hx );
        eval_Ht_5pt_rel_pose_w_depth_recon< N, Ht_max_terms, Ht_max_parts, Ht_max_terms*Ht_max_parts >( tx, s_track, r_cgesvB, d_Ht_idx, s_phc_coeffs_Ht );

        // -- solve k3 --
        cgesv_batched_small_device<N>( tx, r_cgesvA, sipiv, r_cgesvB, sB, sx, dsx, rowid, linfo );
        magmablas_syncwarp();

        // -- compute x for the generation of HxHt for k4 --
        create_x_for_k4( tx, t0, delta_t, one_half_delta_t, s_sols, s_track, s_track_last_success, sB );
        magmablas_syncwarp();

        // -- get HxHt for k4 --
        eval_parameter_homotopy<N, N*2>( tx, t0, s_phc_coeffs_Hx, s_phc_coeffs_Ht, d_const_phc_coeffs_Hx, d_const_phc_coeffs_Ht );
        eval_Jacobian_5pt_rel_pose_w_depth_recon< N, Hx_max_terms, Hx_max_parts, Hx_max_terms*Hx_max_parts, N*Hx_max_terms*Hx_max_parts>( tx, s_track, r_cgesvA, d_Hx_idx, s_phc_coeffs_Hx );
        eval_Ht_5pt_rel_pose_w_depth_recon< N, Ht_max_terms, Ht_max_parts, Ht_max_terms*Ht_max_parts >( tx, s_track, r_cgesvB, d_Ht_idx, s_phc_coeffs_Ht );

        // -- solve k4 --
        cgesv_batched_small_device<N>( tx, r_cgesvA, sipiv, r_cgesvB, sB, sx, dsx, rowid, linfo );
        magmablas_syncwarp();

        // -- make prediction --
        s_sols[tx] += sB[tx] * delta_t * 1.0/6.0;
        s_track[tx] = s_sols[tx];
        __syncthreads();

        // ===================================================================
        // -- Gauss-Newton Corrector --
        // ===================================================================
        //#pragma unroll
        for(int i = 0; i < max_corr_steps; i++) {

          eval_parameter_homotopy<N, N*2>( tx, t0, s_phc_coeffs_Hx, s_phc_coeffs_Ht, d_const_phc_coeffs_Hx, d_const_phc_coeffs_Ht );
          eval_Jacobian_5pt_rel_pose_w_depth_recon< N, Hx_max_terms, Hx_max_parts, Hx_max_terms*Hx_max_parts, N*Hx_max_terms*Hx_max_parts>( tx, s_track, r_cgesvA, d_Hx_idx, s_phc_coeffs_Hx );
          eval_H_5pt_rel_pose_w_depth_recon< N, Ht_max_terms, Ht_max_parts, Ht_max_terms*Ht_max_parts >( tx, s_track, r_cgesvB, d_Ht_idx, s_phc_coeffs_Hx );

          // -- G-N corrector first solve --
          cgesv_batched_small_device<N>( tx, r_cgesvA, sipiv, r_cgesvB, sB, sx, dsx, rowid, linfo );
          magmablas_syncwarp();

          // -- correct the sols --
          s_track[tx] -= sB[tx];
          __syncthreads();

          // -- compute the norms; norm[0] is norm(sB), norm[1] is norm(sol) --
          compute_norm2<N>( tx, sB, s_track, s_sqrt_sols, s_sqrt_corr, s_norm );
          __syncthreads();
          
          s_isSuccessful = s_norm[0] < 0.000001 * s_norm[1];
          __syncthreads();

          if (s_isSuccessful)
	           break;
        }

        // -- stop if the values of the solution is too large --
        if (s_norm[1] > 1e14) {
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
          if (s_pred_success_count >= predSuccessCount) {
            s_pred_success_count = 0;
            delta_t *= 2;
          }
        }

        // ===================================================================
        // -- Decide delta t at end zone --
        // ===================================================================
        end_zone = (!end_zone && fabs(1 - t0) <= 0.0500001);
        if (end_zone)
            if (delta_t > fabs(1 - t0))
                delta_t = fabs(1 - t0);
        else if (delta_t > fabs(1 - 0.05 - t0))
            delta_t = fabs(1 - 0.05 - t0);

        hc_step++;
      }
      else {
        break;
      }
    }

    // -- write back to gm --
    /*#pragma unroll
    for(int i = 0; i < N; i++){
        d_cgesvA[ i * ldda + rowid ] = r_cgesvA[i];
    }*/

    // -- d_cgesvB tells whether the track is finished, if not, stores t0 and delta_t --
    d_cgesvB[tx] = (t0 >= 1.0 || (1.0-t0 < 0.00001)) ? MAGMA_C_MAKE(1.0, hc_step) : MAGMA_C_MAKE(t0, delta_t);
    //d_cgesvB[tx] = r_cgesvB;
    //d_cgesvB[tx] = sB[tx];

    // -- d_startSols tells the pred_success_count and inf_failed for unfinished tracks --
    //d_startSols[tx] = MAGMA_C_MAKE(s_pred_success_count, inf_failed);

    // -- d_track stores the solutions --
    d_track[tx] = s_track[tx];
    //d_track[tx] = s_phc_coeffs_Ht[tx];
  }

  extern "C" real_Double_t
  kernel_HC_Solver_5pt_rel_pose_w_depth_recon(
    magma_int_t N, magma_int_t batchCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array,
    magmaFloatComplex_ptr d_phc_coeffs_Hx, magmaFloatComplex_ptr d_phc_coeffs_Ht,
    magma_int_t numOf_phc_coeffs
  )
  {
    real_Double_t gpu_time;
    const magma_int_t thread_x = N;
    dim3 threads(thread_x, 1, 1);
    dim3 grid(batchCount, 1, 1);
    cudaError_t e = cudaErrorInvalidValue;

    // -- decalre shared memory --
    magma_int_t shmem  = 0;
    shmem += (N+1) * sizeof(magmaFloatComplex);       // startSols
    shmem += (N+1) * sizeof(magmaFloatComplex);       // track
    shmem += (N+1) * sizeof(magmaFloatComplex);       // track_pred_init

    shmem += numOf_phc_coeffs * sizeof(magmaFloatComplex);  // -- s_phc_coeffs_Hx --
    shmem += numOf_phc_coeffs * sizeof(magmaFloatComplex);  // -- s_phc_coeffs_Ht --

    shmem += N * sizeof(magmaFloatComplex); // sB
    shmem += N * sizeof(magmaFloatComplex); // sx
    shmem += N * sizeof(float);            // dsx
    shmem += N * sizeof(int);               // pivot
    shmem += N * sizeof(float);             // s_sqrt for sol norm-2 in G-N corrector
    shmem += N * sizeof(float);             // s_sqrt for corr norm-2 in G-N corrector
    shmem += 2 * sizeof(float);             // s_norm for norm-2 in G-N corrector
    shmem += 1 * sizeof(bool);              // is_successful 
    shmem += 1 * sizeof(int);               // predictor_success counter

    void *kernel_args[] = {&ldda, &d_startSols_array, &d_Track_array, &d_cgesvA_array, &d_cgesvB_array,
                           &d_Hx_idx_array, &d_Ht_idx_array, &d_phc_coeffs_Hx, &d_phc_coeffs_Ht};

    gpu_time = magma_sync_wtime( my_queue );

    // -- <N, numOfCoeffs, max_steps, max_corr_steps, successes_to_incremental_factor, Hx_max_terms, Hx_max_parts, Ht_max_terms, Ht_max_parts> --
    // -- N*Hx_max_terms*Hx_max_parts = 11*3*3=99 --
    // -- Ht_max_terms*Ht_max_parts = 9*4 = 36 --
    e = cudaLaunchKernel((void*)homotopy_continuation_solver_5pt_rel_pose_w_depth_recon< 16, 44, 125, 4, 5, 7, 4, 9, 5 >, grid, threads, kernel_args, shmem, my_queue->cuda_stream());

    gpu_time = magma_sync_wtime( my_queue ) - gpu_time;
    if( e != cudaSuccess ) {
        printf("cudaLaunchKernel of homotopy_continuation_solver_5pt_rel_pose_w_depth_recon is not successful!\n");
    }

    return gpu_time;
  }

}

#endif
