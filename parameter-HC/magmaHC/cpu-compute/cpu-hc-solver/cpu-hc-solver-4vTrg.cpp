#ifndef cpu_hc_solver_4vTrg_cpp
#define cpu_hc_solver_4vTrg_cpp
// =======================================================================
// Solve homotopy continuation on cpu for the 4 view triangulation problem
//
// Modifications
//    Chien  21-12-29:   Initially Created
//
// =======================================================================
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <chrono>
#include <vector>

// -- magma --
#include "flops.h"
#include "magma_v2.h"
#include "magma_lapack.h"

// -- openmp --
#include <omp.h>

// -- cpu compute header files --
#include "../cpu-compute.h"
#include "../cpu-eval-HxHt/cpu-eval-HxHt-4vTrg.h"
#include "../cpu-eval-HxH/cpu-eval-HxH-4vTrg.h"

namespace magmaHCWrapper {

  extern "C"
  real_Double_t cpu_hc_solver_4vTrg(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_Track_gpu,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    bool continue_from_gpu, magma_queue_t my_queue, int max_steps,
    magmaFloatComplex *h_track_numOfPred_s)
  {
    // -- define parameters --
    const int successes_of_increase_factor = 10;
    const int correction_iterations = 5;
    const float end_zone_factor = 0.05;
    real_Double_t     cpu_time;

    const magmaFloatComplex C0 = MAGMA_C_ONE + MAGMA_C_ONE;
    const magmaFloatComplex C1 = MAGMA_C_ZERO;
    const magmaFloatComplex C2 = MAGMA_C_ONE;
    const magmaFloatComplex C3 = MAGMA_C_NEG_ONE;

    int numOfGpuSuccesses = 0;
    int cpu_inf_failed_count = 0;
    int gpu_inf_failed_count = 0;
    int cpu_track_success_count = 0;
    magma_int_t *ipiv;
    magmaFloatComplex *h_Track_last_success;
    //magmaFloatComplex *h_Track_cpu;
    magma_imalloc_cpu( &ipiv, batchCount*N );
    magma_cmalloc_cpu( &h_Track_last_success, batchCount*N );
    //magma_cmalloc_cpu( &h_Track_cpu, batchCount*N );

    // -- declare and allocate coefficients arrays for openmp parallelization --
    magmaFloatComplex *h_startCoefs_p;
    magmaFloatComplex *h_targetCoefs_p;
    magma_cmalloc_cpu( &h_startCoefs_p, batchCount*coefsCount );
    magma_cmalloc_cpu( &h_targetCoefs_p, batchCount*coefsCount );
    for(int sp = 0; sp < batchCount; sp++) {
      for(int cf = 0; cf < coefsCount; cf++) {
        (h_startCoefs_p + sp * coefsCount)[cf] = h_startCoefs[cf];
        (h_targetCoefs_p + sp * coefsCount)[cf] = h_targetCoefs[cf];
      }
    }

    // -- write HC solutions to the file --
    std::ofstream track_sols_file;
    track_sols_file.open("/users/cchien3/data/cchien3/MyBitBucket/issac-parametric-hc/cpu_converged_HC_tracks.txt");
    std::ofstream cpu_hc_steps_file;
    cpu_hc_steps_file.open("/users/cchien3/data/cchien3/MyBitBucket/issac-parametric-hc/cpu_hc_steps.txt");
    if ( !track_sols_file.is_open() )
       std::cout<<"files cannot be opened!"<<std::endl;

    cpu_time = magma_sync_wtime( my_queue );
    #if !defined (BATCHED_DISABLE_PARCPU) && defined(_OPENMP)
      magma_int_t nthreads = magma_get_lapack_numthreads();
      omp_set_num_threads(nthreads);
      openblas_set_num_threads(1);
      std::cout<<"get lapack numthreads: "<<nthreads<<std::endl;
    #pragma omp parallel for schedule(dynamic) reduction(+:numOfGpuSuccesses, cpu_inf_failed_count, gpu_inf_failed_count, cpu_track_success_count)
    #endif
    for (magma_int_t s=0; s < batchCount; s++) {
      // -- local declarations for openmp parallelization to avoid race condition --
      int pred_success_count = 0;
      bool track_success = 0;
      int track_steps = 0;
      float t0, t_step, delta_t;
      bool inf_failed = 0;
      float sqrt_sols = 0.0;
      float sqrt_corr = 0.0;
      bool isSuccessful = 0;
      float one_half_delta_t;   // -- 1/2 \Delta t --
      bool gpu_track_inf_failed = 0;
      bool end_zone = false;

      magma_int_t locinfo;
      magma_int_t nrhs = 1;
      
      // -- define variables --
      pred_success_count = 0;
      track_success = 0;
      t0 = 0.0;
      t_step = 0.0;
      delta_t = 0.01;
      //end_zone = false;
      inf_failed = 0;
      for (int i = 0; i < N; i++) {
        (h_Track_last_success + s * N)[i] = (h_Track_gpu + s * (N+1))[i];
        (h_Track_cpu + s * N)[i] = (h_Track_gpu + s * (N+1))[i];
      }

      // -- proceed the homotopy continuation algorithm --
      for(int step = 0; step <= max_steps; step++) {
        if (t0 < 1.0 && (1.0-t0 > 0.00001)) {
          t_step = t0;
          one_half_delta_t = 0.5 * delta_t;
          // ===================================================================
          // -- Runge-Kutta Predictor --
          // ===================================================================
          // -- get k1 HxHt --
          cpu_eval_HxHt_4vTrg(s, t0, N, h_Track_cpu + s * N, C0, C1, C2, C3, h_startCoefs_p + s * coefsCount, h_targetCoefs_p + s * coefsCount, h_cgesvA + s * N * N, h_cgesvB + s * N );

          /*if (s == 0 && step == 0) {
            magma_cprint(N, 1, h_cgesvB + s * N, N);
            magma_cprint(N, N, h_cgesvA + s * N * N, N);
          }*/

          // -- solve k1 --
          lapackf77_cgesv( &N, &nrhs, h_cgesvA + s * N * N, &N, ipiv + s * N, h_cgesvB + s * N, &N, &locinfo );
          if (locinfo != 0) {
              printf("lapackf77_cgesv matrix %lld returned error %lld: %s.\n", (long long) s, (long long) locinfo, magma_strerror( locinfo ));
          }

          // -- prepare data for k2 --
          for (int i = 0; i < N; i++) {
            (h_Sols_cpu + s * N)[i] += (h_cgesvB + s * N)[i] * delta_t * 1.0/6.0;    // -- s = s + (\Delta t) (k1/6) --
            (h_cgesvB + s * N)[i] *= one_half_delta_t;                                // -- k1 * (\Delta t)/2 --
            (h_Track_cpu + s * N)[i] += (h_cgesvB + s * N)[i];                        // -- x = x + k1 * (\Delta t)/2 --
          }
          t0 += one_half_delta_t;                                                     // -- t = t + (\Delta t)/2 --

          // -- get k2 HxHt --
          cpu_eval_HxHt_4vTrg(s, t0, N, h_Track_cpu + s * N, C0, C1, C2, C3, h_startCoefs_p + s * coefsCount, h_targetCoefs_p + s * coefsCount, h_cgesvA + s * N * N, h_cgesvB + s * N );

          /*if (s == 0 && step == 0) {
            magma_cprint(N, 1, h_cgesvB + s * N, N);
            magma_cprint(N, N, h_cgesvA + s * N * N, N);
          }*/

          // -- solve k2 --
          lapackf77_cgesv( &N, &nrhs, h_cgesvA + s * N * N, &N, ipiv + s * N, h_cgesvB + s * N, &N, &locinfo );
          if (locinfo != 0) {
              printf("lapackf77_cgesv matrix %lld returned error %lld: %s.\n", (long long) s, (long long) locinfo, magma_strerror( locinfo ));
          }

          // -- prepare data for k3 --
          for (int i = 0; i < N; i++) {
            (h_Sols_cpu + s * N)[i] += (h_cgesvB + s * N)[i] * delta_t * 1.0/3.0;     // -- s = s + (\Delta t) (k1/6 + k2/3) --
            (h_Track_cpu + s * N)[i] = (h_Track_last_success + s * N)[i];          // -- copy the initial prior prediction solution --
            (h_cgesvB + s * N)[i] *= one_half_delta_t;                   // -- k2 * (\Delta t)/2 --
            (h_Track_cpu + s * N)[i] += (h_cgesvB + s * N)[i];                        // -- x = x + k2 * (\Delta t)/2 --
          }

          // -- get k3 HxHt --
          cpu_eval_HxHt_4vTrg(s, t0, N, h_Track_cpu + s * N, C0, C1, C2, C3, h_startCoefs_p + s * coefsCount, h_targetCoefs_p + s * coefsCount, h_cgesvA + s * N * N, h_cgesvB + s * N );

          /*if (s == 0 && step == 0) {
            magma_cprint(N, 1, h_cgesvB + s * N, N);
            magma_cprint(N, N, h_cgesvA + s * N * N, N);
          }*/

          // -- solve k3 --
          lapackf77_cgesv( &N, &nrhs, h_cgesvA + s * N * N, &N, ipiv + s * N, h_cgesvB + s * N, &N, &locinfo );
          if (locinfo != 0) {
              printf("lapackf77_cgesv matrix %lld returned error %lld: %s.\n", (long long) s, (long long) locinfo, magma_strerror( locinfo ));
          }

          // -- prepare data for k4 --
          for (int i = 0; i < N; i++) {
            (h_Sols_cpu + s * N)[i] += (h_cgesvB + s * N)[i] * delta_t * 1.0/3.0;     // -- s = s + (\Delta t) (k1/6 + k2/3 + k3/3) --
            (h_Track_cpu + s * N)[i] = (h_Track_last_success + s * N)[i];                            // -- copy the initial prior prediction solution --
            (h_cgesvB + s * N)[i] *= delta_t;                                          // -- k3 * (\Delta t) --
            (h_Track_cpu + s * N)[i] += (h_cgesvB + s * N)[i];                                            // -- x = x + k3 * (\Delta t) --
          }
          t0 += one_half_delta_t;                        // -- now t becomes t = t + (\Delta t) --

          // -- get k4 HxHt --
          cpu_eval_HxHt_4vTrg(s, t0, N, h_Track_cpu + s * N, C0, C1, C2, C3, h_startCoefs_p + s * coefsCount, h_targetCoefs_p + s * coefsCount, h_cgesvA + s * N * N, h_cgesvB + s * N );

          // -- solve k4 --
          lapackf77_cgesv( &N, &nrhs, h_cgesvA + s * N * N, &N, ipiv + s * N, h_cgesvB + s * N, &N, &locinfo );
          if (locinfo != 0) {
              printf("lapackf77_cgesv matrix %lld returned error %lld: %s.\n", (long long) s, (long long) locinfo, magma_strerror( locinfo ));
          }

          // -- make prediction --
          for (int i = 0; i < N; i++) {
            (h_Sols_cpu + s * N)[i] += (h_cgesvB + s * N)[i] * delta_t * 1.0/6.0;
            (h_Track_cpu + s * N)[i] = (h_Sols_cpu + s * N)[i];
          }

          /*if (s == 0 && step == 0) {
            magma_cprint(N, 1, h_Track_cpu + s * N, N);
          }*/

          // ===================================================================
          // -- Gauss-Newton Corrector --
          // ===================================================================
          for(int c = 0; c < correction_iterations; c++) {
            sqrt_sols = 0;
            sqrt_corr = 0;

            // -- create correction term --
            cpu_eval_HxH_4vTrg( s, t0, N, h_Track_cpu + s * N, C0, C1, C2, C3, h_startCoefs_p + s * coefsCount, h_targetCoefs_p + s * coefsCount, h_cgesvA + s * N * N, h_cgesvB + s * N );

            /*if (s == 0 && step == 0 && c == 0) {
              magma_cprint(N, 1, h_cgesvB + s * N, N);
              magma_cprint(N, N, h_cgesvA + s * N * N, N);
            }*/

            // -- G-N solve --
            lapackf77_cgesv( &N, &nrhs, h_cgesvA + s * N * N, &N, ipiv + s * N, h_cgesvB + s * N, &N, &locinfo );
            if (locinfo != 0) {
                printf("lapackf77_cgesv matrix %lld returned error %lld: %s.\n", (long long) s, (long long) locinfo, magma_strerror( locinfo ));
            }

            // -- make correction and compute norm --
            for (int i = 0; i < N; i++) {
              (h_Track_cpu + s * N)[i] -= (h_cgesvB + s * N)[i];
              sqrt_sols += MAGMA_C_REAL((h_cgesvB + s * N)[i])*MAGMA_C_REAL((h_cgesvB + s * N)[i]) + MAGMA_C_IMAG((h_cgesvB + s * N)[i])*MAGMA_C_IMAG((h_cgesvB + s * N)[i]);
              sqrt_corr += MAGMA_C_REAL((h_Track_cpu + s * N)[i])*MAGMA_C_REAL((h_Track_cpu + s * N)[i]) + MAGMA_C_IMAG((h_Track_cpu + s * N)[i])*MAGMA_C_IMAG((h_Track_cpu + s * N)[i]);
            }

            // -- check whether the prediction is successful --
            isSuccessful = sqrt_sols < 0.000001 * sqrt_corr;

            // -- a successful prediction stops the correction --
            if (isSuccessful)
               break;
          }

          if ((sqrt_corr > 1e14) && (t0 < 1.0) && (1.0-t0 > 0.001)) {
            //std::cout<<"infinity failed for track "<<s<<std::endl;
            inf_failed = 1;
            cpu_inf_failed_count++;
            break;
          }

          // ===================================================================
          // -- Decide Track Changes --
          // ===================================================================
          if (!isSuccessful) {
            pred_success_count = 0;
            delta_t *= 0.5;
            t0 = t_step;
            // -- should be the last successful tracked sols --
            for (int i = 0; i < N; i++) {
              (h_Track_cpu + s * N)[i] = (h_Track_last_success + s * N)[i];
              (h_Sols_cpu + s * N)[i] = (h_Track_last_success + s * N)[i];
            }
          }
          else {
            for (int i = 0; i < N; i++) {
              (h_Track_last_success + s * N)[i] = (h_Track_cpu + s * N)[i];
              (h_Sols_cpu + s * N)[i] = (h_Track_cpu + s * N)[i];
            }
            pred_success_count++;
            if (pred_success_count >= successes_of_increase_factor) {
              pred_success_count = 0;
              delta_t *= 2;
            }
          }

          // ===================================================================
          // -- Decide delta t at end zone --
          // ===================================================================
          end_zone = (!end_zone && fabs(1 - t0) <= (end_zone_factor + 0.0000001));
          if (end_zone)
              if (delta_t > fabs(1 - t0))
                  delta_t = fabs(1 - t0);
          else if (delta_t > fabs(1 - end_zone_factor - t0))
              delta_t = fabs(1 - end_zone_factor - t0);

          track_steps++;
        }
        else {
          track_success = (t0 > 1.0 || (1.0-t0 < 0.00001));
          h_Track_Success[s] = track_success;
          cpu_track_success_count++;
          break;
        }
      }

      if (track_success) {
        int imag_zero_counter = 0;
        // -- find real tracks --
        for (int j = 0; j < N; j++) {
          if (fabs(MAGMA_C_IMAG((h_Track_cpu + s * N)[j])) < 0.001)
            imag_zero_counter++;
        }

        // -- if the track is a real solution --
        if (imag_zero_counter == N) {
          for(int sol_idx = 0; sol_idx < N; sol_idx++) {
            track_sols_file << MAGMA_C_REAL((h_Track_cpu + s * N)[sol_idx]) << "\n";
          }
          track_sols_file << "\n";
          cpu_hc_steps_file << track_steps <<"\n";
        }

        //for(int i = 0; i < N; i++) {
          //track_sols_file << MAGMA_C_REAL((h_Track_cpu + s * N)[i]) << "+" << MAGMA_C_IMAG((h_Track_cpu + s * N)[i]) << "i"<< "\n";
        //}
        //track_sols_file << "\n";
      }

      /*if (track_success) {
          //std::cout<<"CPU-only computes track "<<s<<":"<<std::endl;
          //magma_cprint(N, 1, h_Track_cpu + s * N, N);
          std::cout<<track_steps<<std::endl;
      }*/
      /*if (s == 0) {
          std::cout<<"CPU-only computes track "<<s<<":"<<std::endl;
          magma_cprint(N, 1, h_Track_cpu + s * N, N);
          std::cout<<s<<"\t"<<track_steps<<"\t"<<track_success<<std::endl;
      }*/
      //std::cout<<s<<"\t"<<track_steps<<"\t"<<track_success<<std::endl;
      //std::cout<<s<<", t0 = "<<t0<<", track_steps = "<<track_steps<<", inf_failed = "<<inf_failed<<", track successes: "<<track_success<<std::endl;
    }
    #if !defined (BATCHED_DISABLE_PARCPU) && defined(_OPENMP)
      openblas_set_num_threads(nthreads);
    #endif
    cpu_time = magma_sync_wtime( my_queue ) - cpu_time;

    std::cout<<"Number of converged solution tracks: "<<cpu_track_success_count<<std::endl;
    std::cout<<"Number of infinity failed tracks: "<<cpu_inf_failed_count<<std::endl;

    track_sols_file.close();

    magma_free_cpu( ipiv );
    magma_free_cpu( h_Track_last_success );
    //magma_free_cpu( h_Track_cpu );
    magma_free_cpu( h_startCoefs_p );
    magma_free_cpu( h_targetCoefs_p );

    return cpu_time;
  }

} // end of namespace

#endif
