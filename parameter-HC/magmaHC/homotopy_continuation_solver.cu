#ifndef homotopy_continuation_solver_cu
#define homotopy_continuation_solver_cu
// =======================================================================
//
// Modifications
//    Chien  21-10-18:   Originally Created
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

#include "magmaHC-problems.cuh"
#include "gpu-kernels/magmaHC-kernels.h"
#include "cpu-compute/cpu-compute.h"
#include "const-matrices/define-const-matrices.h"
#include "const-matrices/const-matrices-allocations.h"

namespace magmaHCWrapper {

  /*void find_real_tracks(
    magmaFloatComplex *h_Track_sols, magmaFloatComplex *h_Track_converge,
    int batchCount, int numOfVars, std::string HC_problem, bool noise
  )
  {
    double err = 0.0;
    int imag_zero_counter = 0;

    double smallest_imag_val_thr = 0.001;
    double smallest_sol_err_thr = 0.0001;

    if (noise) {
      smallest_imag_val_thr = 0.001;
      smallest_sol_err_thr = 0.01;
    }
    else {
      smallest_imag_val_thr = 0.0001;
      smallest_sol_err_thr = 0.00001;
    }

    // -- write HC solutions to the file --
    std::ofstream track_sols_file;
    std::ofstream hc_steps_file;
    track_sols_file.open("/users/cchien3/data/cchien3/MyBitBucket/magmahc-in-3v-unknownf/real_sols.txt");
    hc_steps_file.open("/users/cchien3/data/cchien3/MyBitBucket/magmahc-in-3v-unknownf/steps.txt");
    if ( !track_sols_file.is_open() || !hc_steps_file.is_open() )
       std::cout<<"files cannot be opened!"<<std::endl;

    if (HC_problem == "3view_unknownf_pHC") {
      std::cout<<"entering ..."<<std::endl;
      for (int i = 0; i < batchCount; i++) {
        err = 0.0;
        imag_zero_counter = 0;
        // -- find succssful converged HC tracks --
        if (MAGMA_C_REAL((h_Track_converge + i * numOfVars)[1]) == 1.0) {

          // -- find real tracks --
          for (int j = 0; j < numOfVars; j++) {
            if (fabs(MAGMA_C_IMAG((h_Track_sols + i * (numOfVars+1))[j])) < smallest_imag_val_thr)
              imag_zero_counter++;
          }

          // -- if the track is a real solution --
          if (imag_zero_counter == numOfVars) {
            for(int sol_idx = 0; sol_idx < numOfVars; sol_idx++) {
              track_sols_file << MAGMA_C_REAL((h_Track_sols + i * (numOfVars+1))[sol_idx]) << "\n";
              //track_sols_file << MAGMA_C_REAL((h_Track_sols + i * (numOfVars+1))[sol_idx]) << "\t" << MAGMA_C_IMAG((h_Track_sols + i * (numOfVars+1))[sol_idx]) << "\n";
            }
            track_sols_file << "\n";
            hc_steps_file << MAGMA_C_IMAG((h_Track_converge + i * numOfVars)[1]) <<"\n";
          }
        }
      }
    }
    else {
      std::cout<<"not doing the find_real_tracks operation!"<<std::endl;
    }

    track_sols_file.close();
    hc_steps_file.close();
  }*/

  void homotopy_continuation_solver(
    magmaFloatComplex *h_startSols, magmaFloatComplex *h_Track,
    magmaFloatComplex *h_startParams, magmaFloatComplex *h_targetParams,
    magma_int_t *h_Hx_idx, magma_int_t *h_Ht_idx,
    magmaFloatComplex *h_phc_coeffs_H, magmaFloatComplex *h_phc_coeffs_Ht,
    problem_params* pp, const_mats *cm, std::string hc_problem, std::ofstream &track_sols_file)
  {
    magma_init();
    magma_print_environment();

    magma_int_t batchCount = pp->numOfTracks;
    magma_int_t coefsCount = pp->numOfParams;
    magma_int_t N = pp->numOfVars;

    real_Double_t     gpu_time;
    real_Double_t     cpu_time;
    real_Double_t     data_h2d_time, data_d2h_time;
    magmaFloatComplex *h_cgesvA, *h_cgesvB;
    magmaFloatComplex *h_cgesvA_verify, *h_cgesvB_verify;
    magmaFloatComplex *h_track_sols, *h_track_numOfPred_s;
    magmaFloatComplex_ptr d_startSols, d_Track;
    magmaFloatComplex_ptr d_startCoefs, d_targetCoefs;
    magmaFloatComplex_ptr d_cgesvA, d_cgesvB;
    magma_int_t lda, ldb, ldda, lddb, ldd_coefs, sizeA, sizeB;
    magma_int_t ldd_const_matrices_Hx, ldd_const_matrices_Ht;
    magma_int_t ldd_const_matrices_Hx_collection, ldd_const_matrices_Ht_collection;
    magma_int_t ione = 1;
    magma_int_t ISEED[4] = {0,0,0,1};

    std::cout<<"check point 0"<<std::endl;

    magmaFloatComplex **d_startSols_array = NULL;
    magmaFloatComplex **d_Track_array = NULL;
    magmaFloatComplex **d_startCoefs_array = NULL;
    magmaFloatComplex **d_targetCoefs_array = NULL;
    magmaFloatComplex **d_cgesvA_array = NULL;
    magmaFloatComplex **d_cgesvB_array = NULL;

    magma_device_t cdev;       // variable to indicate current gpu id
    magma_queue_t my_queue;    // magma queue variable, internally holds a cuda stream and a cublas handle
    magma_getdevice( &cdev );
    magma_queue_create( cdev, &my_queue );     // create a queue on this cdev

    std::cout<<"check point 1"<<std::endl;

    lda    = N;
    ldb    = lda;
    ldda   = magma_roundup( N, 32 );  // multiple of 32 by default
    lddb   = ldda;
    ldd_coefs = magma_roundup( coefsCount, 32 );  // multiple of 32 by default
    sizeA = lda*N*batchCount;
    sizeB = ldb*batchCount;

    std::cout<<"check point 2"<<std::endl;

    // ==================================================================================================
    // -- Hx and Ht index matrices --
    // ==================================================================================================
    magmaInt_ptr d_Hx_idx;
    magmaInt_ptr d_Ht_idx;
    magma_int_t **d_Hx_idx_array = NULL;
    magma_int_t **d_Ht_idx_array = NULL;
    magma_int_t size_Hx = N*N*pp->Hx_maximal_terms*(pp->Hx_maximal_parts);
    magma_int_t size_Ht = N*pp->Ht_maximal_terms*(pp->Ht_maximal_parts);
    magma_int_t rnded_size_Hx = magma_roundup( size_Hx, 32 );
    magma_int_t rnded_size_Ht = magma_roundup( size_Ht, 32 );

    // -- allocate gpu memories --
    magma_imalloc( &d_Hx_idx, size_Hx );
    magma_imalloc( &d_Ht_idx, size_Ht );
    magma_malloc( (void**) &d_Hx_idx_array, sizeof(magma_int_t*) );
    magma_malloc( (void**) &d_Ht_idx_array, sizeof(magma_int_t*) );

    // -- transfer data from cpu to gpu --
    magma_isetmatrix( size_Hx, 1, h_Hx_idx, size_Hx, d_Hx_idx, rnded_size_Hx, my_queue );
    magma_isetmatrix( size_Ht, 1, h_Ht_idx, size_Ht, d_Ht_idx, rnded_size_Ht, my_queue );
    magma_iset_pointer( d_Hx_idx_array, d_Hx_idx, rnded_size_Hx, 0, 0, rnded_size_Hx, 1, my_queue );
    magma_iset_pointer( d_Ht_idx_array, d_Ht_idx, rnded_size_Ht, 0, 0, rnded_size_Ht, 1, my_queue );
    // ==================================================================================================

    // ==================================================================================================
    // -- coefficient p2c of funtion t for Hx and Ht --
    // ==================================================================================================
    magma_int_t phc_coeffs_Hx_size = pp->numOfCoeffsFromParams*(pp->max_orderOf_t+1);
    magma_int_t phc_coeffs_Ht_size = pp->numOfCoeffsFromParams*(pp->max_orderOf_t);
    magma_int_t ldd_phc_coeffs_Hx = magma_roundup( phc_coeffs_Hx_size, 32 );  // multiple of 32 by default
    magma_int_t ldd_phc_coeffs_Ht = magma_roundup( phc_coeffs_Ht_size, 32 );  // multiple of 32 by default
    magmaFloatComplex_ptr d_phc_coeffs_Hx;
    magmaFloatComplex_ptr d_phc_coeffs_Ht;
    // -- allocate GPU memory --
    magma_cmalloc( &d_phc_coeffs_Hx, ldd_phc_coeffs_Hx );
    magma_cmalloc( &d_phc_coeffs_Ht, ldd_phc_coeffs_Ht );
    // -- transfer from CPU to GPU --
    magma_csetmatrix( phc_coeffs_Hx_size, 1, h_phc_coeffs_H, phc_coeffs_Hx_size, d_phc_coeffs_Hx, ldd_phc_coeffs_Hx, my_queue );
    magma_csetmatrix( phc_coeffs_Ht_size, 1, h_phc_coeffs_Ht, phc_coeffs_Ht_size, d_phc_coeffs_Ht, ldd_phc_coeffs_Ht, my_queue );
    // ==================================================================================================

    // -- allocate CPU memory --
    magma_cmalloc_cpu( &h_cgesvA, N*N*batchCount );
    magma_cmalloc_cpu( &h_cgesvB, N*batchCount );
    magma_cmalloc_cpu( &h_cgesvA_verify, N*N*batchCount );
    magma_cmalloc_cpu( &h_cgesvB_verify, N*batchCount );
    magma_cmalloc_cpu( &h_track_sols, (N+1)*batchCount );
    magma_cmalloc_cpu( &h_track_numOfPred_s, (N+1)*batchCount );

    // -- allocate constant matrcies in CPU and GPU memories and define values of constant matrices --
    cm->const_matrices_allocations(pp, h_startParams, h_targetParams, my_queue);

    // -- allocate GPU gm --
    magma_cmalloc( &d_startSols, (N+1)*batchCount );
    magma_cmalloc( &d_Track, (N+1)*batchCount );
    magma_cmalloc( &d_startCoefs, ldd_coefs );
    magma_cmalloc( &d_targetCoefs, ldd_coefs );
    magma_cmalloc( &d_cgesvA, ldda*N*batchCount );
    magma_cmalloc( &d_cgesvB, ldda*batchCount );

    // -- allocate 2d arrays in GPU gm --
    magma_malloc( (void**) &d_startSols_array,  batchCount * sizeof(magmaFloatComplex*) );
    magma_malloc( (void**) &d_Track_array,    batchCount * sizeof(magmaFloatComplex*) );
    magma_malloc( (void**) &d_startCoefs_array,    sizeof(magmaFloatComplex*) );
    magma_malloc( (void**) &d_targetCoefs_array, sizeof(magmaFloatComplex*) );
    magma_malloc( (void**) &d_cgesvA_array,    batchCount * sizeof(magmaFloatComplex*) );
    magma_malloc( (void**) &d_cgesvB_array,    batchCount * sizeof(magmaFloatComplex*) );

    // -- random initialization for h_cgesvA and h_cgesvB (doesn't matter the value) --
    lapackf77_clarnv( &ione, ISEED, &sizeA, h_cgesvA );
    lapackf77_clarnv( &ione, ISEED, &sizeB, h_cgesvB );

    // -- transfer data from CPU memory to GPU memory --
    data_h2d_time = magma_sync_wtime( my_queue );
    magma_csetmatrix( N+1, batchCount, h_startSols, (N+1), d_startSols, (N+1), my_queue );
    magma_csetmatrix( N+1, batchCount, h_Track, (N+1), d_Track, (N+1), my_queue );
    magma_csetmatrix( coefsCount, 1, h_startParams, coefsCount, d_startCoefs, ldd_coefs, my_queue );
    magma_csetmatrix( coefsCount, 1, h_targetParams, coefsCount, d_targetCoefs, ldd_coefs, my_queue );
    magma_csetmatrix( N, N*batchCount, h_cgesvA, lda, d_cgesvA, ldda, my_queue );
    magma_csetmatrix( N, batchCount,   h_cgesvB, ldb, d_cgesvB, lddb, my_queue );

    // -- connect pointer to 2d arrays --
    magma_cset_pointer( d_startSols_array, d_startSols, (N+1), 0, 0, (N+1), batchCount, my_queue );
    magma_cset_pointer( d_Track_array, d_Track, (N+1), 0, 0, (N+1), batchCount, my_queue );
    magma_cset_pointer( d_startCoefs_array, d_startCoefs, ldd_coefs, 0, 0, ldd_coefs, 1, my_queue );
    magma_cset_pointer( d_targetCoefs_array, d_targetCoefs, ldd_coefs, 0, 0, ldd_coefs, 1, my_queue );
    magma_cset_pointer( d_cgesvA_array, d_cgesvA, ldda, 0, 0, ldda*N, batchCount, my_queue );
    magma_cset_pointer( d_cgesvB_array, d_cgesvB, lddb, 0, 0, ldda, batchCount, my_queue );

    data_h2d_time = magma_sync_wtime( my_queue ) - data_h2d_time;
    //std::cout<<"Host to device data transfer time: "<<data_h2d_time*1000<<std::endl;

    int s = 5;
    //std::cout<<"origins"<<std::endl;
    //magma_cprint(N+1, 1, h_startSols + s * (N+1), (N+1));
    //magma_cprint(coefsCount, 1, h_targetParams, ldd_coefs);
    //magma_cprint(coefsCount, 1, cm->h_const_cd, ldd_coefs);
    //magma_cprint(N, 1, h_cgesvB + s * ldb, ldb);

    // ===================================================================
    // magma GPU cgesv batched solver for Homotopy Continuation
    // ===================================================================
    std::cout<<"GPU computing ..."<<std::endl;

    if (hc_problem == "3view_unknownf_pHC") {
      std::cout<<"Solving 3view_unknownf problem using parametric HC ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_3view_unknownf_pHC(N, batchCount, ldda, my_queue, d_startSols_array, d_Track_array,
                                                     d_startCoefs_array, d_targetCoefs_array, d_cgesvA_array, d_cgesvB_array, cm, 
                                                     d_Hx_idx_array, d_Ht_idx_array, d_phc_coeffs_Hx, d_phc_coeffs_Ht,
                                                     pp->numOfCoeffsFromParams);
    }
    else if (hc_problem == "4vTrg") {
      std::cout<<"Solving 4-view triangulation problem using parametric HC ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_4vTrg(N, batchCount, ldda, my_queue, d_startSols_array, d_Track_array, d_cgesvA_array, d_cgesvB_array, 
                                        cm, d_Hx_idx_array, d_Ht_idx_array, d_phc_coeffs_Hx, d_phc_coeffs_Ht, pp->numOfCoeffsFromParams);
    }
    else if (hc_problem == "3vTrg") {
      std::cout<<"Solving 3-view triangulation problem using parametric HC ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_3vTrg(N, batchCount, ldda, my_queue, d_startSols_array, d_Track_array, d_cgesvA_array, d_cgesvB_array, 
                                        cm, d_Hx_idx_array, d_Ht_idx_array, d_phc_coeffs_Hx, d_phc_coeffs_Ht, pp->numOfCoeffsFromParams);
    }
    else if (hc_problem == "3vTrg_relax") {
      std::cout<<"Solving relaxed 3-view triangulation problem using parametric HC ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_3vTrg_relax(N, batchCount, ldda, my_queue, d_startSols_array, d_Track_array, d_cgesvA_array, d_cgesvB_array, 
                                              cm, d_Hx_idx_array, d_Ht_idx_array, d_phc_coeffs_Hx, d_phc_coeffs_Ht, pp->numOfCoeffsFromParams);
    }
    else if (hc_problem == "5pt_rel_pose_w_depth_recon") {
      std::cout<<"Solving 5-points relative pose with depth reconstruction problem using parametric HC ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_5pt_rel_pose_w_depth_recon(N, batchCount, ldda, my_queue, d_startSols_array, d_Track_array, d_cgesvA_array, d_cgesvB_array, 
                                                             cm, d_Hx_idx_array, d_Ht_idx_array, d_phc_coeffs_Hx, d_phc_coeffs_Ht, pp->numOfCoeffsFromParams);
    }
    else if (hc_problem == "optimalPnP_w_quaternion") {
      std::cout<<"Solving optimal PnP with quaternion problem using parametric HC ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_optimalPnP_w_quaternion(N, batchCount, ldda, my_queue, d_startSols_array, d_Track_array, d_cgesvA_array, d_cgesvB_array, 
                                                             cm, d_Hx_idx_array, d_Ht_idx_array, d_phc_coeffs_Hx, d_phc_coeffs_Ht, pp->numOfCoeffsFromParams);
    }
    else if (hc_problem == "3pt_rel_pose_w_homo_constraint") {
      std::cout<<"Solving 3-point relative pose with homography constraint problem using parametric HC ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_3pt_rel_pose_w_homo_constraint(N, batchCount, ldda, my_queue, d_startSols_array, d_Track_array, d_cgesvA_array, d_cgesvB_array, 
                                                                 cm, d_Hx_idx_array, d_Ht_idx_array, d_phc_coeffs_Hx, d_phc_coeffs_Ht, pp->numOfCoeffsFromParams);
    }
    else if (hc_problem == "r6p") {
      std::cout<<"Solving 6-point rolling shutter absolute pose estimation problem using parametric HC ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_r6p(N, batchCount, ldda, my_queue, d_startSols_array, d_Track_array, d_cgesvA_array, d_cgesvB_array, 
                                      cm, d_Hx_idx_array, d_Ht_idx_array, d_phc_coeffs_Hx, d_phc_coeffs_Ht, pp->numOfCoeffsFromParams);
    }
    else if (hc_problem == "refractive_p5p") {
      std::cout<<"Solving refractive P5P problem using parametric HC ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_refractive_p5p(N, batchCount, ldda, my_queue, d_startSols_array, d_Track_array, d_cgesvA_array, d_cgesvB_array, 
                                                 cm, d_Hx_idx_array, d_Ht_idx_array, d_phc_coeffs_Hx, d_phc_coeffs_Ht, pp->numOfCoeffsFromParams);
    }
    else if (hc_problem == "refractive_p6p") {
      std::cout<<"Solving refractive P6P problem using parametric HC ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_refractive_p6p(N, batchCount, ldda, my_queue, d_startSols_array, d_Track_array, d_cgesvA_array, d_cgesvB_array, 
                                                 cm, d_Hx_idx_array, d_Ht_idx_array, d_phc_coeffs_Hx, d_phc_coeffs_Ht, pp->numOfCoeffsFromParams);
    }

    // -- check returns from the kernel --
    data_d2h_time = magma_sync_wtime( my_queue );
    magma_cgetmatrix( (N+1), batchCount, d_Track, (N+1), h_track_sols, (N+1), my_queue );
    magma_cgetmatrix( (N+1), batchCount, d_startSols, (N+1), h_track_numOfPred_s, (N+1), my_queue );
    magma_cgetmatrix( N, batchCount, d_cgesvB, lddb, h_cgesvB_verify, ldb, my_queue );
    magma_cgetmatrix( N, N*batchCount, d_cgesvA, ldda, h_cgesvA_verify, lda, my_queue );
    data_d2h_time = magma_sync_wtime( my_queue ) - data_d2h_time;
    std::cout<<"results:"<<std::endl;
    magma_cprint(N+1, 1, h_track_sols + s * (N+1), N+1);
    //magma_cprint(N+1, 1, h_track_numOfPred_s + s * (N+1), N+1);
    magma_cprint(N, 1, h_cgesvB_verify + s * ldb, ldb);
    magma_cprint(N, N, h_cgesvA_verify + s * lda, lda);

    // -- print out all solution tracks --
    /*for (int ps = 0; ps < batchCount; ps++) {
      magma_cprint(N+1, 1, h_track_sols + ps * (N+1), N+1);
      magma_cprint(N, 1, h_cgesvB_verify + ps * ldb, ldb);
    }*/

    //find_real_tracks(h_track_sols, h_cgesvB_verify, batchCount, N, hc_problem, 0);

    // ===================================================================
    // CPU-HC
    // ===================================================================
    magmaFloatComplex *h_Track_cpu;
    magma_int_t *h_Track_Success;
    magma_cmalloc_cpu( &h_Track_cpu, batchCount*N );
    magma_imalloc_cpu( &h_Track_Success, batchCount );
    for (int i = 0; i < batchCount; i++) {
      h_Track_Success[i] = 0;
    }

    // -- cpu computation on the HC problem --
    if (hc_problem == "3view_unknownf_pHC") {
      cpu_time = cpu_hc_solver_3view_unknownf_pHC(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startParams, h_targetParams, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, 0, my_queue, 291, h_track_numOfPred_s);
    }
    else if (hc_problem == "4vTrg") {
      cpu_time = cpu_hc_solver_4vTrg(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startParams, h_targetParams, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, 0, my_queue, 141, h_track_numOfPred_s);
    }
    else if (hc_problem == "3vTrg") {
      cpu_time = cpu_hc_solver_3vTrg(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startParams, h_targetParams, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, 0, my_queue, 158, h_track_numOfPred_s);
    }
    else if (hc_problem == "3vTrg_relax") {
      cpu_time = cpu_hc_solver_3vTrg_relax(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startParams, h_targetParams, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, 0, my_queue, 124, h_track_numOfPred_s);
    }
    else if (hc_problem == "5pt_rel_pose_w_depth_recon") {
      cpu_time = cpu_hc_solver_5pt_rel_pose_w_depth_recon(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startParams, h_targetParams, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, 0, my_queue, 500, h_track_numOfPred_s);
    }
    else if (hc_problem == "optimalPnP_w_quaternion") {
      cpu_time = cpu_hc_solver_optimalPnP_w_quaternion(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startParams, h_targetParams, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, 0, my_queue, 500, h_track_numOfPred_s);
    }
    else if (hc_problem == "3pt_rel_pose_w_homo_constraint") {
      cpu_time = cpu_hc_solver_3pt_rel_pose_w_homo_constraint(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startParams, h_targetParams, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, 0, my_queue, 500, h_track_numOfPred_s);
    }
    else if (hc_problem == "r6p") {
      cpu_time = cpu_hc_solver_r6p(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startParams, h_targetParams, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, 0, my_queue, 500, h_track_numOfPred_s);
    }
    else if (hc_problem == "refractive_p5p") {
      cpu_time = cpu_hc_solver_refractive_p5p(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startParams, h_targetParams, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, 0, my_queue, 47, h_track_numOfPred_s);
    }
    else if (hc_problem == "refractive_p6p") {
      cpu_time = cpu_hc_solver_refractive_p6p(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startParams, h_targetParams, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, 0, my_queue, 62, h_track_numOfPred_s);
    }

    printf("%% Problem   # of eqs   # of sols    GPU+CPU time (ms)    CPU time (ms)\n");
    printf("%%===========================================================================\n");
    printf("   %s     %5lld   %10lld           %7.2f            %7.2f\n",
            hc_problem.c_str(), (long long) N, (long long) batchCount, (gpu_time)*1000, cpu_time*1000);
    
    // -- write converged HC tracks to a file --
    for (int sol_idx = 0; sol_idx < batchCount; sol_idx++) {
      if (h_Track_Success[sol_idx] == 1) {
        for(int i = 0; i < N; i++) {
          track_sols_file << MAGMA_C_REAL((h_Track_cpu + sol_idx * N)[i]) << "+" << MAGMA_C_IMAG((h_Track_cpu + sol_idx * N)[i]) << "i"<< "\n";
        }
      }
    }

    magma_queue_destroy( my_queue );

    magma_free_cpu( h_cgesvA );
    magma_free_cpu( h_cgesvB );
    magma_free_cpu( h_cgesvA_verify );
    magma_free_cpu( h_cgesvB_verify );
    magma_free_cpu( h_track_sols );
    magma_free_cpu( h_track_numOfPred_s );
    magma_free_cpu( h_Track_cpu );

    magma_free( d_startSols );
    magma_free( d_Track );
    magma_free( d_startCoefs );
    magma_free( d_targetCoefs );
    magma_free( d_cgesvA );
    magma_free( d_cgesvB );
    magma_free( d_startSols_array );
    magma_free( d_Track_array );
    magma_free( d_startCoefs_array );
    magma_free( d_targetCoefs_array );
    magma_free( d_cgesvA_array );
    magma_free( d_cgesvB_array );

    magma_free(d_phc_coeffs_Hx);
    magma_free(d_phc_coeffs_Ht);

    cm->free_const_matrices();

    fflush( stdout );
    printf( "\n" );
    magma_finalize();
  }

} // end of namespace

#endif
