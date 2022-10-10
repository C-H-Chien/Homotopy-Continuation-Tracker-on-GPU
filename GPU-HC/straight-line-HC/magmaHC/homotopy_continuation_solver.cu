#ifndef homotopy_continuation_solver_cu
#define homotopy_continuation_solver_cu
// =======================================================================
//
// Modifications
//    Chien  21-05-06:   Originally Created
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

  void homotopy_continuation_solver(
    magmaFloatComplex *h_startSols, magmaFloatComplex *h_Track,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magma_int_t *h_Hx_idx, magma_int_t *h_Ht_idx,
    problem_params* pp, const_mats *cm, std::string hc_problem,
    std::ofstream &track_sols_file, std::ofstream &tracks_success_file)
  {
    magma_int_t batchCount = pp->numOfTracks;
    magma_int_t coefsCount = pp->numOfCoeffs;
    magma_int_t N = pp->numOfVars;
    
    magma_init();
    magma_print_environment();

    real_Double_t     gpu_time, cpu_continual_time;
    real_Double_t     cpu_time;
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

    lda    = N;
    ldb    = lda;
    ldda   = magma_roundup( N, 32 );  // multiple of 32 by default
    lddb   = ldda;
    ldd_coefs = magma_roundup( coefsCount, 32 );  // multiple of 32 by default
    ldd_const_matrices_Hx = magma_roundup( N*N*pp->Hx_maximal_terms, 32 );  // multiple of 32 by default
    ldd_const_matrices_Ht = magma_roundup( N*pp->Ht_maximal_terms, 32 );  // multiple of 32 by default
    ldd_const_matrices_Hx_collection = magma_roundup( N*N*pp->Hx_maximal_terms*pp->Hx_maximal_parts, 32 );  // multiple of 32 by default
    ldd_const_matrices_Ht_collection = magma_roundup( N*pp->Ht_maximal_terms*pp->Ht_maximal_parts, 32 );  // multiple of 32 by default
    sizeA = lda*N*batchCount;
    sizeB = ldb*batchCount;
    
    // ==================================================================================================
    // -- Hx and Ht index matrices --
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

    // -- allocate CPU memory --
    magma_cmalloc_cpu( &h_cgesvA, N*N*batchCount );
    magma_cmalloc_cpu( &h_cgesvB, N*batchCount );
    magma_cmalloc_cpu( &h_cgesvA_verify, N*N*batchCount );
    magma_cmalloc_cpu( &h_cgesvB_verify, N*batchCount );
    magma_cmalloc_cpu( &h_track_sols, (N+1)*batchCount );
    magma_cmalloc_cpu( &h_track_numOfPred_s, (N+1)*batchCount );

    // -- allocate constant matrcies in CPU and GPU memories and define values of constant matrices --
    cm->const_matrices_allocations(pp, ldd_const_matrices_Hx_collection, ldd_const_matrices_Ht_collection, ldd_const_matrices_Hx, ldd_const_matrices_Ht, h_startCoefs, h_targetCoefs, ldd_coefs, my_queue, hc_problem);

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

    // -- transfer data from CPU memory to GPU memory --
    magma_csetmatrix( N+1, batchCount, h_startSols, (N+1), d_startSols, (N+1), my_queue );
    magma_csetmatrix( N+1, batchCount, h_Track, (N+1), d_Track, (N+1), my_queue );
    magma_csetmatrix( coefsCount, 1, h_startCoefs, coefsCount, d_startCoefs, ldd_coefs, my_queue );
    magma_csetmatrix( coefsCount, 1, h_targetCoefs, coefsCount, d_targetCoefs, ldd_coefs, my_queue );
    magma_csetmatrix( N, N*batchCount, h_cgesvA, lda, d_cgesvA, ldda, my_queue );
    magma_csetmatrix( N, batchCount,   h_cgesvB, ldb, d_cgesvB, lddb, my_queue );

    // -- connect pointer to 2d arrays --
    magma_cset_pointer( d_startSols_array, d_startSols, (N+1), 0, 0, (N+1), batchCount, my_queue );
    magma_cset_pointer( d_Track_array, d_Track, (N+1), 0, 0, (N+1), batchCount, my_queue );
    magma_cset_pointer( d_startCoefs_array, d_startCoefs, ldd_coefs, 0, 0, ldd_coefs, 1, my_queue );
    magma_cset_pointer( d_targetCoefs_array, d_targetCoefs, ldd_coefs, 0, 0, ldd_coefs, 1, my_queue );
    magma_cset_pointer( d_cgesvA_array, d_cgesvA, ldda, 0, 0, ldda*N, batchCount, my_queue );
    magma_cset_pointer( d_cgesvB_array, d_cgesvB, lddb, 0, 0, ldda, batchCount, my_queue );

    std::cout<<"check point 4"<<std::endl;

    int s = 0;
    //std::cout<<"origins"<<std::endl;
    //magma_cprint(N+1, 1, h_startSols + s * (N+1), (N+1));
    //magma_cprint(coefsCount, 1, h_startCoefs, ldd_coefs);
    //magma_cprint(N, 1, h_cgesvB + s * ldb, ldb);

    // ===================================================================
    // magma GPU cgesv batched solver for Homotopy Continuation
    // ===================================================================
    std::cout<<"GPU computing ..."<<std::endl;

    if (hc_problem == "alea6") {
      std::cout<<"Solving alea6 problem ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_alea6(N, batchCount, coefsCount, ldda, my_queue, d_startSols_array, d_Track_array,
                                     d_startCoefs_array, d_targetCoefs_array, d_cgesvA_array, d_cgesvB_array, cm, d_Hx_idx_array, d_Ht_idx_array);
    }
    else if (hc_problem == "eco12") {
      std::cout<<"Solving eco12 problem ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_eco12(N, batchCount, coefsCount, ldda, my_queue, d_startSols_array, d_Track_array,
                                        d_startCoefs_array, d_targetCoefs_array, d_cgesvA_array, d_cgesvB_array, cm, d_Hx_idx_array, d_Ht_idx_array);
    }
    else if (hc_problem == "d1") {
      std::cout<<"Solving d1 problem ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_d1(N, batchCount, coefsCount, ldda, my_queue, d_startSols_array, d_Track_array,
                                     d_startCoefs_array, d_targetCoefs_array, d_cgesvA_array, d_cgesvB_array, cm, d_Hx_idx_array, d_Ht_idx_array);
    }
    else if (hc_problem == "katsura6") {
      std::cout<<"Solving katsura6 problem ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_katsura6(N, batchCount, coefsCount, ldda, my_queue, d_startSols_array, d_Track_array,
                                          d_startCoefs_array, d_targetCoefs_array, d_cgesvA_array, d_cgesvB_array, cm, d_Hx_idx_array, d_Ht_idx_array);
    }
    else if (hc_problem == "katsura7") {
      std::cout<<"Solving katsura7 problem ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_katsura7(N, batchCount, coefsCount, ldda, my_queue, d_startSols_array, d_Track_array,
                                          d_startCoefs_array, d_targetCoefs_array, d_cgesvA_array, d_cgesvB_array, cm, d_Hx_idx_array, d_Ht_idx_array);
    }
    else if (hc_problem == "katsura8") {
      std::cout<<"Solving katsura8 problem ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_katsura8(N, batchCount, coefsCount, ldda, my_queue, d_startSols_array, d_Track_array,
                                          d_startCoefs_array, d_targetCoefs_array, d_cgesvA_array, d_cgesvB_array, cm, d_Hx_idx_array, d_Ht_idx_array);
    }
    else if (hc_problem == "katsura9") {
      std::cout<<"Solving katsura9 problem ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_katsura9(N, batchCount, coefsCount, ldda, my_queue, d_startSols_array, d_Track_array,
                                          d_startCoefs_array, d_targetCoefs_array, d_cgesvA_array, d_cgesvB_array, cm, d_Hx_idx_array, d_Ht_idx_array);
    }
    else if (hc_problem == "katsura10") {
      std::cout<<"Solving katsura10 problem ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_katsura10(N, batchCount, coefsCount, ldda, my_queue, d_startSols_array, d_Track_array,
                                          d_startCoefs_array, d_targetCoefs_array, d_cgesvA_array, d_cgesvB_array, cm, d_Hx_idx_array, d_Ht_idx_array);
    }
    else if (hc_problem == "katsura11") {
      std::cout<<"Solving katsura11 problem ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_katsura11(N, batchCount, coefsCount, ldda, my_queue, d_startSols_array, d_Track_array,
                                          d_startCoefs_array, d_targetCoefs_array, d_cgesvA_array, d_cgesvB_array, cm, d_Hx_idx_array, d_Ht_idx_array);
    }
    else if (hc_problem == "katsura12") {
      std::cout<<"Solving katsura12 problem ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_katsura12(N, batchCount, coefsCount, ldda, my_queue, d_startSols_array, d_Track_array,
                                          d_startCoefs_array, d_targetCoefs_array, d_cgesvA_array, d_cgesvB_array, cm, d_Hx_idx_array, d_Ht_idx_array);
    }
    else if (hc_problem == "katsura13") {
      std::cout<<"Solving katsura13 problem ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_katsura13(N, batchCount, coefsCount, ldda, my_queue, d_startSols_array, d_Track_array,
                                          d_startCoefs_array, d_targetCoefs_array, d_cgesvA_array, d_cgesvB_array, cm, d_Hx_idx_array, d_Ht_idx_array);
    }
    else if (hc_problem == "katsura14") {
      std::cout<<"Solving katsura14 problem ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_katsura14(N, batchCount, coefsCount, ldda, my_queue, d_startSols_array, d_Track_array,
                                          d_startCoefs_array, d_targetCoefs_array, d_cgesvA_array, d_cgesvB_array, cm, d_Hx_idx_array, d_Ht_idx_array);
    }
    else if (hc_problem == "katsura15") {
      std::cout<<"Solving katsura15 problem ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_katsura15(N, batchCount, coefsCount, ldda, my_queue, d_startSols_array, d_Track_array,
                                          d_startCoefs_array, d_targetCoefs_array, d_cgesvA_array, d_cgesvB_array, cm, d_Hx_idx_array, d_Ht_idx_array);
    }
    /*else if (hc_problem == "cyclic7") {
      std::cout<<"Solving cyclic7 problem ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_cyclic7(N, batchCount, coefsCount, ldda, my_queue, d_startSols_array, d_Track_array,
                                          d_startCoefs_array, d_targetCoefs_array, d_cgesvA_array, d_cgesvB_array, cm, d_Hx_idx_array, d_Ht_idx_array);
    }
    else if (hc_problem == "cyclic8") {
      std::cout<<"Solving cyclic8 problem ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_cyclic8(N, batchCount, coefsCount, ldda, my_queue, d_startSols_array, d_Track_array,
                                          d_startCoefs_array, d_targetCoefs_array, d_cgesvA_array, d_cgesvB_array, cm, d_Hx_idx_array, d_Ht_idx_array);
    }
    else if (hc_problem == "cyclic9") {
      std::cout<<"Solving cyclic9 problem ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_cyclic9(N, batchCount, coefsCount, ldda, my_queue, d_startSols_array, d_Track_array,
                                          d_startCoefs_array, d_targetCoefs_array, d_cgesvA_array, d_cgesvB_array, cm, d_Hx_idx_array, d_Ht_idx_array);
    }
    else if (hc_problem == "game6two") {
      std::cout<<"Solving game6two problem ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_game6two(N, batchCount, coefsCount, ldda, my_queue, d_startSols_array, d_Track_array,
                                          d_startCoefs_array, d_targetCoefs_array, d_cgesvA_array, d_cgesvB_array, cm, d_Hx_idx_array, d_Ht_idx_array);
    }
    else if (hc_problem == "game7two") {
      std::cout<<"Solving game7two problem ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_game7two(N, batchCount, coefsCount, ldda, my_queue, d_startSols_array, d_Track_array,
                                          d_startCoefs_array, d_targetCoefs_array, d_cgesvA_array, d_cgesvB_array, cm, d_Hx_idx_array, d_Ht_idx_array);
    }
    else if (hc_problem == "pole28sys") {
      std::cout<<"Solving pole28sys problem ..."<<std::endl<<std::endl;
      gpu_time = kernel_HC_Solver_pole28sys(N, batchCount, coefsCount, ldda, my_queue, d_startSols_array, d_Track_array,
                                          d_startCoefs_array, d_targetCoefs_array, d_cgesvA_array, d_cgesvB_array, cm, d_Hx_idx_array, d_Ht_idx_array);
    }*/

    // -- check returns from the kernel --
    //data_d2h_time = magma_sync_wtime( my_queue );
    magma_cgetmatrix( (N+1), batchCount, d_Track, (N+1), h_track_sols, (N+1), my_queue );
    magma_cgetmatrix( (N+1), batchCount, d_startSols, (N+1), h_track_numOfPred_s, (N+1), my_queue );
    magma_cgetmatrix( N, batchCount, d_cgesvB, lddb, h_cgesvB_verify, ldb, my_queue );
    magma_cgetmatrix( N, N*batchCount, d_cgesvA, ldda, h_cgesvA_verify, lda, my_queue );
    //data_d2h_time = magma_sync_wtime( my_queue ) - data_d2h_time;
    std::cout<<"results:"<<std::endl;
    magma_cprint(N+1, 1, h_track_sols + s * (N+1), N+1);
    //magma_cprint(N+1, 1, h_track_numOfPred_s + s * (N+1), N+1);
    magma_cprint(N, 1, h_cgesvB_verify + s * ldb, ldb);
    //magma_cprint(N, N, h_cgesvA_verify + s * lda, lda);

    // ============================================================================================
    // -- CPU-HC --
    // ============================================================================================
    magmaFloatComplex *h_Track_cpu;
    magma_int_t *h_Track_Success;
    magma_cmalloc_cpu( &h_Track_cpu, batchCount*N );
    magma_imalloc_cpu( &h_Track_Success, batchCount );
    for (int i = 0; i < batchCount; i++) {
      h_Track_Success[i] = 0;
    }

    // -- cpu computation on the HC problem --
    std::cout<<"CPU-only computing ..."<<std::endl;
    if (hc_problem == "alea6") {
      cpu_time = cpu_hc_solver_alea6(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startCoefs, h_targetCoefs, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, my_queue, 120, tracks_success_file);
    }
    else if (hc_problem == "eco12") {
      cpu_time = cpu_hc_solver_echo12(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startCoefs, h_targetCoefs, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, my_queue, 70, tracks_success_file);
    }
    else if (hc_problem == "d1") {
      cpu_time = cpu_hc_solver_d1(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startCoefs, h_targetCoefs, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, my_queue, 50, tracks_success_file);
    }
    else if (hc_problem == "katsura6") {
      cpu_time = cpu_hc_solver_katsura6(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startCoefs, h_targetCoefs, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, my_queue, 400, tracks_success_file);
    }
    else if (hc_problem == "katsura7") {
      cpu_time = cpu_hc_solver_katsura7(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startCoefs, h_targetCoefs, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, my_queue, 400, tracks_success_file);
    }
    else if (hc_problem == "katsura8") {
      cpu_time = cpu_hc_solver_katsura8(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startCoefs, h_targetCoefs, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, my_queue, 400, tracks_success_file);
    }
    else if (hc_problem == "katsura9") {
      cpu_time = cpu_hc_solver_katsura9(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startCoefs, h_targetCoefs, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, my_queue, 400, tracks_success_file);
    }
    else if (hc_problem == "katsura10") {
      cpu_time = cpu_hc_solver_katsura10(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startCoefs, h_targetCoefs, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, my_queue, 400, tracks_success_file);
    }
    else if (hc_problem == "katsura11") {
      cpu_time = cpu_hc_solver_katsura11(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startCoefs, h_targetCoefs, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, my_queue, 400, tracks_success_file);
    }
    else if (hc_problem == "katsura12") {
      cpu_time = cpu_hc_solver_katsura12(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startCoefs, h_targetCoefs, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, my_queue, 400, tracks_success_file);
    }
    else if (hc_problem == "katsura13") {
      cpu_time = cpu_hc_solver_katsura13(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startCoefs, h_targetCoefs, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, my_queue, 400, tracks_success_file);
    }
    else if (hc_problem == "katsura14") {
      cpu_time = cpu_hc_solver_katsura14(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startCoefs, h_targetCoefs, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, my_queue, 400, tracks_success_file);
    }
    else if (hc_problem == "katsura15") {
      cpu_time = cpu_hc_solver_katsura15(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startCoefs, h_targetCoefs, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, my_queue, 400, tracks_success_file);
    }
    /*else if (hc_problem == "cyclic7") {
      cpu_time = cpu_hc_solver_cyclic7(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startCoefs, h_targetCoefs, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, my_queue, 90, tracks_success_file);
    }
    else if (hc_problem == "cyclic8") {
      cpu_time = cpu_hc_solver_cyclic8(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startCoefs, h_targetCoefs, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, my_queue, 90, tracks_success_file);
    }
    else if (hc_problem == "cyclic9") {
      cpu_time = cpu_hc_solver_cyclic9(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startCoefs, h_targetCoefs, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, my_queue, 90, tracks_success_file);
    }
    else if (hc_problem == "game6two") {
      cpu_time = cpu_hc_solver_game6two(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startCoefs, h_targetCoefs, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, my_queue, 400, tracks_success_file);
    }
    else if (hc_problem == "game7two") {
      cpu_time = cpu_hc_solver_game7two(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startCoefs, h_targetCoefs, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, my_queue, 400, tracks_success_file);
    }
    else if (hc_problem == "pole28sys") {
      cpu_time = cpu_hc_solver_pole28sys(h_Track_cpu, h_Track_Success, h_startSols, h_Track, h_startCoefs, h_targetCoefs, h_cgesvA, h_cgesvB, batchCount, coefsCount, N, my_queue, 400, tracks_success_file);
    }*/

    printf("%% Problem   # of eqs   # of sols    GPU+CPU time (ms)    CPU time (ms)\n");
    printf("%%===========================================================================\n");
    printf("   %s     %5lld   %10lld           %7.2f            %7.2f\n",
            hc_problem.c_str(), (long long) N, (long long) batchCount, (gpu_time+cpu_continual_time)*1000, cpu_time*1000);

    // -- write converged HC tracks to a file --
    for (int sol_idx = 0; sol_idx < batchCount; sol_idx++) {
      if (h_Track_Success[sol_idx] == 1) {
        for(int i = 0; i < N; i++) {
          track_sols_file << MAGMA_C_REAL((h_Track_cpu + sol_idx * N)[i]) << "+" << MAGMA_C_IMAG((h_Track_cpu + sol_idx * N)[i]) << "i"<< "\n";
        }
        //track_sols_file<<"\n";
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
    magma_free_cpu( h_Track_Success );

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

    cm->free_const_matrices();

    fflush( stdout );
    printf( "\n" );
    magma_finalize();
  }

} // end of namespace

#endif
