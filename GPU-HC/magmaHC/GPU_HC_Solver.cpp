#ifndef GPU_HC_Solver_cpp
#define GPU_HC_Solver_cpp
// =============================================================================================================================
//
// ChangLogs
//    22-10-18:   Initially Created (Copied from other repos)
//    23-12-28:   Use macros and organize this file as definitions of GPU_HC_Solver class functions
//    23-12-29:   Change the file name to GPU_HC_Solver.cpp as a pool of defining member functions in class GPU_HC_Solver.hpp
//
//> (c) LEMS, Brown University
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// =============================================================================================================================
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <chrono>
#include <vector>

#include <cuda.h>
#include <cuda_runtime.h>

//> MAGMA
#include "magma_v2.h"
#include "magma_lapack.h"
#include "magma_internal.h"

#include "GPU_HC_Solver.hpp"
#include "definitions.hpp"
#include "gpu-kernels/magmaHC-kernels.hpp"


//> Constructor
GPU_HC_Solver::GPU_HC_Solver( magmaFloatComplex *h_startSols,   magmaFloatComplex *h_Track, \
                              magmaFloatComplex *h_startParams, magmaFloatComplex *h_targetParams, \
                              magma_int_t       *h_Hx_idx,      magma_int_t       *h_Ht_idx, \
                              magmaFloatComplex *h_phc_coeffs_Hx, magmaFloatComplex *h_phc_coeffs_Ht ): \
    h_startSols(h_startSols), h_Track(h_Track), h_startParams(h_startParams), h_targetParams(h_targetParams), \
    h_Hx_idx(h_Hx_idx), h_Ht_idx(h_Ht_idx), h_phc_coeffs_Hx(h_phc_coeffs_Hx), h_phc_coeffs_Ht(h_phc_coeffs_Ht) {

    magma_init();
    magma_print_environment();

    magma_getdevice( &cdev );
    magma_queue_create( cdev, &my_queue );     // create a queue on this cdev

    //> Define the array sizes
    ldda                = magma_roundup( NUM_OF_VARS, 32 );     //> multiple of 32 by default
    lddb                = ldda;
    ldd_params          = magma_roundup( NUM_OF_PARAMS, 32 );   //> multiple of 32 by default
    size_Hx             = NUM_OF_VARS*NUM_OF_VARS*HX_MAXIMAL_TERMS*HX_MAXIMAL_PARTS;
    size_Ht             = NUM_OF_VARS*HT_MAXIMAL_TERMS*HT_MAXIMAL_PARTS;
    phc_coeffs_Hx_size  = (NUM_OF_COEFFS_FROM_PARAMS+1)*(MAX_ORDER_OF_T+1);
    phc_coeffs_Ht_size  = (NUM_OF_COEFFS_FROM_PARAMS+1)*(MAX_ORDER_OF_T);
    ldd_phc_Params_Hx   = magma_roundup( phc_coeffs_Hx_size, 32 );  // multiple of 32 by default
    ldd_phc_Params_Ht   = magma_roundup( phc_coeffs_Ht_size, 32 );  // multiple of 32 by default

    magmaFloatComplex **d_startSols_array = NULL;
    magmaFloatComplex **d_Track_array     = NULL;

    h_is_GPU_HC_Sol_Converge = new bool[ NUM_OF_TRACKS ];
    h_is_GPU_HC_Sol_Infinity = new bool[ NUM_OF_TRACKS ];
}

// void GPU_HC_Solver::Prepare_Files_for_Write() {
    

// }

void GPU_HC_Solver::Allocate_Arrays() {
    //> CPU Allocations
    magma_cmalloc_cpu( &h_params_diff,        NUM_OF_PARAMS+1 );
    magma_cmalloc_cpu( &h_GPU_HC_Track_Sols,  (NUM_OF_VARS+1)*NUM_OF_TRACKS );

    //> GPU Allocations
    magma_imalloc( &d_Hx_idx,             size_Hx );
    magma_imalloc( &d_Ht_idx,             size_Ht );
    magma_cmalloc( &d_phc_coeffs_Hx,      ldd_phc_Params_Hx );
    magma_cmalloc( &d_phc_coeffs_Ht,      ldd_phc_Params_Ht );
    magma_cmalloc( &d_diffParams,         ldd_params );
    magma_cmalloc( &d_startSols,          (NUM_OF_VARS+1)*NUM_OF_TRACKS );
    magma_cmalloc( &d_Track,              (NUM_OF_VARS+1)*NUM_OF_TRACKS );
    magma_cmalloc( &d_startParams,        ldd_params );
    magma_cmalloc( &d_targetParams,       ldd_params );
    cudacheck( cudaMalloc( &d_is_GPU_HC_Sol_Converge, NUM_OF_TRACKS * sizeof(bool) ));
    cudacheck( cudaMalloc( &d_is_GPU_HC_Sol_Infinity, NUM_OF_TRACKS * sizeof(bool) ));

    magma_malloc( (void**) &d_startSols_array,  (NUM_OF_TRACKS) * sizeof(magmaFloatComplex*) );
    magma_malloc( (void**) &d_Track_array,      (NUM_OF_TRACKS) * sizeof(magmaFloatComplex*) );

    magma_cmalloc_cpu( &h_Debug_Purpose,  NUM_OF_TRACKS );
    magma_cmalloc( &d_Debug_Purpose,      NUM_OF_TRACKS );
}

void GPU_HC_Solver::Data_Transfer_From_Host_To_Device() {

    //> Compute difference of start and target parameters
    //  (This is used only in the trifocal relative pose problem, i.e., TRIFOCAL_2OP1P_30X30)
    for(int i = 0; i < NUM_OF_PARAMS; i++) (h_params_diff)[i] = (h_targetParams)[i] - (h_startParams)[i];
    h_params_diff[NUM_OF_PARAMS] = MAGMA_C_ZERO;

    data_h2d_time = magma_sync_wtime( my_queue );

    magma_isetmatrix( size_Hx,            (1), h_Hx_idx,        size_Hx,            d_Hx_idx,        size_Hx,           my_queue );
    magma_isetmatrix( size_Ht,            (1), h_Ht_idx,        size_Ht,            d_Ht_idx,        size_Ht,           my_queue );
    magma_csetmatrix( phc_coeffs_Hx_size, (1), h_phc_coeffs_Hx, phc_coeffs_Hx_size, d_phc_coeffs_Hx, ldd_phc_Params_Hx, my_queue );
    magma_csetmatrix( phc_coeffs_Ht_size, (1), h_phc_coeffs_Ht, phc_coeffs_Ht_size, d_phc_coeffs_Ht, ldd_phc_Params_Ht, my_queue );
    magma_csetmatrix( NUM_OF_PARAMS+1,    (1), h_params_diff,   NUM_OF_PARAMS+1,    d_diffParams,    ldd_params,        my_queue );

    magma_csetmatrix( NUM_OF_VARS+1, NUM_OF_TRACKS, h_startSols,    (NUM_OF_VARS+1), d_startSols,    (NUM_OF_VARS+1), my_queue );
    magma_csetmatrix( NUM_OF_VARS+1, NUM_OF_TRACKS, h_Track,        (NUM_OF_VARS+1), d_Track,        (NUM_OF_VARS+1), my_queue );
    magma_csetmatrix( NUM_OF_PARAMS, (1),           h_startParams,  NUM_OF_PARAMS,   d_startParams,  ldd_params,      my_queue );
    magma_csetmatrix( NUM_OF_PARAMS, (1),           h_targetParams, NUM_OF_PARAMS,   d_targetParams, ldd_params,      my_queue );
    
    //> connect pointer to 2d arrays
    magma_cset_pointer( d_startSols_array, d_startSols, (NUM_OF_VARS+1), 0, 0, (NUM_OF_VARS+1), NUM_OF_TRACKS, my_queue );
    magma_cset_pointer( d_Track_array,     d_Track,     (NUM_OF_VARS+1), 0, 0, (NUM_OF_VARS+1), NUM_OF_TRACKS, my_queue );
    data_h2d_time = magma_sync_wtime( my_queue ) - data_h2d_time;
}

void GPU_HC_Solver::Solve_by_GPU_HC() {
    std::cout << "GPU computing ..." << std::endl;

#if REL_POSE_5PT_GEO_FORM_QUAT
    gpu_time = kernel_HC_Solver_5pt_rel_pos_geo_form_quat
               (my_queue, d_startSols_array, d_Track_array, \
                d_Hx_idx, d_Ht_idx, d_phc_coeffs_Hx, d_phc_coeffs_Ht, \
                d_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Infinity, \
                d_Debug_Purpose);

#elif REL_POSE_5PT_ALG_FORM_QUAT
    gpu_time = kernel_HC_Solver_5pt_rel_pos_alg_form_quat
               (my_queue, d_startSols_array, d_Track_array, \
                d_Hx_idx, d_Ht_idx, d_phc_coeffs_Hx, d_phc_coeffs_Ht, \
                d_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Infinity, \
                d_Debug_Purpose);
#elif TRIFOCAL_2OP1P_30X30
    gpu_time = kernel_HC_Solver_trifocal_2op1p_30x30 \
               (my_queue, d_startSols_array,  d_Track_array, \
                d_startParams, d_targetParams, d_diffParams, \
                d_Hx_idx, d_Ht_idx, \
                d_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Infinity, \
                d_Debug_Purpose );
#endif

    //> Check returns from the GPU kernel
    data_d2h_time = magma_sync_wtime( my_queue );
    magma_cgetmatrix( (NUM_OF_VARS+1), NUM_OF_TRACKS, d_Track,  (NUM_OF_VARS+1), h_GPU_HC_Track_Sols,    (NUM_OF_VARS+1), my_queue );
    cudacheck( cudaMemcpy( h_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Converge, NUM_OF_TRACKS*sizeof(bool), cudaMemcpyDeviceToHost) );
    cudacheck( cudaMemcpy( h_is_GPU_HC_Sol_Infinity, d_is_GPU_HC_Sol_Infinity, NUM_OF_TRACKS*sizeof(bool), cudaMemcpyDeviceToHost) );
    data_d2h_time = magma_sync_wtime( my_queue ) - data_d2h_time;

#if GPU_DEBUG
    magma_cgetmatrix( NUM_OF_TRACKS, (1), d_Debug_Purpose, NUM_OF_TRACKS, h_Debug_Purpose, NUM_OF_VARS, my_queue );
#endif
    
    // int num_of_convergence = 0;
    // int num_of_infinity_sols = 0;
    // for (int bs = 0; bs < NUM_OF_TRACKS; bs++) {
    //   track_sols_file << std::setprecision(10);
    //   int num_of_real_sols = 0;

    //   track_sols_file << h_is_GPU_HC_Sol_Converge[ bs ] << "\n";
    //   for (int vs = 0; vs < NUM_OF_VARS; vs++) {
    //     track_sols_file << std::setprecision(20) << MAGMA_C_REAL((h_GPU_HC_Track_Sols + bs * (NUM_OF_VARS+1))[vs]) << "\t" \
    //                     << std::setprecision(20) << MAGMA_C_IMAG((h_GPU_HC_Track_Sols + bs * (NUM_OF_VARS+1))[vs]) << "\n";
    //   }
    //   track_sols_file << "\n";

    //   if ( h_is_GPU_HC_Sol_Converge[ bs ] ) num_of_convergence++;
    //   if ( h_is_GPU_HC_Sol_Infinity[ bs ] ) num_of_infinity_sols++;
    // }
    
    //> Object for the Evaluations class
    Evaluations Evaluate_GPUHC_Sols;
    Evaluate_GPUHC_Sols.Write_Converged_Sols( h_GPU_HC_Track_Sols, h_is_GPU_HC_Sol_Converge );
    Evaluate_GPUHC_Sols.Evaluate_Sols( h_GPU_HC_Track_Sols, h_is_GPU_HC_Sol_Converge, h_is_GPU_HC_Sol_Infinity );

    std::cout << "## Evaluation of GPU-HC Solutions: " << std::endl;
    std::cout << " - Number of Converged Solutions:       " << Evaluate_GPUHC_Sols.Num_Of_Coverged_Sols << std::endl;
    std::cout << " - Number of Real Solutions:            " << Evaluate_GPUHC_Sols.Num_Of_Real_Sols << std::endl;
    std::cout << " - Number of Infinity Failed Solutions: " << Evaluate_GPUHC_Sols.Num_Of_Inf_Sols << std::endl;
    
    printf("## GPU time = %7.2f (ms)\n", (gpu_time)*1000);
    // printf("- Number of convergence: %d\n", num_of_convergence);
    // printf("- Number of infinity failed solutions: %d\n", num_of_infinity_sols);
}

GPU_HC_Solver::~GPU_HC_Solver() {

    magma_queue_destroy( my_queue );

    delete [] h_is_GPU_HC_Sol_Converge;
    delete [] h_is_GPU_HC_Sol_Infinity;

    magma_free_cpu( h_GPU_HC_Track_Sols );
    magma_free_cpu( h_params_diff );
    magma_free( d_diffParams );
    magma_free( d_is_GPU_HC_Sol_Converge );
    magma_free( d_is_GPU_HC_Sol_Infinity );
    
    magma_free( d_startSols );
    magma_free( d_Track );
    magma_free( d_startParams );
    magma_free( d_targetParams );
    magma_free( d_Hx_idx );
    magma_free( d_Ht_idx );
    magma_free( d_phc_coeffs_Hx );
    magma_free( d_phc_coeffs_Ht );
    magma_free_cpu( h_Debug_Purpose );
    magma_free( d_Debug_Purpose );

    fflush( stdout );
    printf( "\n" );
    magma_finalize();

    // //> Close all files
    // track_sols_file.close();
}


#endif
