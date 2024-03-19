#ifndef GPU_HC_Solver_cpp
#define GPU_HC_Solver_cpp
// =============================================================================================================================
//
// ChangLogs
//    22-10-18:   Initially Created (Copied from other repos)
//    23-12-28:   Use macros and organize this file as definitions of GPU_HC_Solver class functions
//    23-12-29:   Change the file name to GPU_HC_Solver.cpp as a pool of defining member functions in class GPU_HC_Solver.hpp
//    24-02-26:   Add Data Reader to clean up the main code.
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

#include "Data_Reader.hpp"
#include "GPU_HC_Solver.hpp"
#include "definitions.hpp"
#include "Evaluations.hpp"
#include "gpu-kernels/magmaHC-kernels.hpp"

//> Parameter Homotopy Coefficients
#include "PHC_Coeffs/p2c-trifocal_2op1p_30x30.h"
#include "PHC_Coeffs/p2c-5pt_rel_pos_alg_form_quat.h"
#include "PHC_Coeffs/p2c-5pt_rel_pos_geo_form_quat.h"

//> Constructor
GPU_HC_Solver::GPU_HC_Solver(std::string Input_Repo_Path): REPO_PATH(Input_Repo_Path) {

    //> Initialization
    magma_init();
    magma_print_environment();

    magma_getdevice( &cdev );
    magma_queue_create( cdev, &my_queue );

    //> Define the array sizes
    ldda                  = magma_roundup( NUM_OF_VARS, 32 );     //> multiple of 32 by default
    lddb                  = ldda;
    ldd_params            = magma_roundup( NUM_OF_PARAMS, 32 );   //> multiple of 32 by default
    dHdx_Index_Size       = NUM_OF_VARS*NUM_OF_VARS*HX_MAXIMAL_TERMS*HX_MAXIMAL_PARTS;
    dHdt_Index_Size       = NUM_OF_VARS*HT_MAXIMAL_TERMS*HT_MAXIMAL_PARTS;
    dHdx_PHC_Coeffs_Size  = (NUM_OF_COEFFS_FROM_PARAMS+1)*(MAX_ORDER_OF_T+1);
    dHdt_PHC_Coeffs_Size  = (NUM_OF_COEFFS_FROM_PARAMS+1)*(MAX_ORDER_OF_T);
    ldd_phc_Params_Hx     = magma_roundup( dHdx_PHC_Coeffs_Size, 32 );  // multiple of 32 by default
    ldd_phc_Params_Ht     = magma_roundup( dHdt_PHC_Coeffs_Size, 32 );  // multiple of 32 by default

    printf("dHdx_Index_Size      = %5.2f KB\n", (double)(dHdx_Index_Size     *sizeof(magma_int_t))       / 1024.);
    printf("dHdt_Index_Size      = %5.2f KB\n", (double)(dHdx_Index_Size     *sizeof(magma_int_t))       / 1024.);
    printf("dHdx_PHC_Coeffs_Size = %5.2f KB\n", (double)(dHdx_PHC_Coeffs_Size*sizeof(magmaFloatComplex)) / 1024.);
    printf("dHdt_PHC_Coeffs_Size = %5.2f KB\n", (double)(dHdt_PHC_Coeffs_Size*sizeof(magmaFloatComplex)) / 1024.);

    magmaFloatComplex **d_Start_Sols_array      = NULL;
    magmaFloatComplex **d_Homotopy_Sols_array   = NULL;

    //> Define problem file path for problem data reader and output file path for results evaluations
    Problem_File_Path = REPO_PATH + std::string("problems/") + HC_PROBLEM;
    Write_Files_Path = REPO_PATH + WRITE_FILES_FOLDER;
}

void GPU_HC_Solver::Allocate_Arrays() {
    //> CPU Allocations
    magma_cmalloc_cpu( &h_Start_Sols,           NUM_OF_TRACKS*(NUM_OF_VARS+1) );
    magma_cmalloc_cpu( &h_Homotopy_Sols,        NUM_OF_TRACKS*(NUM_OF_VARS+1) );
    magma_cmalloc_cpu( &h_Start_Params,         NUM_OF_PARAMS );
    magma_cmalloc_cpu( &h_Target_Params,        NUM_OF_PARAMS );
    magma_cmalloc_cpu( &h_dHdx_PHC_Coeffs,      dHdx_PHC_Coeffs_Size );
    magma_cmalloc_cpu( &h_dHdt_PHC_Coeffs,      dHdt_PHC_Coeffs_Size );
    magma_imalloc_cpu( &h_dHdx_Index,           dHdx_Index_Size );
    magma_imalloc_cpu( &h_dHdt_Index,           dHdt_Index_Size );
    magma_cmalloc_cpu( &h_diffParams,           NUM_OF_PARAMS+1 );
    magma_cmalloc_cpu( &h_GPU_HC_Track_Sols,    (NUM_OF_VARS+1)*NUM_OF_TRACKS );
    magma_cmalloc_cpu( &h_Debug_Purpose,        NUM_OF_TRACKS );

    h_is_GPU_HC_Sol_Converge = new bool[ NUM_OF_TRACKS ];
    h_is_GPU_HC_Sol_Infinity = new bool[ NUM_OF_TRACKS ];

    //> GPU Allocations
    magma_cmalloc( &d_Start_Sols,               NUM_OF_TRACKS*(NUM_OF_VARS+1) );
    magma_cmalloc( &d_Homotopy_Sols,            NUM_OF_TRACKS*(NUM_OF_VARS+1) );
    magma_cmalloc( &d_Start_Params,             ldd_params );
    magma_cmalloc( &d_Target_Params,            ldd_params );
    magma_cmalloc( &d_dHdx_PHC_Coeffs,          ldd_phc_Params_Hx );
    magma_cmalloc( &d_dHdt_PHC_Coeffs,          ldd_phc_Params_Ht );
    magma_imalloc( &d_dHdx_Index,               dHdx_Index_Size );
    magma_imalloc( &d_dHdt_Index,               dHdt_Index_Size );
    magma_cmalloc( &d_diffParams,               ldd_params );
    magma_cmalloc( &d_Debug_Purpose,            NUM_OF_TRACKS );

    magma_malloc( (void**) &d_Start_Sols_array,     (NUM_OF_TRACKS) * sizeof(magmaFloatComplex*) );
    magma_malloc( (void**) &d_Homotopy_Sols_array,  (NUM_OF_TRACKS) * sizeof(magmaFloatComplex*) );

    cudacheck( cudaMalloc( &d_is_GPU_HC_Sol_Converge, NUM_OF_TRACKS * sizeof(bool) ));
    cudacheck( cudaMalloc( &d_is_GPU_HC_Sol_Infinity, NUM_OF_TRACKS * sizeof(bool) ));
}

bool GPU_HC_Solver::Read_Problem_Data() {

    //> Load problem data to arrays
    bool is_Data_Read_Successfully = false;
    Data_Reader Load_Problem_Data(Problem_File_Path);

    //> (1) Start parameters
    is_Data_Read_Successfully = Load_Problem_Data.Read_Start_Params( h_Start_Params );
    if (!is_Data_Read_Successfully) { LOG_DATA_LOAD_ERROR("Start Parameters"); return false; }

    //> (2) Target parameters
    is_Data_Read_Successfully = Load_Problem_Data.Read_Target_Params( h_Target_Params );
    if (!is_Data_Read_Successfully) { LOG_DATA_LOAD_ERROR("Target Parameters"); return false; }

    //> (3) Start solutions
    is_Data_Read_Successfully = Load_Problem_Data.Read_Start_Sols( h_Start_Sols, h_Homotopy_Sols );
    if (!is_Data_Read_Successfully) { LOG_DATA_LOAD_ERROR("Start Solutions"); return false; }

    //> (4) dH/dx evaluation indices
    is_Data_Read_Successfully = Load_Problem_Data.Read_dHdx_Indices( h_dHdx_Index );
    if (!is_Data_Read_Successfully) { LOG_DATA_LOAD_ERROR("dH/dx Evaluation Indices"); return false; }

    //> (5) dH/dt evaluation indices
    is_Data_Read_Successfully = Load_Problem_Data.Read_dHdt_Indices( h_dHdt_Index );
    if (!is_Data_Read_Successfully) { LOG_DATA_LOAD_ERROR("dH/dt Evaluation Indices"); return false; }

    return true;
}

void GPU_HC_Solver::Construct_Coeffs_From_Params() {

#if REL_POSE_5PT_GEO_FORM_QUAT
    magmaHCWrapper::p2c_5pt_rel_pos_geo_form_quat(h_Target_Params, h_Start_Params, h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs);
#elif REL_POSE_5PT_ALG_FORM_QUAT
    magmaHCWrapper::p2c_5pt_rel_pos_alg_form_quat(h_Target_Params, h_Start_Params, h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs);
#elif TRIFOCAL_2OP1P_30X30
    magmaHCWrapper::p2c_trifocal_2op1p_30x30(h_Target_Params, h_Start_Params, h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs);
#endif
}

void GPU_HC_Solver::Data_Transfer_From_Host_To_Device() {

    //> Compute difference of start and target parameters
    //  (This is used only in the trifocal relative pose problem, i.e., TRIFOCAL_2OP1P_30X30)
    for(int i = 0; i < NUM_OF_PARAMS; i++) (h_diffParams)[i] = (h_Target_Params)[i] - (h_Start_Params)[i];
    h_diffParams[NUM_OF_PARAMS] = MAGMA_C_ZERO;

    transfer_h2d_time = magma_sync_wtime( my_queue );
    magma_csetmatrix( NUM_OF_VARS+1,        NUM_OF_TRACKS, h_Start_Sols,        (NUM_OF_VARS+1),      d_Start_Sols,       (NUM_OF_VARS+1),   my_queue );
    magma_csetmatrix( NUM_OF_VARS+1,        NUM_OF_TRACKS, h_Homotopy_Sols,     (NUM_OF_VARS+1),      d_Homotopy_Sols,    (NUM_OF_VARS+1),   my_queue );
    magma_isetmatrix( dHdx_Index_Size,      (1),           h_dHdx_Index,        dHdx_Index_Size,      d_dHdx_Index,       dHdx_Index_Size,   my_queue );
    magma_isetmatrix( dHdt_Index_Size,      (1),           h_dHdt_Index,        dHdt_Index_Size,      d_dHdt_Index,       dHdt_Index_Size,   my_queue );
    magma_csetmatrix( dHdx_PHC_Coeffs_Size, (1),           h_dHdx_PHC_Coeffs,   dHdx_PHC_Coeffs_Size, d_dHdx_PHC_Coeffs,  ldd_phc_Params_Hx, my_queue );
    magma_csetmatrix( dHdt_PHC_Coeffs_Size, (1),           h_dHdt_PHC_Coeffs,   dHdt_PHC_Coeffs_Size, d_dHdt_PHC_Coeffs,  ldd_phc_Params_Ht, my_queue );
    magma_csetmatrix( NUM_OF_PARAMS+1,      (1),           h_diffParams,        NUM_OF_PARAMS+1,      d_diffParams,       ldd_params,        my_queue );
    magma_csetmatrix( NUM_OF_PARAMS,        (1),           h_Start_Params,      NUM_OF_PARAMS,        d_Start_Params,     ldd_params,        my_queue );
    magma_csetmatrix( NUM_OF_PARAMS,        (1),           h_Target_Params,     NUM_OF_PARAMS,        d_Target_Params,    ldd_params,        my_queue );

    //> connect pointer to 2d arrays
    magma_cset_pointer( d_Start_Sols_array, d_Start_Sols, (NUM_OF_VARS+1), 0, 0, (NUM_OF_VARS+1), NUM_OF_TRACKS, my_queue );
    magma_cset_pointer( d_Homotopy_Sols_array,     d_Homotopy_Sols,     (NUM_OF_VARS+1), 0, 0, (NUM_OF_VARS+1), NUM_OF_TRACKS, my_queue );
    transfer_h2d_time = magma_sync_wtime( my_queue ) - transfer_h2d_time;
}

void GPU_HC_Solver::Solve_by_GPU_HC() {
    std::cout << "GPU computing ..." << std::endl << std::endl;

#if REL_POSE_5PT_GEO_FORM_QUAT
    gpu_time = kernel_HC_Solver_5pt_rel_pos_geo_form_quat
               (my_queue, d_Start_Sols_array, d_Homotopy_Sols_array, \
                d_dHdx_Index, d_dHdt_Index, d_dHdx_PHC_Coeffs, d_dHdt_PHC_Coeffs, \
                d_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Infinity, \
                d_Debug_Purpose);

#elif REL_POSE_5PT_ALG_FORM_QUAT
    gpu_time = kernel_HC_Solver_5pt_rel_pos_alg_form_quat
               (my_queue, d_Start_Sols_array, d_Homotopy_Sols_array, \
                d_dHdx_Index, d_dHdt_Index, d_dHdx_PHC_Coeffs, d_dHdt_PHC_Coeffs, \
                d_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Infinity, \
                d_Debug_Purpose);
#elif TRIFOCAL_2OP1P_30X30
    gpu_time = kernel_HC_Solver_trifocal_2op1p_30x30 \
               (my_queue, d_Start_Sols_array,  d_Homotopy_Sols_array, \
                d_Start_Params, d_Target_Params, d_diffParams, \
                d_dHdx_Index, d_dHdt_Index, \
                d_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Infinity, \
                d_Debug_Purpose );
#endif

    //> Check returns from the GPU kernel
    transfer_d2h_time = magma_sync_wtime( my_queue );
    magma_cgetmatrix( (NUM_OF_VARS+1), NUM_OF_TRACKS, d_Homotopy_Sols,  (NUM_OF_VARS+1), h_GPU_HC_Track_Sols,    (NUM_OF_VARS+1), my_queue );
    cudacheck( cudaMemcpy( h_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Converge, NUM_OF_TRACKS*sizeof(bool), cudaMemcpyDeviceToHost) );
    cudacheck( cudaMemcpy( h_is_GPU_HC_Sol_Infinity, d_is_GPU_HC_Sol_Infinity, NUM_OF_TRACKS*sizeof(bool), cudaMemcpyDeviceToHost) );
    transfer_d2h_time = magma_sync_wtime( my_queue ) - transfer_d2h_time;

#if GPU_DEBUG
    magma_cgetmatrix( NUM_OF_TRACKS, (1), d_Debug_Purpose, NUM_OF_TRACKS, h_Debug_Purpose, NUM_OF_VARS, my_queue );
#endif

    //> Object for the Evaluations class
    Evaluations Evaluate_GPUHC_Sols( Write_Files_Path );
    Evaluate_GPUHC_Sols.Write_Converged_Sols( h_GPU_HC_Track_Sols, h_is_GPU_HC_Sol_Converge );
    Evaluate_GPUHC_Sols.Evaluate_Sols( h_GPU_HC_Track_Sols, h_is_GPU_HC_Sol_Converge, h_is_GPU_HC_Sol_Infinity );
    Evaluate_GPUHC_Sols.Find_Unique_Sols( h_GPU_HC_Track_Sols, h_is_GPU_HC_Sol_Converge );

    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << "## Solving " << PRINT_OUT_PROBLEM_NAME << std::endl;

    //> Print out timings
    printf("## Timings:\n");
    printf(" - GPU Computation Time = %7.2f (ms)\n", (gpu_time)*1000);
    printf(" - H2D Transfer Time    = %7.2f (ms)\n", (transfer_h2d_time)*1000);
    printf(" - D2H Transfer Time    = %7.2f (ms)\n", (transfer_d2h_time)*1000);

    //> Print out evaluation results
    std::cout << "## Evaluation of GPU-HC Solutions: "      << std::endl;
    std::cout << " - Number of Converged Solutions:       " << Evaluate_GPUHC_Sols.Num_Of_Coverged_Sols << std::endl;
    std::cout << " - Number of Real Solutions:            " << Evaluate_GPUHC_Sols.Num_Of_Real_Sols << std::endl;
    std::cout << " - Number of Infinity Failed Solutions: " << Evaluate_GPUHC_Sols.Num_Of_Inf_Sols << std::endl;
    std::cout << " - Number of Unique Solutions:          " << Evaluate_GPUHC_Sols.Num_Of_Unique_Sols << std::endl;


}

GPU_HC_Solver::~GPU_HC_Solver() {

    magma_queue_destroy( my_queue );

    delete [] h_is_GPU_HC_Sol_Converge;
    delete [] h_is_GPU_HC_Sol_Infinity;

    magma_free_cpu( h_Start_Sols );
    magma_free_cpu( h_Homotopy_Sols );
    magma_free_cpu( h_Start_Params );
    magma_free_cpu( h_Target_Params );
    magma_free_cpu( h_dHdx_PHC_Coeffs );
    magma_free_cpu( h_dHdt_PHC_Coeffs );
    magma_free_cpu( h_dHdx_Index );
    magma_free_cpu( h_dHdt_Index );

    magma_free_cpu( h_GPU_HC_Track_Sols );
    magma_free_cpu( h_Debug_Purpose );
    magma_free_cpu( h_diffParams );

    magma_free( d_diffParams );
    magma_free( d_is_GPU_HC_Sol_Converge );
    magma_free( d_is_GPU_HC_Sol_Infinity );

    magma_free( d_Start_Sols );
    magma_free( d_Homotopy_Sols );
    magma_free( d_Start_Params );
    magma_free( d_Target_Params );
    magma_free( d_dHdx_Index );
    magma_free( d_dHdt_Index );
    magma_free( d_dHdx_PHC_Coeffs );
    magma_free( d_dHdt_PHC_Coeffs );
    magma_free( d_Debug_Purpose );

    fflush( stdout );
    printf( "\n" );
    magma_finalize();
}


#endif
