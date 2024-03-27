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
#undef max
#undef min

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
GPU_HC_Solver::GPU_HC_Solver(YAML::Node Problem_Settings_File): Problem_Setting_YAML_File(Problem_Settings_File) {

    //> Parse data from the YAML file
    //> (1) Problem Name and GPU-HC Type
    HC_problem                      = Problem_Setting_YAML_File["problem_name"].as<std::string>();
    HC_print_problem_name           = Problem_Setting_YAML_File["problem_print_out_name"].as<std::string>();
    GPUHC_type                      = Problem_Setting_YAML_File["GPUHC_Type"].as<std::string>();
    //> (2) GPU-HC Configurations
    GPUHC_Max_Steps                 = Problem_Setting_YAML_File["GPUHC_Max_Steps"].as<int>();
    GPUHC_Max_Correction_Steps      = Problem_Setting_YAML_File["GPUHC_Max_Correction_Steps"].as<int>();
    GPUHC_delta_t_incremental_steps = Problem_Setting_YAML_File["GPUHC_Num_Of_Steps_to_Increase_Delta_t"].as<int>();
    //> (3) Problem Specifications
    Num_Of_Vars                     = Problem_Setting_YAML_File["Num_Of_Vars"].as<int>();
    Num_Of_Params                   = Problem_Setting_YAML_File["Num_Of_Params"].as<int>();
    Num_Of_Tracks                   = Problem_Setting_YAML_File["Num_Of_Tracks"].as<int>();
    dHdx_Max_Terms                  = Problem_Setting_YAML_File["dHdx_Max_Terms"].as<int>();
    dHdx_Max_Parts                  = Problem_Setting_YAML_File["dHdx_Max_Parts"].as<int>();
    dHdt_Max_Terms                  = Problem_Setting_YAML_File["dHdt_Max_Terms"].as<int>();
    dHdt_Max_Parts                  = Problem_Setting_YAML_File["dHdt_Max_Parts"].as<int>();
    Max_Order_Of_T                  = Problem_Setting_YAML_File["Max_Order_Of_T"].as<int>();

    //> Check if conversion from problem parameters to polynomial coefficients is needed
    if (GPUHC_type == std::string("P2C")) {
        Use_P2C = true;
        Num_Of_Coeffs_From_Params   = Problem_Setting_YAML_File["Num_Of_Coeffs_From_Params"].as<int>();
    }
    else Use_P2C = false;

    //> Initialization
    magma_init();
    magma_print_environment();

    magma_getdevice( &cdev );
    magma_queue_create( cdev, &my_queue );

    //> Define the array sizes
    ldda                  = magma_roundup( Num_Of_Vars, 32 );     //> multiple of 32 by default
    lddb                  = ldda;
    ldd_params            = magma_roundup( Num_Of_Params, 32 );   //> multiple of 32 by default
    dHdx_Index_Size       = Num_Of_Vars*Num_Of_Vars*dHdx_Max_Terms*dHdx_Max_Parts;
    dHdt_Index_Size       = Num_Of_Vars*dHdt_Max_Terms*dHdt_Max_Parts;
    dHdx_PHC_Coeffs_Size  = (Use_P2C) ? ((Num_Of_Coeffs_From_Params+1)*(Max_Order_Of_T+1)) : 0;
    dHdt_PHC_Coeffs_Size  = (Use_P2C) ? ((Num_Of_Coeffs_From_Params+1)*(Max_Order_Of_T)) : 0;
    ldd_phc_Params_Hx     = magma_roundup( dHdx_PHC_Coeffs_Size, 32 );  // multiple of 32 by default
    ldd_phc_Params_Ht     = magma_roundup( dHdt_PHC_Coeffs_Size, 32 );  // multiple of 32 by default

    printf("dHdx_Index_Size      = %5.2f KB\n", (double)(dHdx_Index_Size     *sizeof(magma_int_t))       / 1024.);
    printf("dHdt_Index_Size      = %5.2f KB\n", (double)(dHdt_Index_Size     *sizeof(magma_int_t))       / 1024.);
    printf("dHdx_PHC_Coeffs_Size = %5.2f KB\n", (double)(dHdx_PHC_Coeffs_Size*sizeof(magmaFloatComplex)) / 1024.);
    printf("dHdt_PHC_Coeffs_Size = %5.2f KB\n", (double)(dHdt_PHC_Coeffs_Size*sizeof(magmaFloatComplex)) / 1024.);

    magmaFloatComplex **d_Start_Sols_array      = NULL;
    magmaFloatComplex **d_Homotopy_Sols_array   = NULL;

    //> Define problem file path for problem data reader and output file path for results evaluations
    Problem_File_Path = std::string("../../problems/") + HC_problem;
    Write_Files_Path = std::string("../../") + WRITE_FILES_FOLDER;
}

void GPU_HC_Solver::Allocate_Arrays() {
    //> CPU Allocations
    magma_cmalloc_cpu( &h_Start_Sols,           Num_Of_Tracks*(Num_Of_Vars+1) );
    magma_cmalloc_cpu( &h_Homotopy_Sols,        Num_Of_Tracks*(Num_Of_Vars+1) );
    magma_cmalloc_cpu( &h_Start_Params,         Num_Of_Params );
    magma_cmalloc_cpu( &h_Target_Params,        Num_Of_Params );
    magma_cmalloc_cpu( &h_dHdx_PHC_Coeffs,      dHdx_PHC_Coeffs_Size );
    magma_cmalloc_cpu( &h_dHdt_PHC_Coeffs,      dHdt_PHC_Coeffs_Size );
    magma_imalloc_cpu( &h_dHdx_Index,           dHdx_Index_Size );
    magma_imalloc_cpu( &h_dHdt_Index,           dHdt_Index_Size );
    magma_cmalloc_cpu( &h_diffParams,           Num_Of_Params+1 );
    magma_cmalloc_cpu( &h_GPU_HC_Track_Sols,    (Num_Of_Vars+1)*Num_Of_Tracks );
    magma_cmalloc_cpu( &h_Debug_Purpose,        Num_Of_Tracks );

    h_is_GPU_HC_Sol_Converge = new bool[ Num_Of_Tracks ];
    h_is_GPU_HC_Sol_Infinity = new bool[ Num_Of_Tracks ];

    //> GPU Allocations
    magma_cmalloc( &d_Start_Sols,               Num_Of_Tracks*(Num_Of_Vars+1) );
    magma_cmalloc( &d_Homotopy_Sols,            Num_Of_Tracks*(Num_Of_Vars+1) );
    magma_cmalloc( &d_Start_Params,             ldd_params );
    magma_cmalloc( &d_Target_Params,            ldd_params );
    magma_cmalloc( &d_dHdx_PHC_Coeffs,          ldd_phc_Params_Hx );
    magma_cmalloc( &d_dHdt_PHC_Coeffs,          ldd_phc_Params_Ht );
    magma_imalloc( &d_dHdx_Index,               dHdx_Index_Size );
    magma_imalloc( &d_dHdt_Index,               dHdt_Index_Size );
    magma_cmalloc( &d_diffParams,               ldd_params );
    magma_cmalloc( &d_Debug_Purpose,            Num_Of_Tracks );    

    magma_malloc( (void**) &d_Start_Sols_array,     (Num_Of_Tracks) * sizeof(magmaFloatComplex*) );
    magma_malloc( (void**) &d_Homotopy_Sols_array,  (Num_Of_Tracks) * sizeof(magmaFloatComplex*) );

    cudacheck( cudaMalloc( &d_is_GPU_HC_Sol_Converge, Num_Of_Tracks * sizeof(bool) ));
    cudacheck( cudaMalloc( &d_is_GPU_HC_Sol_Infinity, Num_Of_Tracks * sizeof(bool) ));
}

bool GPU_HC_Solver::Read_Problem_Data() {

    //> Load problem data to arrays
    bool is_Data_Read_Successfully = false;
    Data_Reader Load_Problem_Data(Problem_File_Path, Num_Of_Tracks, Num_Of_Vars);

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
    if (HC_problem == "5pt_rel_pos_geo_form_quat")      magmaHCWrapper::p2c_5pt_rel_pos_geo_form_quat(h_Target_Params, h_Start_Params, h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs);
    else if (HC_problem == "5pt_rel_pos_alg_form_quat") magmaHCWrapper::p2c_5pt_rel_pos_alg_form_quat(h_Target_Params, h_Start_Params, h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs);
    else if (HC_problem == "trifocal_2op1p_30x30_P2C")  magmaHCWrapper::p2c_trifocal_2op1p_30x30(h_Target_Params, h_Start_Params, h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs);
}

void GPU_HC_Solver::Data_Transfer_From_Host_To_Device() {

    //> Compute difference of start and target parameters
    //  (This is used only in the trifocal relative pose problem, i.e., TRIFOCAL_2OP1P_30X30)
    for(int i = 0; i < Num_Of_Params; i++) (h_diffParams)[i] = (h_Target_Params)[i] - (h_Start_Params)[i];
    h_diffParams[Num_Of_Params] = MAGMA_C_ZERO;

    transfer_h2d_time = magma_sync_wtime( my_queue );
    magma_csetmatrix( Num_Of_Vars+1,        Num_Of_Tracks, h_Start_Sols,        (Num_Of_Vars+1),      d_Start_Sols,       (Num_Of_Vars+1),   my_queue );
    magma_csetmatrix( Num_Of_Vars+1,        Num_Of_Tracks, h_Homotopy_Sols,     (Num_Of_Vars+1),      d_Homotopy_Sols,    (Num_Of_Vars+1),   my_queue );
    magma_isetmatrix( dHdx_Index_Size,      (1),           h_dHdx_Index,        dHdx_Index_Size,      d_dHdx_Index,       dHdx_Index_Size,   my_queue );
    magma_isetmatrix( dHdt_Index_Size,      (1),           h_dHdt_Index,        dHdt_Index_Size,      d_dHdt_Index,       dHdt_Index_Size,   my_queue );
    magma_csetmatrix( Num_Of_Params+1,      (1),           h_diffParams,        Num_Of_Params+1,      d_diffParams,       ldd_params,        my_queue );
    magma_csetmatrix( Num_Of_Params,        (1),           h_Start_Params,      Num_Of_Params,        d_Start_Params,     ldd_params,        my_queue );
    magma_csetmatrix( Num_Of_Params,        (1),           h_Target_Params,     Num_Of_Params,        d_Target_Params,    ldd_params,        my_queue );
    if (Use_P2C) {
        magma_csetmatrix( dHdx_PHC_Coeffs_Size, (1),       h_dHdx_PHC_Coeffs,   dHdx_PHC_Coeffs_Size, d_dHdx_PHC_Coeffs,  ldd_phc_Params_Hx, my_queue );
        magma_csetmatrix( dHdt_PHC_Coeffs_Size, (1),       h_dHdt_PHC_Coeffs,   dHdt_PHC_Coeffs_Size, d_dHdt_PHC_Coeffs,  ldd_phc_Params_Ht, my_queue );
    }

    //> connect pointer to 2d arrays
    magma_cset_pointer( d_Start_Sols_array, d_Start_Sols, (Num_Of_Vars+1), 0, 0, (Num_Of_Vars+1), Num_Of_Tracks, my_queue );
    magma_cset_pointer( d_Homotopy_Sols_array,     d_Homotopy_Sols,     (Num_Of_Vars+1), 0, 0, (Num_Of_Vars+1), Num_Of_Tracks, my_queue );
    transfer_h2d_time = magma_sync_wtime( my_queue ) - transfer_h2d_time;
}

void GPU_HC_Solver::Solve_by_GPU_HC() {
    std::cout << "GPU computing ..." << std::endl << std::endl;

    if (HC_problem == "5pt_rel_pos_geo_form_quat") {
        gpu_time = kernel_HC_Solver_5pt_rel_pos_geo_form_quat
                   (my_queue, GPUHC_Max_Steps, GPUHC_Max_Correction_Steps, GPUHC_delta_t_incremental_steps, \
                    d_Start_Sols_array, d_Homotopy_Sols_array, \
                    d_dHdx_Index, d_dHdt_Index, d_dHdx_PHC_Coeffs, d_dHdt_PHC_Coeffs, \
                    d_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Infinity, d_Debug_Purpose);
    }
    else if (HC_problem == "5pt_rel_pos_alg_form_quat") {
        gpu_time = kernel_HC_Solver_5pt_rel_pos_alg_form_quat
                   (my_queue, GPUHC_Max_Steps, GPUHC_Max_Correction_Steps, GPUHC_delta_t_incremental_steps, \
                    d_Start_Sols_array, d_Homotopy_Sols_array, \
                    d_dHdx_Index, d_dHdt_Index, d_dHdx_PHC_Coeffs, d_dHdt_PHC_Coeffs, \
                    d_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Infinity, d_Debug_Purpose);
    }
    else if (HC_problem == "trifocal_2op1p_30x30_P2C") {
        gpu_time = kernel_HC_Solver_trifocal_2op1p_30x30_P2C \
                   (my_queue, GPUHC_Max_Steps, GPUHC_Max_Correction_Steps, GPUHC_delta_t_incremental_steps, \
                    d_Start_Sols_array, d_Homotopy_Sols_array, \
                    d_dHdx_Index, d_dHdt_Index, d_dHdx_PHC_Coeffs, d_dHdt_PHC_Coeffs, \
                    d_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Infinity, d_Debug_Purpose);
    }
    else if (HC_problem == "trifocal_2op1p_30x30") {
        gpu_time = kernel_HC_Solver_trifocal_2op1p_30x30 \
                   (my_queue, GPUHC_Max_Steps, GPUHC_Max_Correction_Steps, GPUHC_delta_t_incremental_steps, \
                    d_Start_Sols_array,  d_Homotopy_Sols_array, \
                    d_Start_Params, d_Target_Params, d_diffParams, \
                    d_dHdx_Index, d_dHdt_Index, \
                    d_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Infinity, d_Debug_Purpose );
    }

    //> Check returns from the GPU kernel
    transfer_d2h_time = magma_sync_wtime( my_queue );
    magma_cgetmatrix( (Num_Of_Vars+1), Num_Of_Tracks, d_Homotopy_Sols,  (Num_Of_Vars+1), h_GPU_HC_Track_Sols,    (Num_Of_Vars+1), my_queue );
    cudacheck( cudaMemcpy( h_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Converge, Num_Of_Tracks*sizeof(bool), cudaMemcpyDeviceToHost) );
    cudacheck( cudaMemcpy( h_is_GPU_HC_Sol_Infinity, d_is_GPU_HC_Sol_Infinity, Num_Of_Tracks*sizeof(bool), cudaMemcpyDeviceToHost) );
    transfer_d2h_time = magma_sync_wtime( my_queue ) - transfer_d2h_time;

#if GPU_DEBUG
    magma_cgetmatrix( Num_Of_Tracks, (1), d_Debug_Purpose, Num_Of_Tracks, h_Debug_Purpose, Num_Of_Vars, my_queue );
#endif

    //> Object for the Evaluations class
    Evaluations Evaluate_GPUHC_Sols( Write_Files_Path, Num_Of_Tracks, Num_Of_Vars );
    Evaluate_GPUHC_Sols.Write_Converged_Sols( h_GPU_HC_Track_Sols, h_is_GPU_HC_Sol_Converge );
    Evaluate_GPUHC_Sols.Evaluate_Sols( h_GPU_HC_Track_Sols, h_is_GPU_HC_Sol_Converge, h_is_GPU_HC_Sol_Infinity );
    Evaluate_GPUHC_Sols.Find_Unique_Sols( h_GPU_HC_Track_Sols, h_is_GPU_HC_Sol_Converge );

    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << "## Solving " << HC_print_problem_name << std::endl;

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
