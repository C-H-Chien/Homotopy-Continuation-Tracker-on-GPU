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
#include <random>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <chrono>
#include <vector>
// #include <type_traits>

#include <cuda.h>
#include <cuda_runtime.h>

//> MAGMA
#include "magma_v2.h"
#include "magma_lapack.h"
#include "magma_internal.h"
#undef max
#undef min

#include "GPU_HC_Solver.hpp"
#include "definitions.hpp"
#include "gpu-kernels/magmaHC-kernels.hpp"

//> Constructor
template< typename T_index_mat >
GPU_HC_Solver<T_index_mat>::GPU_HC_Solver(YAML::Node Problem_Settings_File, int Data_Size_for_Indices_)
    : Problem_Setting_YAML_File(Problem_Settings_File), Data_Size_for_Indices(Data_Size_for_Indices_) 
{
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
    //> (4) GPU Kernel Settings
    Use_Runge_Kutta_in_a_Loop       = Problem_Setting_YAML_File["Use_Runge_Kutta_in_a_Loop"].as<bool>();
    Mem_for_Indices                 = Problem_Setting_YAML_File["Mem_for_Indices"].as<std::string>();
    Inline_Eval_Functions           = Problem_Setting_YAML_File["Inline_Eval_Functions"].as<bool>();
    Limit_Loop_Unroll               = Problem_Setting_YAML_File["Limit_Loop_Unroll"].as<bool>();
    //> (5) Algorithmic Settings
    Truncate_HC_Path_by_Positive_Depths = Problem_Setting_YAML_File["Truncate_HC_Path_by_Positive_Depths"].as<bool>();
    //> (6) RANSAC data
    RANSAC_Dataset_Name             = Problem_Setting_YAML_File["RANSAC_Dataset"].as<std::string>();
    
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
    // ldd_params            = magma_roundup( Num_Of_Params, 32 );   //> multiple of 32 by default
    dHdx_Index_Size       = Num_Of_Vars*Num_Of_Vars*dHdx_Max_Terms*dHdx_Max_Parts;
    dHdt_Index_Size       = Num_Of_Vars*dHdt_Max_Terms*dHdt_Max_Parts;
    dHdx_PHC_Coeffs_Size  = (Use_P2C) ? ((Num_Of_Coeffs_From_Params+1)*(Max_Order_Of_T+1)) : 0;
    dHdt_PHC_Coeffs_Size  = (Use_P2C) ? ((Num_Of_Coeffs_From_Params+1)*(Max_Order_Of_T)) : 0;
    ldd_phc_Params_Hx     = magma_roundup( dHdx_PHC_Coeffs_Size, 32 );  // multiple of 32 by default
    ldd_phc_Params_Ht     = magma_roundup( dHdt_PHC_Coeffs_Size, 32 );  // multiple of 32 by default
#if SHOW_EVAL_INDX_DATA_SIZE
    printf("dHdx_Index_Size      = %5.2f KB\n", (double)(dHdx_Index_Size     *sizeof(T_index_mat))       / 1024.);
    printf("dHdt_Index_Size      = %5.2f KB\n", (double)(dHdt_Index_Size     *sizeof(T_index_mat))       / 1024.);
    if (GPUHC_type == std::string("P2C")) {
        printf("dHdx_PHC_Coeffs_Size = %5.2f KB\n", (double)(dHdx_PHC_Coeffs_Size*sizeof(magmaFloatComplex)) / 1024.);
        printf("dHdt_PHC_Coeffs_Size = %5.2f KB\n", (double)(dHdt_PHC_Coeffs_Size*sizeof(magmaFloatComplex)) / 1024.);
    }
#endif
    magmaFloatComplex **d_Start_Sols_array      = NULL;
    magmaFloatComplex **d_Homotopy_Sols_array   = NULL;

    //> Define problem file path for problem data reader and output file path for results evaluations
    Problem_File_Path = std::string("../../problems/") + HC_problem;
    RANSAC_Data_File_Path = std::string("../../RANSAC_Data/") + HC_problem + "/" + RANSAC_Dataset_Name;
    Write_Files_Path = std::string("../../") + WRITE_FILES_FOLDER;

    //> Evaluations
    Evaluate_GPUHC_Sols = std::shared_ptr<Evaluations>(new Evaluations(Write_Files_Path, Num_Of_Tracks, Num_Of_Vars));
}

template< typename T_index_mat >
void GPU_HC_Solver<T_index_mat>::Allocate_Arrays() {

    //> CPU Allocations
    magma_cmalloc_cpu( &h_Start_Sols,           Num_Of_Tracks*(Num_Of_Vars+1) );
    magma_cmalloc_cpu( &h_Homotopy_Sols,        Num_Of_Tracks*(Num_Of_Vars+1)*NUM_OF_RANSAC_ITERATIONS );
    magma_cmalloc_cpu( &h_Start_Params,         Num_Of_Params+1 );
    magma_cmalloc_cpu( &h_Target_Params,        (Num_Of_Params+1)*NUM_OF_RANSAC_ITERATIONS );
    magma_cmalloc_cpu( &h_dHdx_PHC_Coeffs,      dHdx_PHC_Coeffs_Size );
    magma_cmalloc_cpu( &h_dHdt_PHC_Coeffs,      dHdt_PHC_Coeffs_Size );
    magma_cmalloc_cpu( &h_diffParams,           (Num_Of_Params+1)*NUM_OF_RANSAC_ITERATIONS );
    magma_cmalloc_cpu( &h_GPU_HC_Track_Sols,    Num_Of_Tracks*(Num_Of_Vars+1)*NUM_OF_RANSAC_ITERATIONS ); //> Use to store GPU results from the CPU side
    magma_cmalloc_cpu( &h_Debug_Purpose,        Num_Of_Tracks*NUM_OF_RANSAC_ITERATIONS );

    h_dHdx_Index              = new T_index_mat[ dHdx_Index_Size ];
    h_dHdt_Index              = new T_index_mat[ dHdt_Index_Size ];
    h_is_GPU_HC_Sol_Converge  = new bool[ Num_Of_Tracks*NUM_OF_RANSAC_ITERATIONS ];
    h_is_GPU_HC_Sol_Infinity  = new bool[ Num_Of_Tracks*NUM_OF_RANSAC_ITERATIONS ];

    //> GPU Allocations
    magma_cmalloc( &d_Start_Sols,               Num_Of_Tracks*(Num_Of_Vars+1) );
    magma_cmalloc( &d_Homotopy_Sols,            Num_Of_Tracks*(Num_Of_Vars+1)*NUM_OF_RANSAC_ITERATIONS );
    magma_cmalloc( &d_Start_Params,             Num_Of_Params+1 );
    magma_cmalloc( &d_Target_Params,            (Num_Of_Params+1)*NUM_OF_RANSAC_ITERATIONS );
    magma_cmalloc( &d_dHdx_PHC_Coeffs,          ldd_phc_Params_Hx );
    magma_cmalloc( &d_dHdt_PHC_Coeffs,          ldd_phc_Params_Ht );
    magma_cmalloc( &d_diffParams,               (Num_Of_Params+1)*NUM_OF_RANSAC_ITERATIONS );
    magma_cmalloc( &d_Debug_Purpose,            Num_Of_Tracks*NUM_OF_RANSAC_ITERATIONS );

    cudacheck( cudaMalloc( &d_dHdx_Index,       dHdx_Index_Size *sizeof(T_index_mat)) );
    cudacheck( cudaMalloc( &d_dHdt_Index,       dHdt_Index_Size *sizeof(T_index_mat)) );

    magma_malloc( (void**) &d_Start_Sols_array,     (Num_Of_Tracks)                              * sizeof(magmaFloatComplex*) );
    magma_malloc( (void**) &d_Homotopy_Sols_array,  (Num_Of_Tracks+1)*(NUM_OF_RANSAC_ITERATIONS) * sizeof(magmaFloatComplex*) );

    cudacheck( cudaMalloc( &d_is_GPU_HC_Sol_Converge, (Num_Of_Tracks*NUM_OF_RANSAC_ITERATIONS) * sizeof(bool) ));
    cudacheck( cudaMalloc( &d_is_GPU_HC_Sol_Infinity, (Num_Of_Tracks*NUM_OF_RANSAC_ITERATIONS) * sizeof(bool) ));
}

template< typename T_index_mat >
bool GPU_HC_Solver<T_index_mat>::Read_Problem_Data() {

    //> Load problem data to arrays
    bool is_Data_Read_Successfully = false;

    //> Data reader
    Load_Problem_Data = std::shared_ptr<Data_Reader>(new Data_Reader(Problem_File_Path, RANSAC_Data_File_Path, Num_Of_Tracks, Num_Of_Vars, Num_Of_Params, Use_P2C));

    //> (1) Start parameters
    is_Data_Read_Successfully = Load_Problem_Data->Read_Start_Params( h_Start_Params );
    if (!is_Data_Read_Successfully) { LOG_DATA_LOAD_ERROR("Start Parameters"); return false; }

    // //> (2) Target parameters
    // is_Data_Read_Successfully = Load_Problem_Data->Read_Target_Params( h_Target_Params );
    // if (!is_Data_Read_Successfully) { LOG_DATA_LOAD_ERROR("Target Parameters"); return false; }

    //> (3) Start solutions
    is_Data_Read_Successfully = Load_Problem_Data->Read_Start_Sols( h_Start_Sols, h_Homotopy_Sols );
    if (!is_Data_Read_Successfully) { LOG_DATA_LOAD_ERROR("Start Solutions"); return false; }

    //> (4) dH/dx evaluation indices
    is_Data_Read_Successfully = Load_Problem_Data->Read_dHdx_Indices<T_index_mat>( h_dHdx_Index );
    if (!is_Data_Read_Successfully) { LOG_DATA_LOAD_ERROR("dH/dx Evaluation Indices"); return false; }

    //> (5) dH/dt evaluation indices
    is_Data_Read_Successfully = Load_Problem_Data->Read_dHdt_Indices<T_index_mat>( h_dHdt_Index );
    if (!is_Data_Read_Successfully) { LOG_DATA_LOAD_ERROR("dH/dt Evaluation Indices"); return false; }

    //> (6) Coefficients from parameters, if required
    if (Use_P2C) {
        if (!Load_Problem_Data->Construct_Coeffs_From_Params( HC_problem, h_Target_Params, h_Start_Params, h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs ))
            return false;
    }

    return true;
}

template< typename T_index_mat >
bool GPU_HC_Solver<T_index_mat>::Read_RANSAC_Data( int tp_index ) {
    //> Load problem data to arrays
    bool is_Data_Read_Successfully = false;

    //> (0) get number of triplet edgels 
    Num_Of_Triplet_Edgels = Load_Problem_Data->get_Num_Of_Triplet_Edgels( tp_index );
    if (Num_Of_Triplet_Edgels == 0)
        return false;
    // LOG_INFO_MESG("Number of triplet edgels = " + std::to_string(Num_Of_Triplet_Edgels));

    //> Allocate triplet edgel locations and tangents arrays
    h_Triplet_Edge_Locations  = new float[ Num_Of_Triplet_Edgels*2*3 ];
    h_Triplet_Edge_Tangents   = new float[ Num_Of_Triplet_Edgels*2*3 ];

    //> (1) Camera intrinsic/extrinsic matrices
    is_Data_Read_Successfully = Load_Problem_Data->Read_Camera_Matrices( h_Camera_Pose21, h_Camera_Pose31, h_Camera_Intrinsic_Matrix, tp_index );
    if (!is_Data_Read_Successfully) { LOG_DATA_LOAD_ERROR("Camera Matrices"); return false; }

    //> (2) Triplet edgel correspondences, convert data to edgel locations and tangents arrays
    Load_Problem_Data->Read_Triplet_Edgels( h_Triplet_Edge_Locations, h_Triplet_Edge_Tangents );

    return true;
}

template< typename T_index_mat >
void GPU_HC_Solver<T_index_mat>::Prepare_Target_Params( unsigned rand_seed_ ) {
    unsigned Edgel_Indx[3] = {0};
    const int offset_for_tangents = 18;
    std::array<int, 3> picked_samples;
    std::vector< std::array<int, 3> > target_params_match_indices;
#if FEED_RANDOM_SEED
    srand (time(NULL));
#else
    srand (rand_seed_);
#endif    
    for (int ti = 0; ti < NUM_OF_RANSAC_ITERATIONS; ti++) {
        //> (1) Converting from triplet edgels to target parameters
        //> Pick 3 triplet edgels
        while(1) {
            for (int ri = 0; ri < 3; ri++) Edgel_Indx[ri] = rand() % Num_Of_Triplet_Edgels;
            if ( (Edgel_Indx[0] != Edgel_Indx[1]) && (Edgel_Indx[0] != Edgel_Indx[1]) && (Edgel_Indx[1] != Edgel_Indx[2]) ) break;
        }
        for (int i = 0; i < 3; i++) picked_samples[i] = Edgel_Indx[i];
        target_params_match_indices.push_back(picked_samples);

        //> Locations of the edgels
        for (int i = 0; i < 3; i++) {
            int ei = Edgel_Indx[i];
            for (int j = 0; j < 6; j++) {
                (h_Target_Params + ti*(Num_Of_Params+1))[i*6 + j] = MAGMA_C_MAKE(h_Triplet_Edge_Locations(ei, j), 0.0);
            }
        }
        //> Tangents of the edgels
        for (int i = 0; i < 2; i++) {
            int ei = Edgel_Indx[i];
            for (int j = 0; j < 6; j++) {
                (h_Target_Params + ti*(Num_Of_Params+1))[i*6 + j + offset_for_tangents] = MAGMA_C_MAKE(h_Triplet_Edge_Tangents(ei, j), 0.0);
            }
        }
        (h_Target_Params + ti*(Num_Of_Params+1))[30] = MAGMA_C_MAKE(1.0, 0.0);
        (h_Target_Params + ti*(Num_Of_Params+1))[31] = MAGMA_C_MAKE(0.5, 0.0);
        (h_Target_Params + ti*(Num_Of_Params+1))[32] = MAGMA_C_MAKE(-1.0, 0.0);
        (h_Target_Params + ti*(Num_Of_Params+1))[33] = MAGMA_C_ONE;

        //> (2) Compute the difference of the start and target parameters in h_diffParams
        for (int i = 0; i <= Num_Of_Params; i++) 
            (h_diffParams + ti*(Num_Of_Params+1))[i] = (h_Target_Params + ti*(Num_Of_Params+1))[i] - (h_Start_Params)[i];
    }
#if RANSAC_DEBUG
    int sample_index = 0;
    if (NUM_OF_RANSAC_ITERATIONS > sample_index) Load_Problem_Data->Print_Out_Target_Params_from_Triplet_Edgels(sample_index, target_params_match_indices, h_Target_Params);
    else Load_Problem_Data->Print_Out_Target_Params_from_Triplet_Edgels(0, target_params_match_indices, h_Target_Params);
#endif
}

template< typename T_index_mat >
void GPU_HC_Solver<T_index_mat>::Data_Transfer_From_Host_To_Device() {

    transfer_h2d_time = magma_sync_wtime( my_queue );
    magma_csetmatrix( Num_Of_Vars+1,   Num_Of_Tracks*NUM_OF_RANSAC_ITERATIONS,  h_Homotopy_Sols,  (Num_Of_Vars+1),  d_Homotopy_Sols, Num_Of_Vars+1,     my_queue );
    magma_csetmatrix( Num_Of_Vars+1,   Num_Of_Tracks,                           h_Start_Sols,     (Num_Of_Vars+1),  d_Start_Sols,    Num_Of_Vars+1,     my_queue );
    magma_csetmatrix( Num_Of_Params+1, (1),                                     h_Start_Params,   Num_Of_Params+1,  d_Start_Params,  Num_Of_Params+1,   my_queue );
    magma_csetmatrix( (Num_Of_Params+1)*NUM_OF_RANSAC_ITERATIONS, (1), h_diffParams,    (Num_Of_Params+1)*NUM_OF_RANSAC_ITERATIONS, d_diffParams,    (Num_Of_Params+1)*NUM_OF_RANSAC_ITERATIONS, my_queue );
    magma_csetmatrix( (Num_Of_Params+1)*NUM_OF_RANSAC_ITERATIONS, (1), h_Target_Params, (Num_Of_Params+1)*NUM_OF_RANSAC_ITERATIONS, d_Target_Params, (Num_Of_Params+1)*NUM_OF_RANSAC_ITERATIONS, my_queue );
    cudacheck( cudaMemcpy( d_dHdx_Index, h_dHdx_Index,     dHdx_Index_Size * sizeof(T_index_mat), cudaMemcpyHostToDevice) );
    cudacheck( cudaMemcpy( d_dHdt_Index, h_dHdt_Index,     dHdt_Index_Size * sizeof(T_index_mat), cudaMemcpyHostToDevice) );

    if (Use_P2C) {
        magma_csetmatrix( dHdx_PHC_Coeffs_Size, (1),       h_dHdx_PHC_Coeffs,   dHdx_PHC_Coeffs_Size, d_dHdx_PHC_Coeffs,  ldd_phc_Params_Hx, my_queue );
        magma_csetmatrix( dHdt_PHC_Coeffs_Size, (1),       h_dHdt_PHC_Coeffs,   dHdt_PHC_Coeffs_Size, d_dHdt_PHC_Coeffs,  ldd_phc_Params_Ht, my_queue );
    }

    //> connect pointer to 2d arrays
    magma_cset_pointer( d_Start_Sols_array,    d_Start_Sols,     (Num_Of_Vars+1), 0, 0, (Num_Of_Vars+1), Num_Of_Tracks, my_queue );
    magma_cset_pointer( d_Homotopy_Sols_array, d_Homotopy_Sols,  (Num_Of_Vars+1), 0, 0, (Num_Of_Vars+1), Num_Of_Tracks*NUM_OF_RANSAC_ITERATIONS, my_queue );
    transfer_h2d_time = magma_sync_wtime( my_queue ) - transfer_h2d_time;
}

template< >
void GPU_HC_Solver<int>::Solve_by_GPU_HC() {
    std::cout << "GPU computing ..." << std::endl << std::endl;

    if (HC_problem == "trifocal_2op1p_30x30") {

        if ( Use_P2C ) {
            //> (#1) Naive approach
            LOG_INFO_MESG("kernel_GPUHC_trifocal_2op1p_30x30, Naive GPU-HC approach");
            gpu_time = kernel_GPUHC_trifocal_2op1p_30x30_P2C \
                       (my_queue, GPUHC_Max_Steps, GPUHC_Max_Correction_Steps, GPUHC_delta_t_incremental_steps, \
                        d_Start_Sols_array, d_Homotopy_Sols_array, \
                        d_dHdx_Index, d_dHdt_Index, d_dHdx_PHC_Coeffs, d_dHdt_PHC_Coeffs, \
                        d_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Infinity, d_Debug_Purpose);
        }
        else if ( Mem_for_Indices == "GM" && !Use_Runge_Kutta_in_a_Loop && !Inline_Eval_Functions && !Limit_Loop_Unroll && !Truncate_HC_Path_by_Positive_Depths ) {
            //> (#2) 32-bit indices in global memory + no Runge-Kutta in a loop + no inlined device functions + no limited loop unroll + no truncated HC paths
            LOG_INFO_MESG("kernel_GPUHC_trifocal_2op1p_30x30_GM32b");
            gpu_time = kernel_GPUHC_trifocal_2op1p_30x30_GM32b \
                       (my_queue, GPUHC_Max_Steps, GPUHC_Max_Correction_Steps, GPUHC_delta_t_incremental_steps, d_Start_Sols_array,  d_Homotopy_Sols_array, \
                        d_Start_Params, d_Target_Params, d_diffParams, d_dHdx_Index, d_dHdt_Index, \
                        d_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Infinity, d_Debug_Purpose );
        }
        else if ( Mem_for_Indices == "GM" && Use_Runge_Kutta_in_a_Loop && !Inline_Eval_Functions && !Limit_Loop_Unroll && !Truncate_HC_Path_by_Positive_Depths) {
            //> (#3) 32-bit indices in global memory + Runge-Kutta in a loop + inlined device functions + no limited loop unroll + no truncated HC paths
            LOG_INFO_MESG("kernel_GPUHC_trifocal_2op1p_30x30_GM32b_RKL");
            gpu_time = kernel_GPUHC_trifocal_2op1p_30x30_GM32b_RKL \
                       (my_queue, GPUHC_Max_Steps, GPUHC_Max_Correction_Steps, GPUHC_delta_t_incremental_steps, d_Start_Sols_array,  d_Homotopy_Sols_array, \
                        d_Start_Params, d_Target_Params, d_diffParams, d_dHdx_Index, d_dHdt_Index, \
                        d_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Infinity, d_Debug_Purpose );
        }
        else if ( Mem_for_Indices == "GM" && Use_Runge_Kutta_in_a_Loop && Inline_Eval_Functions && !Limit_Loop_Unroll && !Truncate_HC_Path_by_Positive_Depths) {
            //> (#4) 32-bit indices in global memory + Runge-Kutta in a loop + inlined device functions + no limited loop unroll + no truncated HC paths
            LOG_INFO_MESG("kernel_GPUHC_trifocal_2op1p_30x30_GM32b_RKL_inline");
            gpu_time = kernel_GPUHC_trifocal_2op1p_30x30_GM32b_RKL_inline \
                       (my_queue, GPUHC_Max_Steps, GPUHC_Max_Correction_Steps, GPUHC_delta_t_incremental_steps, d_Start_Sols_array,  d_Homotopy_Sols_array, \
                        d_Start_Params, d_Target_Params, d_diffParams, d_dHdx_Index, d_dHdt_Index, \
                        d_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Infinity, d_Debug_Purpose );
        }
        else if ( Mem_for_Indices == "GM" && Use_Runge_Kutta_in_a_Loop && !Inline_Eval_Functions && Limit_Loop_Unroll && !Truncate_HC_Path_by_Positive_Depths) {
            //> (#5) 32-bit indices in global memory + Runge-Kutta in a loop + no inlined device functions + limited loop unroll + no truncated HC paths
            LOG_INFO_MESG("kernel_GPUHC_trifocal_2op1p_30x30_GM32b_RKL_LimUnroll");
            gpu_time = kernel_GPUHC_trifocal_2op1p_30x30_GM32b_RKL_LimUnroll \
                       (my_queue, GPUHC_Max_Steps, GPUHC_Max_Correction_Steps, GPUHC_delta_t_incremental_steps, d_Start_Sols_array,  d_Homotopy_Sols_array, \
                        d_Start_Params, d_Target_Params, d_diffParams, d_dHdx_Index, d_dHdt_Index, \
                        d_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Infinity, d_Debug_Purpose );
        }
        else if ( Mem_for_Indices == "GM" && Use_Runge_Kutta_in_a_Loop && Inline_Eval_Functions && Limit_Loop_Unroll && !Truncate_HC_Path_by_Positive_Depths) {
            //> (#6) 32-bit indices in global memory + Runge-Kutta in a loop + inlined device functions + limited loop unroll + no truncated HC paths
            LOG_INFO_MESG("kernel_GPUHC_trifocal_2op1p_30x30_GM32b_RKL_inline_LimUnroll");
            gpu_time = kernel_GPUHC_trifocal_2op1p_30x30_GM32b_RKL_inline_LimUnroll \
                       (my_queue, GPUHC_Max_Steps, GPUHC_Max_Correction_Steps, GPUHC_delta_t_incremental_steps, d_Start_Sols_array,  d_Homotopy_Sols_array, \
                        d_Start_Params, d_Target_Params, d_diffParams, d_dHdx_Index, d_dHdt_Index, \
                        d_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Infinity, d_Debug_Purpose );
        }
        else if ( Mem_for_Indices == "GM" && Use_Runge_Kutta_in_a_Loop && Inline_Eval_Functions && Limit_Loop_Unroll && Truncate_HC_Path_by_Positive_Depths) {
            //> (#7) 32-bit indices in global memory + Runge-Kutta in a loop + inlined device functions + limited loop unroll + truncated HC paths
            LOG_INFO_MESG("kernel_GPUHC_trifocal_2op1p_30x30_GM32b_RKL_inline_LimUnroll_TrunPaths");
            gpu_time = kernel_GPUHC_trifocal_2op1p_30x30_GM32b_RKL_inline_LimUnroll_TrunPaths \
                       (my_queue, GPUHC_Max_Steps, GPUHC_Max_Correction_Steps, GPUHC_delta_t_incremental_steps, d_Start_Sols_array,  d_Homotopy_Sols_array, \
                        d_Start_Params, d_Target_Params, d_diffParams, d_dHdx_Index, d_dHdt_Index, \
                        d_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Infinity, d_Debug_Purpose );
        }
        else {
            LOG_ERROR("Invalid configurations on the GPU settings / Algorithmic settings!");
            return;
        }
    }

    //> Check returns from the GPU kernel
    transfer_d2h_time = magma_sync_wtime( my_queue );
    magma_cgetmatrix( (Num_Of_Vars+1), Num_Of_Tracks*NUM_OF_RANSAC_ITERATIONS, d_Homotopy_Sols,  (Num_Of_Vars+1), h_GPU_HC_Track_Sols,    (Num_Of_Vars+1), my_queue );
    cudacheck( cudaMemcpy( h_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Converge, Num_Of_Tracks*NUM_OF_RANSAC_ITERATIONS*sizeof(bool), cudaMemcpyDeviceToHost) );
    cudacheck( cudaMemcpy( h_is_GPU_HC_Sol_Infinity, d_is_GPU_HC_Sol_Infinity, Num_Of_Tracks*NUM_OF_RANSAC_ITERATIONS*sizeof(bool), cudaMemcpyDeviceToHost) );
    transfer_d2h_time = magma_sync_wtime( my_queue ) - transfer_d2h_time;

#if GPU_DEBUG
    magma_cgetmatrix( Num_Of_Tracks, NUM_OF_RANSAC_ITERATIONS, d_Debug_Purpose, Num_Of_Tracks, h_Debug_Purpose, Num_Of_Tracks, my_queue );
#endif

    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << "## Solving " << HC_print_problem_name << std::endl;

    //> Print out timings
    printf("## Timings:\n");
    printf(" - GPU Computation Time = %7.2f (ms)\n", (gpu_time)*1000);
    printf(" - H2D Transfer Time    = %7.2f (ms)\n", (transfer_h2d_time)*1000);
    printf(" - D2H Transfer Time    = %7.2f (ms)\n", (transfer_d2h_time)*1000);

    //> Object for the Evaluations class
    // Evaluate_GPUHC_Sols->Write_Converged_Sols( h_GPU_HC_Track_Sols, h_is_GPU_HC_Sol_Converge );
    Evaluate_GPUHC_Sols->Flush_Out_Data();
    Evaluate_GPUHC_Sols->Evaluate_RANSAC_GPUHC_Sols( h_GPU_HC_Track_Sols, h_is_GPU_HC_Sol_Converge, h_is_GPU_HC_Sol_Infinity );
    // Evaluate_GPUHC_Sols->Find_Unique_Sols( h_GPU_HC_Track_Sols, h_is_GPU_HC_Sol_Converge );

    //> Print out evaluation results
    std::cout << "## Evaluation of GPU-HC Solutions: "      << std::endl;
    std::cout << " - Number of Converged Solutions:       " << Evaluate_GPUHC_Sols->Num_Of_Coverged_Sols << std::endl;
    std::cout << " - Number of Real Solutions:            " << Evaluate_GPUHC_Sols->Num_Of_Real_Sols << std::endl;
    std::cout << " - Number of Infinity Failed Solutions: " << Evaluate_GPUHC_Sols->Num_Of_Inf_Sols << std::endl;
    // std::cout << " - Number of Unique Solutions:          " << Evaluate_GPUHC_Sols->Num_Of_Unique_Sols << std::endl;

    Evaluate_GPUHC_Sols->Transform_GPUHC_Sols_to_Trifocal_Relative_Pose( h_GPU_HC_Track_Sols, h_is_GPU_HC_Sol_Converge, h_Camera_Intrinsic_Matrix );
    bool find_good_sol = Evaluate_GPUHC_Sols->get_Solution_with_Maximal_Support( Num_Of_Triplet_Edgels, h_Triplet_Edge_Locations, h_Triplet_Edge_Tangents, h_Camera_Intrinsic_Matrix );
    if (find_good_sol)
        Evaluate_GPUHC_Sols->Measure_Relative_Pose_Error( h_Camera_Pose21, h_Camera_Pose31, h_Debug_Purpose );

    // Evaluate_GPUHC_Sols->get_HC_Steps_of_Actual_Sols( h_Debug_Purpose );
    // Evaluate_GPUHC_Sols->Flush_Out_Data();
    // Evaluate_GPUHC_Sols->Transform_GPUHC_Sols_to_Trifocal_Relative_Pose( h_GPU_HC_Track_Sols, h_is_GPU_HC_Sol_Converge, h_Camera_Intrinsic_Matrix );
    // Evaluate_GPUHC_Sols->Measure_Relative_Pose_Error_from_All_Real_Sols( h_Camera_Pose21, h_Camera_Pose31, h_Debug_Purpose );
    
    // for (int i = 0; i < Evaluate_GPUHC_Sols->HC_steps_of_actual_solutions.size(); i++) {
    //     std::cout << Evaluate_GPUHC_Sols->HC_steps_of_actual_solutions[i] << ", ";
    //     GPUHC_Actual_Sols_Steps_Collections.push_back( Evaluate_GPUHC_Sols->HC_steps_of_actual_solutions[i] );
    // }
    // std::cout << std::endl;

    // if (Evaluate_GPUHC_Sols->success_flag) {
    //     std::cout << "## Found solution matched with GT: " << std::endl;
    //     std::cout << " - Residual of R21: " << Evaluate_GPUHC_Sols->Min_Residual_R21 << " (rad)" << std::endl;
    //     std::cout << " - Residual of R31: " << Evaluate_GPUHC_Sols->Min_Residual_R31 << " (rad)" << std::endl;
    //     std::cout << " - Residual of t21: " << Evaluate_GPUHC_Sols->Min_Residual_t21 << " (m)" << std::endl;
    //     std::cout << " - Residual of t31: " << Evaluate_GPUHC_Sols->Min_Residual_t31 << " (m)" << std::endl;
    // }
    // else {
    //     std::cout << "## Not found a solution matched with GT: " << std::endl;
    //     std::cout << " - Residual of R21: " << Evaluate_GPUHC_Sols->Min_Residual_R21 << " (rad)" << std::endl;
    //     std::cout << " - Residual of R31: " << Evaluate_GPUHC_Sols->Min_Residual_R31 << " (rad)" << std::endl;
    //     std::cout << " - Residual of t21: " << Evaluate_GPUHC_Sols->Min_Residual_t21 << " (m)" << std::endl;
    //     std::cout << " - Residual of t31: " << Evaluate_GPUHC_Sols->Min_Residual_t31 << " (m)" << std::endl;
    // }
}

template< >
void GPU_HC_Solver<char>::Solve_by_GPU_HC() {
    std::cout << "GPU computing ..." << std::endl << std::endl;

    if (HC_problem == "trifocal_2op1p_30x30") {

        if ( Mem_for_Indices == "GM" && Use_Runge_Kutta_in_a_Loop && Inline_Eval_Functions && !Limit_Loop_Unroll && !Truncate_HC_Path_by_Positive_Depths) {
            //> (#3) 8-bit indices in global memory + Runge-Kutta in a loop + inlined device functions + no limited loop unroll + no truncated HC paths
            LOG_INFO_MESG("kernel_GPUHC_trifocal_2op1p_30x30_GM8b_RKL_inline");
            gpu_time = kernel_GPUHC_trifocal_2op1p_30x30_GM8b_RKL_inline \
                       (my_queue, GPUHC_Max_Steps, GPUHC_Max_Correction_Steps, GPUHC_delta_t_incremental_steps, d_Start_Sols_array,  d_Homotopy_Sols_array, \
                        d_Start_Params, d_Target_Params, d_diffParams, d_dHdx_Index, d_dHdt_Index, \
                        d_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Infinity, d_Debug_Purpose );
        }
        else if ( Mem_for_Indices == "SM" && Use_Runge_Kutta_in_a_Loop && Inline_Eval_Functions && !Limit_Loop_Unroll && !Truncate_HC_Path_by_Positive_Depths) {
            //> (#4) 8-bit indices in shared memory + Runge-Kutta in a loop + inlined device functions + no limited loop unroll + no truncated HC paths
            LOG_INFO_MESG("kernel_GPUHC_trifocal_2op1p_30x30_SM8b_RKL_inline");
            gpu_time = kernel_GPUHC_trifocal_2op1p_30x30_SM8b_RKL_inline \
                       (my_queue, GPUHC_Max_Steps, GPUHC_Max_Correction_Steps, GPUHC_delta_t_incremental_steps, d_Start_Sols_array,  d_Homotopy_Sols_array, \
                        d_Start_Params, d_Target_Params, d_diffParams, d_dHdx_Index, d_dHdt_Index, \
                        d_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Infinity, d_Debug_Purpose );
        }
        else if ( Mem_for_Indices == "SM" && Use_Runge_Kutta_in_a_Loop && !Inline_Eval_Functions && !Limit_Loop_Unroll && !Truncate_HC_Path_by_Positive_Depths) {
            //> (#5) 8-bit indices in shared memory + Runge-Kutta in a loop + no inlined device functions + no limited loop unroll + no truncated HC paths
            LOG_INFO_MESG("kernel_GPUHC_trifocal_2op1p_30x30_SM8b_RKL");
            gpu_time = kernel_GPUHC_trifocal_2op1p_30x30_SM8b_RKL \
                       (my_queue, GPUHC_Max_Steps, GPUHC_Max_Correction_Steps, GPUHC_delta_t_incremental_steps, d_Start_Sols_array,  d_Homotopy_Sols_array, \
                        d_Start_Params, d_Target_Params, d_diffParams, d_dHdx_Index, d_dHdt_Index, \
                        d_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Infinity, d_Debug_Purpose );
        }
        else if ( Mem_for_Indices == "SM" && Use_Runge_Kutta_in_a_Loop && !Inline_Eval_Functions && Limit_Loop_Unroll && !Truncate_HC_Path_by_Positive_Depths) {
            //> (#6) 8-bit indices in shared memory + Runge-Kutta in a loop + no inlined device functions + limited loop unroll + no truncated HC paths
            LOG_INFO_MESG("kernel_GPUHC_trifocal_2op1p_30x30_SM8b_RKL_LimUnroll");
            gpu_time = kernel_GPUHC_trifocal_2op1p_30x30_SM8b_RKL_LimUnroll \
                       (my_queue, GPUHC_Max_Steps, GPUHC_Max_Correction_Steps, GPUHC_delta_t_incremental_steps, d_Start_Sols_array,  d_Homotopy_Sols_array, \
                        d_Start_Params, d_Target_Params, d_diffParams, d_dHdx_Index, d_dHdt_Index, \
                        d_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Infinity, d_Debug_Purpose );
        }
        else {
            LOG_ERROR("Invalid configurations on the GPU settings / Algorithmic settings!");
            return;
        }
    }

    //> Check returns from the GPU kernel
    transfer_d2h_time = magma_sync_wtime( my_queue );
    magma_cgetmatrix( (Num_Of_Vars+1), Num_Of_Tracks*NUM_OF_RANSAC_ITERATIONS, d_Homotopy_Sols,  (Num_Of_Vars+1), h_GPU_HC_Track_Sols,    (Num_Of_Vars+1), my_queue );
    cudacheck( cudaMemcpy( h_is_GPU_HC_Sol_Converge, d_is_GPU_HC_Sol_Converge, Num_Of_Tracks*NUM_OF_RANSAC_ITERATIONS*sizeof(bool), cudaMemcpyDeviceToHost) );
    cudacheck( cudaMemcpy( h_is_GPU_HC_Sol_Infinity, d_is_GPU_HC_Sol_Infinity, Num_Of_Tracks*NUM_OF_RANSAC_ITERATIONS*sizeof(bool), cudaMemcpyDeviceToHost) );
    transfer_d2h_time = magma_sync_wtime( my_queue ) - transfer_d2h_time;

#if GPU_DEBUG
    magma_cgetmatrix( Num_Of_Tracks, NUM_OF_RANSAC_ITERATIONS, d_Debug_Purpose, Num_Of_Tracks, h_Debug_Purpose, Num_Of_Tracks, my_queue );
#endif

    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << "## Solving " << HC_print_problem_name << std::endl;

    //> Print out timings
    printf("## Timings:\n");
    printf(" - GPU Computation Time = %7.2f (ms)\n", (gpu_time)*1000);
    printf(" - H2D Transfer Time    = %7.2f (ms)\n", (transfer_h2d_time)*1000);
    printf(" - D2H Transfer Time    = %7.2f (ms)\n", (transfer_d2h_time)*1000);

    //> Object for the Evaluations class
    // Evaluate_GPUHC_Sols->Write_Converged_Sols( h_GPU_HC_Track_Sols, h_is_GPU_HC_Sol_Converge );
    Evaluate_GPUHC_Sols->Flush_Out_Data();
    Evaluate_GPUHC_Sols->Evaluate_RANSAC_GPUHC_Sols( h_GPU_HC_Track_Sols, h_is_GPU_HC_Sol_Converge, h_is_GPU_HC_Sol_Infinity );
    // Evaluate_GPUHC_Sols->Find_Unique_Sols( h_GPU_HC_Track_Sols, h_is_GPU_HC_Sol_Converge );

    //> Print out evaluation results
    std::cout << "## Evaluation of GPU-HC Solutions: "      << std::endl;
    std::cout << " - Number of Converged Solutions:       " << Evaluate_GPUHC_Sols->Num_Of_Coverged_Sols << std::endl;
    std::cout << " - Number of Real Solutions:            " << Evaluate_GPUHC_Sols->Num_Of_Real_Sols << std::endl;
    std::cout << " - Number of Infinity Failed Solutions: " << Evaluate_GPUHC_Sols->Num_Of_Inf_Sols << std::endl;
    // std::cout << " - Number of Unique Solutions:          " << Evaluate_GPUHC_Sols->Num_Of_Unique_Sols << std::endl;

    Evaluate_GPUHC_Sols->Transform_GPUHC_Sols_to_Trifocal_Relative_Pose( h_GPU_HC_Track_Sols, h_is_GPU_HC_Sol_Converge, h_Camera_Intrinsic_Matrix );
    bool find_good_sol = Evaluate_GPUHC_Sols->get_Solution_with_Maximal_Support( Num_Of_Triplet_Edgels, h_Triplet_Edge_Locations, h_Triplet_Edge_Tangents, h_Camera_Intrinsic_Matrix );
    if (find_good_sol)
        Evaluate_GPUHC_Sols->Measure_Relative_Pose_Error( h_Camera_Pose21, h_Camera_Pose31, h_Debug_Purpose );

    // Evaluate_GPUHC_Sols->get_HC_Steps_of_Actual_Sols( h_Debug_Purpose );
    // Evaluate_GPUHC_Sols->Flush_Out_Data();
    // Evaluate_GPUHC_Sols->Transform_GPUHC_Sols_to_Trifocal_Relative_Pose( h_GPU_HC_Track_Sols, h_is_GPU_HC_Sol_Converge, h_Camera_Intrinsic_Matrix );
    // Evaluate_GPUHC_Sols->Measure_Relative_Pose_Error_from_All_Real_Sols( h_Camera_Pose21, h_Camera_Pose31, h_Debug_Purpose );
    
    // for (int i = 0; i < Evaluate_GPUHC_Sols->HC_steps_of_actual_solutions.size(); i++) {
    //     std::cout << Evaluate_GPUHC_Sols->HC_steps_of_actual_solutions[i] << ", ";
    //     GPUHC_Actual_Sols_Steps_Collections.push_back( Evaluate_GPUHC_Sols->HC_steps_of_actual_solutions[i] );
    // }
    // std::cout << std::endl;

    // if (Evaluate_GPUHC_Sols->success_flag) {
    //     std::cout << "## Found solution matched with GT: " << std::endl;
    //     std::cout << " - Residual of R21: " << Evaluate_GPUHC_Sols->Min_Residual_R21 << " (rad)" << std::endl;
    //     std::cout << " - Residual of R31: " << Evaluate_GPUHC_Sols->Min_Residual_R31 << " (rad)" << std::endl;
    //     std::cout << " - Residual of t21: " << Evaluate_GPUHC_Sols->Min_Residual_t21 << " (m)" << std::endl;
    //     std::cout << " - Residual of t31: " << Evaluate_GPUHC_Sols->Min_Residual_t31 << " (m)" << std::endl;
    // }
    // else {
    //     std::cout << "## Not found a solution matched with GT: " << std::endl;
    //     std::cout << " - Residual of R21: " << Evaluate_GPUHC_Sols->Min_Residual_R21 << " (rad)" << std::endl;
    //     std::cout << " - Residual of R31: " << Evaluate_GPUHC_Sols->Min_Residual_R31 << " (rad)" << std::endl;
    //     std::cout << " - Residual of t21: " << Evaluate_GPUHC_Sols->Min_Residual_t21 << " (m)" << std::endl;
    //     std::cout << " - Residual of t31: " << Evaluate_GPUHC_Sols->Min_Residual_t31 << " (m)" << std::endl;
    // }
}

// template< typename T_index_mat >
// void GPU_HC_Solver<T_index_mat>::Export_Data() {
//     Evaluate_GPUHC_Sols->Write_HC_Steps_of_Actual_Solutions( GPUHC_Actual_Sols_Steps_Collections );
// }

template< typename T_index_mat >
void GPU_HC_Solver<T_index_mat>::Free_Triplet_Edgels_Mem() {
    delete [] h_Triplet_Edge_Locations;
    delete [] h_Triplet_Edge_Tangents;
}

template< typename T_index_mat >
GPU_HC_Solver<T_index_mat>::~GPU_HC_Solver() {

    magma_queue_destroy( my_queue );

    delete [] h_is_GPU_HC_Sol_Converge;
    delete [] h_is_GPU_HC_Sol_Infinity;
    delete [] h_dHdx_Index;
    delete [] h_dHdt_Index;

    magma_free_cpu( h_Start_Sols );
    magma_free_cpu( h_Homotopy_Sols );
    magma_free_cpu( h_Start_Params );
    magma_free_cpu( h_Target_Params );
    magma_free_cpu( h_dHdx_PHC_Coeffs );
    magma_free_cpu( h_dHdt_PHC_Coeffs );
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
    magma_free( d_dHdx_PHC_Coeffs );
    magma_free( d_dHdt_PHC_Coeffs );
    magma_free( d_Debug_Purpose );

    cudacheck( cudaFree( d_dHdx_Index ) );
    cudacheck( cudaFree( d_dHdt_Index ) );

    fflush( stdout );
    printf( "\n" );
    magma_finalize();
}

template class GPU_HC_Solver<int>;
template class GPU_HC_Solver<char>;

#endif
