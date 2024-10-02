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
GPU_HC_Solver::GPU_HC_Solver(YAML::Node Problem_Settings_File)
    : Problem_Setting_YAML_File(Problem_Settings_File)
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
    Use_Merge_Code_Optimization     = Problem_Setting_YAML_File["Code_Merged_Optimization"].as<bool>();
    //> (5) Algorithmic Settings
    Truncate_HC_Path_by_Positive_Depths = Problem_Setting_YAML_File["Truncate_HC_Path_by_Positive_Depths"].as<bool>();
    Abort_RANSAC_by_Good_Sol            = Problem_Setting_YAML_File["Abort_RANSAC_by_Good_Sol"].as<bool>();
    //> (6) RANSAC data
    RANSAC_Dataset_Name             = Problem_Setting_YAML_File["RANSAC_Dataset"].as<std::string>();
    //> (7) Multiple GPUs
    Num_Of_GPUs                     = Problem_Setting_YAML_File["Num_Of_GPUs"].as<int>();
    
    //> Check if conversion from problem parameters to polynomial coefficients is needed
    if (GPUHC_type == std::string("P2C")) {
        Use_P2C = true;
        Num_Of_Coeffs_From_Params   = Problem_Setting_YAML_File["Num_Of_Coeffs_From_Params"].as<int>();
    }
    else Use_P2C = false;

    //> Initialization
    magma_init();
    magma_print_environment();
    cudaGetDeviceCount( &device_count );
    check_multiGPUs();
    
    // magma_getdevice( &cdev );
    // magma_getdevices( devices, Num_Of_GPUs, &Num_of_Devices_from_magma );

    //> MAGMA queues on multiple GPUs, one queue per GPU
    for (int gpu_id = 0; gpu_id < Num_Of_GPUs; gpu_id++) {
        dev_id = (Num_Of_GPUs == 1) ? (SET_GPU_DEVICE_ID) : (gpu_id);
        magma_setdevice( dev_id );
        magma_queue_create( dev_id, &gpu_queues[gpu_id] );
    }

    //> Calculate the share of each GPU: divide the number of RANSAC iterations by the number of requested GPUs
    for(magma_int_t gpu_id = 0; gpu_id < Num_Of_GPUs; gpu_id++) {
        sub_RANSAC_iters[gpu_id] = NUM_OF_RANSAC_ITERATIONS / Num_Of_GPUs + ((gpu_id < (NUM_OF_RANSAC_ITERATIONS % Num_Of_GPUs)) ? 1: 0);
		printf("GPU %2d computes %2d RANSAC iterations\n", gpu_id, sub_RANSAC_iters[gpu_id]);
	}
    
    //> Define the array sizes
    dHdx_Index_Size                 = Num_Of_Vars*Num_Of_Vars*dHdx_Max_Terms*dHdx_Max_Parts;
    dHdt_Index_Size                 = Num_Of_Vars*dHdt_Max_Terms*dHdt_Max_Parts;
    unified_dHdx_dHdt_Index_Size    = dHdx_Index_Size + dHdt_Index_Size;
    dHdx_PHC_Coeffs_Size            = (Use_P2C) ? ((Num_Of_Coeffs_From_Params+1)*(Max_Order_Of_T+1)) : 0;
    dHdt_PHC_Coeffs_Size            = (Use_P2C) ? ((Num_Of_Coeffs_From_Params+1)*(Max_Order_Of_T)) : 0;
#if SHOW_EVAL_INDX_DATA_SIZE
    printf("dHdx_Index_Size      = %5.2f KB\n", (double)(dHdx_Index_Size     *sizeof(int))       / 1024.);
    printf("dHdt_Index_Size      = %5.2f KB\n", (double)(dHdt_Index_Size     *sizeof(int))       / 1024.);
    if (GPUHC_type == std::string("P2C")) {
        printf("dHdx_PHC_Coeffs_Size = %5.2f KB\n", (double)(dHdx_PHC_Coeffs_Size*sizeof(magmaFloatComplex)) / 1024.);
        printf("dHdt_PHC_Coeffs_Size = %5.2f KB\n", (double)(dHdt_PHC_Coeffs_Size*sizeof(magmaFloatComplex)) / 1024.);
    }
#endif

    //> (1) First check GPU is Ampere or higher
    magma_int_t arch = magma_getdevice_arch();
    if (arch >= 800) {
        GPU_arch_Ampere_and_above = true;

        //> Try to get the L2-cache size and maximal persisting cache size
        for (int gpu_id = 0; gpu_id < Num_Of_GPUs; gpu_id++) {
            dev_id = (Num_Of_GPUs == 1) ? (SET_GPU_DEVICE_ID) : (gpu_id);
            cudaDeviceProp device_prop;
            cudacheck(cudaGetDeviceProperties( &device_prop, dev_id ));
            std::cout << "L2 Cache Size: " << device_prop.l2CacheSize / 1024 / 1024 << " MB" << std::endl;
            std::cout << "Max Persistent L2 Cache Size: " << device_prop.persistingL2CacheMaxSize / 1024 / 1024 << " MB" << std::endl;
        }

        double Num_Of_MBytes_dHdx      = (double)((dHdx_Index_Size*sizeof(int)) / 1024.) / 1024.;
        double Num_Of_MBytes_dHdt      = (double)((dHdt_Index_Size*sizeof(int)) / 1024.) / 1024.;
        Num_Of_MBytes_Persistent_Data  = Num_Of_MBytes_dHdx + Num_Of_MBytes_dHdt;
        Num_Of_MBytes_Persistent_Cache = Num_Of_MBytes_Persistent_Data + (1); //> 1 MB extra size to be cached
        cudacheck( cudaDeviceSetLimit(cudaLimitPersistingL2CacheSize, Num_Of_MBytes_Persistent_Cache * 1024 * 1024));
#if DEBUG_L2_PERSISTENT_CACHE
        printf("Size of dH/dx + dH/dt = %5.2f MB\n", Num_Of_MBytes_Persistent_Data);
        printf("L2 Persistent Cache   = %5.2f MB\n", Num_Of_MBytes_Persistent_Cache);
#endif
    }
    else {
        GPU_arch_Ampere_and_above = false;
    }

    //> Define problem file path for problem data reader and output file path for results evaluations
    Problem_File_Path = std::string("../../problems/") + HC_problem;
    RANSAC_Data_File_Path = std::string("../../RANSAC_Data/") + HC_problem + "/" + RANSAC_Dataset_Name;
    Write_Files_Path = std::string("../../") + WRITE_FILES_FOLDER;

    //> Evaluations
    // std::string 
    Evaluate_GPUHC_Sols = std::shared_ptr<Evaluations>(new Evaluations(Write_Files_Path, "GPU-HC", Num_Of_Tracks, Num_Of_Vars));
}

void GPU_HC_Solver::Allocate_Arrays() {

    //> CPU Allocations
    magma_cmalloc_cpu( &h_Start_Sols,           Num_Of_Tracks*(Num_Of_Vars+1) );
    magma_cmalloc_cpu( &h_Start_Params,         (Num_Of_Params+1) );

    h_unified_dHdx_dHdt_Index = new int[ unified_dHdx_dHdt_Index_Size ];
    h_dHdx_Index = new int[ dHdx_Index_Size ];
    h_dHdt_Index = new int[ dHdt_Index_Size ];
    
    magma_cmalloc_cpu( &h_GPU_HC_Track_Sols_Stack, Num_Of_Tracks*(Num_Of_Vars+1)*NUM_OF_RANSAC_ITERATIONS ); //> Use to store GPU results from the CPU side
    h_is_GPU_HC_Sol_Converge_Stack = new bool[ Num_Of_Tracks*NUM_OF_RANSAC_ITERATIONS ];
    h_is_GPU_HC_Sol_Infinity_Stack = new bool[ Num_Of_Tracks*NUM_OF_RANSAC_ITERATIONS ];
    h_Camera_Intrinsic_Matrix = new float[9];

    //> Allocate sizes according to sub RANSAC iterations for multiple GPUs, if necessary
    for (int gpu_id = 0; gpu_id < Num_Of_GPUs; gpu_id++) {
        dev_id = (Num_Of_GPUs == 1) ? (SET_GPU_DEVICE_ID) : (gpu_id);
        magma_setdevice(dev_id);
        magma_cmalloc_cpu( &h_Homotopy_Sols[gpu_id],            Num_Of_Tracks*(Num_Of_Vars+1)*sub_RANSAC_iters[gpu_id] );
        magma_cmalloc_cpu( &h_Target_Params[gpu_id],            (Num_Of_Params+1)*sub_RANSAC_iters[gpu_id] );
        magma_cmalloc_cpu( &h_diffParams[gpu_id],               (Num_Of_Params+1)*sub_RANSAC_iters[gpu_id] );
        magma_cmalloc_cpu( &h_GPU_HC_Track_Sols[gpu_id],        Num_Of_Tracks*(Num_Of_Vars+1)*sub_RANSAC_iters[gpu_id] ); //> Use to store GPU results from the CPU side
        magma_cmalloc_cpu( &h_Debug_Purpose[gpu_id],            Num_Of_Tracks*sub_RANSAC_iters[gpu_id] );
        magma_cmalloc_cpu( &h_dHdx_PHC_Coeffs[gpu_id],          dHdx_PHC_Coeffs_Size*sub_RANSAC_iters[gpu_id] );
        magma_cmalloc_cpu( &h_dHdt_PHC_Coeffs[gpu_id],          dHdt_PHC_Coeffs_Size*sub_RANSAC_iters[gpu_id] );
        h_is_GPU_HC_Sol_Converge[gpu_id]        = new bool[ Num_Of_Tracks*sub_RANSAC_iters[gpu_id] ];
        h_is_GPU_HC_Sol_Infinity[gpu_id]        = new bool[ Num_Of_Tracks*sub_RANSAC_iters[gpu_id] ];
    }
    
    //> GPU Allocations
    for (int gpu_id = 0; gpu_id < Num_Of_GPUs; gpu_id++) {
        dev_id = (Num_Of_GPUs == 1) ? (SET_GPU_DEVICE_ID) : (gpu_id);
        magma_setdevice(dev_id);
        magma_cmalloc( &d_Homotopy_Sols[gpu_id],                        Num_Of_Tracks*(Num_Of_Vars+1)*sub_RANSAC_iters[gpu_id] );
        magma_cmalloc( &d_Target_Params[gpu_id],                        (Num_Of_Params+1)*sub_RANSAC_iters[gpu_id] );
        magma_cmalloc( &d_diffParams[gpu_id],                           (Num_Of_Params+1)*sub_RANSAC_iters[gpu_id] );
        magma_cmalloc( &d_Debug_Purpose[gpu_id],                        Num_Of_Tracks*sub_RANSAC_iters[gpu_id] );
        magma_malloc( (void**) &d_Homotopy_Sols_array[gpu_id],          (Num_Of_Tracks+1)*(sub_RANSAC_iters[gpu_id]) * sizeof(magmaFloatComplex*) );
        magma_malloc( (void**) &d_Start_Sols_array[gpu_id],             Num_Of_Tracks * sizeof(magmaFloatComplex*) );
        cudacheck( cudaMalloc( &d_is_GPU_HC_Sol_Converge[gpu_id],       (Num_Of_Tracks*sub_RANSAC_iters[gpu_id]) * sizeof(bool) ));
        cudacheck( cudaMalloc( &d_is_GPU_HC_Sol_Infinity[gpu_id],       (Num_Of_Tracks*sub_RANSAC_iters[gpu_id]) * sizeof(bool) ));

        magma_cmalloc( &d_Start_Sols[gpu_id],                           Num_Of_Tracks*(Num_Of_Vars+1) );
        magma_cmalloc( &d_Start_Params[gpu_id],                         (Num_Of_Params+1) );
        magma_cmalloc( &d_dHdx_PHC_Coeffs[gpu_id],                      dHdx_PHC_Coeffs_Size*sub_RANSAC_iters[gpu_id] );
        magma_cmalloc( &d_dHdt_PHC_Coeffs[gpu_id],                      dHdt_PHC_Coeffs_Size*sub_RANSAC_iters[gpu_id] );
        cudacheck( cudaMalloc( &d_dHdx_Index[gpu_id],                   (dHdx_Index_Size) * sizeof(int)) );
        cudacheck( cudaMalloc( &d_dHdt_Index[gpu_id],                   (dHdt_Index_Size) * sizeof(int)) );
        cudacheck( cudaMalloc( &d_unified_dHdx_dHdt_Index[gpu_id],      (unified_dHdx_dHdt_Index_Size) * sizeof(int) ) );
    }
}

bool GPU_HC_Solver::Read_Problem_Data() {

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
    is_Data_Read_Successfully = Load_Problem_Data->Read_Start_Sols( h_Start_Sols );
    if (!is_Data_Read_Successfully) { LOG_DATA_LOAD_ERROR("Start Solutions"); return false; }

    //> (4) Start solutions used for intermediate homotopy solutions
    for (int gpu_id = 0; gpu_id < Num_Of_GPUs; gpu_id++) {
        is_Data_Read_Successfully = Load_Problem_Data->Feed_Start_Sols_for_Intermediate_Homotopy(h_Start_Sols, h_Homotopy_Sols[gpu_id], sub_RANSAC_iters[gpu_id]);
        if (!is_Data_Read_Successfully) { LOG_DATA_LOAD_ERROR("Start solutions fed to the intermediate HC solutions"); return false; }
    }
    
    //> (5) dH/dx and dH/dt evaluation indices
    is_Data_Read_Successfully = Load_Problem_Data->Read_dHdx_Indices<int>( h_dHdx_Index );
    if (!is_Data_Read_Successfully) { LOG_DATA_LOAD_ERROR("dH/dx Evaluation Indices"); return false; }

    is_Data_Read_Successfully = Load_Problem_Data->Read_dHdt_Indices<int>( h_dHdt_Index );
    if (!is_Data_Read_Successfully) { LOG_DATA_LOAD_ERROR("dH/dt Evaluation Indices"); return false; }

    Load_Problem_Data->Read_unified_dHdx_dHdt_Indices<int>( h_unified_dHdx_dHdt_Index, h_dHdx_Index, h_dHdt_Index, dHdx_Index_Size, dHdt_Index_Size );
    
    return true;
}

bool GPU_HC_Solver::Read_RANSAC_Data( int tp_index ) {
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

    //> (1) Camera extrinsic matrices
    is_Data_Read_Successfully = Load_Problem_Data->Read_Camera_Poses( h_Camera_Pose21, h_Camera_Pose31, tp_index );
    if (!is_Data_Read_Successfully) { LOG_DATA_LOAD_ERROR("Camera Extrinsic Matrices"); return false; }

    //> (2) Camera intrinsic matrix
    is_Data_Read_Successfully = Load_Problem_Data->Read_Intrinsic_Matrix( h_Camera_Intrinsic_Matrix );
    if (!is_Data_Read_Successfully) { LOG_DATA_LOAD_ERROR("Camera Intrinsic Matrices"); return false; }

    //> (3) Triplet edgel correspondences, convert data to edgel locations and tangents arrays
    Load_Problem_Data->Read_Triplet_Edgels( h_Triplet_Edge_Locations, h_Triplet_Edge_Tangents );

    return true;
}

void GPU_HC_Solver::Prepare_Target_Params( unsigned rand_seed_ ) {
    unsigned Edgel_Indx[3] = {0};
    const int offset_for_tangents = 18;
    std::array<int, 3> picked_samples;
    std::vector< std::array<int, 3> > target_params_match_indices;
#if FEED_RANDOM_SEED 
    srand (time(NULL));
#else
    srand (rand_seed_);
#endif
    //> Loop over all requested GPUs
    for (int gpu_id = 0; gpu_id < Num_Of_GPUs; gpu_id++) {

        for(int ti = 0; ti < sub_RANSAC_iters[gpu_id]; ti++) {
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
                    (h_Target_Params[gpu_id] + ti*(Num_Of_Params+1))[i*6 + j] = MAGMA_C_MAKE(h_Triplet_Edge_Locations(ei, j), 0.0);
                }
            }
            //> Tangents of the edgels
            for (int i = 0; i < 2; i++) {
                int ei = Edgel_Indx[i];
                for (int j = 0; j < 6; j++) {
                    (h_Target_Params[gpu_id] + ti*(Num_Of_Params+1))[i*6 + j + offset_for_tangents] = MAGMA_C_MAKE(h_Triplet_Edge_Tangents(ei, j), 0.0);
                }
            }
            (h_Target_Params[gpu_id] + ti*(Num_Of_Params+1))[30] = MAGMA_C_MAKE(1.0, 0.0);
            (h_Target_Params[gpu_id] + ti*(Num_Of_Params+1))[31] = MAGMA_C_MAKE(0.5, 0.0);
            (h_Target_Params[gpu_id] + ti*(Num_Of_Params+1))[32] = MAGMA_C_MAKE(1.0, 0.0);
            (h_Target_Params[gpu_id] + ti*(Num_Of_Params+1))[33] = MAGMA_C_ONE;

            //> (2) Compute the difference of the start and target parameters in h_diffParams
            for (int i = 0; i <= Num_Of_Params; i++) 
                (h_diffParams[gpu_id] + ti*(Num_Of_Params+1))[i] = (h_Target_Params[gpu_id] + ti*(Num_Of_Params+1))[i] - (h_Start_Params)[i];

            //> (3) Coefficients from parameters, if required (only used for P2C mode)
            if (Use_P2C) {
                if (!Load_Problem_Data->Construct_Coeffs_From_Params
                    ( HC_problem, h_Target_Params[gpu_id] + ti*(Num_Of_Params+1), h_Start_Params, \
                                  h_dHdx_PHC_Coeffs[gpu_id] + ti*(dHdx_PHC_Coeffs_Size), \
                                  h_dHdt_PHC_Coeffs[gpu_id] + ti*(dHdt_PHC_Coeffs_Size)) ) {
                    LOG_ERROR("Failed to prepare coefficients from parameters (Use_P2C)");
                    exit(1);
                }
            }
        }
    }

//> TODO: Reflect multiple GPUs for this debugging part 
#if RANSAC_DEBUG
    int sample_index = 0;
    if (NUM_OF_RANSAC_ITERATIONS > sample_index) Load_Problem_Data->Print_Out_Target_Params_from_Triplet_Edgels(sample_index, target_params_match_indices, h_Target_Params);
    else Load_Problem_Data->Print_Out_Target_Params_from_Triplet_Edgels(0, target_params_match_indices, h_Target_Params);
#endif
}

void GPU_HC_Solver::Set_RANSAC_Abort_Arrays() {

    //> Allocate arrays and assign values to arrays for initialization
    if (Abort_RANSAC_by_Good_Sol) {
        for (int gpu_id = 0; gpu_id < Num_Of_GPUs; gpu_id++) {
            dev_id = (Num_Of_GPUs == 1) ? (SET_GPU_DEVICE_ID) : (gpu_id);
            magma_setdevice(dev_id);

            //> CPU side
            h_Found_Trifocal_Sols[gpu_id]           = new bool[ 1 ];
            h_Trifocal_Sols_Batch_Index[gpu_id]     = new int[ Num_Of_Tracks * sub_RANSAC_iters[gpu_id] ];

            //> Assign zeros for h_Found_Trifocal_Sols
            h_Found_Trifocal_Sols[gpu_id][0] = 0;

            //> Assign -1 to batch index of trifocal solutions
            initialize_trifocal_sols_batch_index( Num_Of_Tracks*sub_RANSAC_iters[gpu_id], h_Trifocal_Sols_Batch_Index[gpu_id] );

            //> GPU side
            cudacheck( cudaMalloc( &d_Intrinsic_Matrix[gpu_id],             (9) * sizeof(float) ) );
            cudacheck( cudaMalloc( &d_Triplet_Edge_Locations[gpu_id],       (Num_Of_Triplet_Edgels*2*3) * sizeof(float) ) );
            cudacheck( cudaMalloc( &d_Found_Trifocal_Sols[gpu_id],          (1) * sizeof(bool) ) );
            cudacheck( cudaMalloc( &d_Trifocal_Sols_Batch_Index[gpu_id],    (Num_Of_Tracks*sub_RANSAC_iters[gpu_id]) * sizeof(int) ) );
        }
    }
    
}

void GPU_HC_Solver::Data_Transfer_From_Host_To_Device() {

    //> Loop over all requested GPUs
    for (int gpu_id = 0; gpu_id < Num_Of_GPUs; gpu_id++) {
        dev_id = (Num_Of_GPUs == 1) ? (SET_GPU_DEVICE_ID) : (gpu_id);
        magma_setdevice(dev_id);

        magma_csetmatrix( Num_Of_Vars+1,   Num_Of_Tracks*sub_RANSAC_iters[gpu_id],  h_Homotopy_Sols[gpu_id],  (Num_Of_Vars+1),  d_Homotopy_Sols[gpu_id], Num_Of_Vars+1,     gpu_queues[gpu_id] );
        magma_csetmatrix( Num_Of_Vars+1,   Num_Of_Tracks,                           h_Start_Sols,     (Num_Of_Vars+1),  d_Start_Sols[gpu_id],    Num_Of_Vars+1,     gpu_queues[gpu_id] );
        magma_csetmatrix( Num_Of_Params+1, (1),                                     h_Start_Params,   Num_Of_Params+1,  d_Start_Params[gpu_id],  Num_Of_Params+1,   gpu_queues[gpu_id] );
        magma_csetmatrix( (Num_Of_Params+1)*sub_RANSAC_iters[gpu_id], (1), h_diffParams[gpu_id],    (Num_Of_Params+1)*sub_RANSAC_iters[gpu_id], d_diffParams[gpu_id],    (Num_Of_Params+1)*sub_RANSAC_iters[gpu_id], gpu_queues[gpu_id] );
        magma_csetmatrix( (Num_Of_Params+1)*sub_RANSAC_iters[gpu_id], (1), h_Target_Params[gpu_id], (Num_Of_Params+1)*sub_RANSAC_iters[gpu_id], d_Target_Params[gpu_id], (Num_Of_Params+1)*sub_RANSAC_iters[gpu_id], gpu_queues[gpu_id] );
        cudacheck( cudaMemcpy( d_dHdx_Index[gpu_id], h_dHdx_Index,     dHdx_Index_Size * sizeof(int), cudaMemcpyHostToDevice) );
        cudacheck( cudaMemcpy( d_dHdt_Index[gpu_id], h_dHdt_Index,     dHdt_Index_Size * sizeof(int), cudaMemcpyHostToDevice) );
        cudacheck( cudaMemcpy( d_unified_dHdx_dHdt_Index[gpu_id], h_unified_dHdx_dHdt_Index, unified_dHdx_dHdt_Index_Size * sizeof(int), cudaMemcpyHostToDevice) );
        
        //> connect pointer to 2d arrays
        magma_cset_pointer( d_Start_Sols_array[gpu_id],    d_Start_Sols[gpu_id],     (Num_Of_Vars+1), 0, 0, (Num_Of_Vars+1), Num_Of_Tracks, gpu_queues[gpu_id] );
        magma_cset_pointer( d_Homotopy_Sols_array[gpu_id], d_Homotopy_Sols[gpu_id],  (Num_Of_Vars+1), 0, 0, (Num_Of_Vars+1), Num_Of_Tracks*sub_RANSAC_iters[gpu_id], gpu_queues[gpu_id] );

        if (Abort_RANSAC_by_Good_Sol) {
            cudacheck( cudaMemcpy( d_Triplet_Edge_Locations[gpu_id], h_Triplet_Edge_Locations, (Num_Of_Triplet_Edgels*2*3) * sizeof(float), cudaMemcpyHostToDevice ) );
            cudacheck( cudaMemcpy( d_Intrinsic_Matrix[gpu_id], h_Camera_Intrinsic_Matrix, (9) * sizeof(float), cudaMemcpyHostToDevice ) );
            cudacheck( cudaMemcpy( d_Trifocal_Sols_Batch_Index[gpu_id], h_Trifocal_Sols_Batch_Index[gpu_id], (Num_Of_Tracks*sub_RANSAC_iters[gpu_id])*sizeof(int), cudaMemcpyHostToDevice ) );
            cudacheck( cudaMemcpy( d_Found_Trifocal_Sols[gpu_id], h_Found_Trifocal_Sols[gpu_id], (1)*sizeof(bool), cudaMemcpyHostToDevice ) );

            magma_csetmatrix( dHdx_PHC_Coeffs_Size*sub_RANSAC_iters[gpu_id], (1), h_dHdx_PHC_Coeffs[gpu_id], dHdx_PHC_Coeffs_Size*sub_RANSAC_iters[gpu_id], d_dHdx_PHC_Coeffs[gpu_id], dHdx_PHC_Coeffs_Size*sub_RANSAC_iters[gpu_id], gpu_queues[gpu_id] );
            magma_csetmatrix( dHdt_PHC_Coeffs_Size*sub_RANSAC_iters[gpu_id], (1), h_dHdt_PHC_Coeffs[gpu_id], dHdt_PHC_Coeffs_Size*sub_RANSAC_iters[gpu_id], d_dHdt_PHC_Coeffs[gpu_id], dHdt_PHC_Coeffs_Size*sub_RANSAC_iters[gpu_id], gpu_queues[gpu_id] );
        }
    }
}

void GPU_HC_Solver::Set_CUDA_Stream_Attributes() {
    if (GPU_arch_Ampere_and_above) {
        //> Set stream attribute
        for (int gpu_id = 0; gpu_id < Num_Of_GPUs; gpu_id++) {
            cudaStreamAttrValue stream_attribute_non_thrashing;
            stream_attribute_non_thrashing.accessPolicyWindow.base_ptr  = reinterpret_cast<void*>(d_unified_dHdx_dHdt_Index[gpu_id]);
            stream_attribute_non_thrashing.accessPolicyWindow.num_bytes = Num_Of_MBytes_Persistent_Data * 1024 * 1024;
            stream_attribute_non_thrashing.accessPolicyWindow.hitRatio  = 1.0;
            stream_attribute_non_thrashing.accessPolicyWindow.hitProp   = cudaAccessPropertyPersisting;
            stream_attribute_non_thrashing.accessPolicyWindow.missProp  = cudaAccessPropertyStreaming;

            cudacheck( cudaStreamSetAttribute( gpu_queues[gpu_id]->cuda_stream(), cudaStreamAttributeAccessPolicyWindow, &stream_attribute_non_thrashing) );
        }
    }
}

void GPU_HC_Solver::Solve_by_GPU_HC() {
    std::cout << "GPU computing ..." << std::endl << std::endl;

    //> Get kernel version number (coresponding to the google spreadsheet)
    kernel_version = get_kernel_version_number();
    
    //> start the timer
    multi_GPUs_time = magma_wtime();

    //> For now we only solve the trifocal 2op1p 30x30 problem
    if (HC_problem == "trifocal_2op1p_30x30") {

        //> Distribute workload across multiple GPUs
        for (int gpu_id = 0; gpu_id < Num_Of_GPUs; gpu_id++) {
            dev_id = (Num_Of_GPUs == 1) ? (SET_GPU_DEVICE_ID) : (gpu_id);
            magma_setdevice(dev_id);

            switch(kernel_version) {
                case 1: //> (#1) PH-(x), CodeOpt-(x), TrunPaths-(x), TrunRANSAC-(x) (Naive GPU-HC approach)
                    gpu_time[gpu_id] = kernel_GPUHC_trifocal_2op1p_30x30_P2C \
                        (gpu_queues[gpu_id], sub_RANSAC_iters[gpu_id], \
                         GPUHC_Max_Steps, GPUHC_Max_Correction_Steps, GPUHC_delta_t_incremental_steps, \
                         d_Start_Sols_array[gpu_id], d_Homotopy_Sols_array[gpu_id], \
                         d_dHdx_Index[gpu_id],       d_dHdt_Index[gpu_id], \
                         d_dHdx_PHC_Coeffs[gpu_id],  d_dHdt_PHC_Coeffs[gpu_id], \
                         d_is_GPU_HC_Sol_Converge[gpu_id], d_is_GPU_HC_Sol_Infinity[gpu_id], d_Debug_Purpose[gpu_id] );
                    break;
                case 2: //> (#2) PH-(v), CodeOpt-(x), TrunPaths-(x), TrunRANSAC-(x)
                    gpu_time[gpu_id] = kernel_GPUHC_trifocal_2op1p_30x30_PH \
                        (gpu_queues[gpu_id], sub_RANSAC_iters[gpu_id], \
                         GPUHC_Max_Steps, GPUHC_Max_Correction_Steps, GPUHC_delta_t_incremental_steps, \
                         d_Start_Sols_array[gpu_id],  d_Homotopy_Sols_array[gpu_id], \
                         d_Start_Params[gpu_id],      d_Target_Params[gpu_id], d_diffParams[gpu_id], \
                         d_dHdx_Index[gpu_id],        d_dHdt_Index[gpu_id], \
                         d_is_GPU_HC_Sol_Converge[gpu_id], d_is_GPU_HC_Sol_Infinity[gpu_id], d_Debug_Purpose[gpu_id] );
                    break;
                case 3: //> (#3) PH-(v), CodeOpt-(v), TrunPaths-(x), TrunRANSAC-(x)
                    if (GPU_arch_Ampere_and_above) {
                        gpu_time[gpu_id] = kernel_GPUHC_trifocal_2op1p_30x30_PH_CodeOpt \
                            (gpu_queues[gpu_id], sub_RANSAC_iters[gpu_id], \
                            GPUHC_Max_Steps, GPUHC_Max_Correction_Steps, GPUHC_delta_t_incremental_steps, \
                            d_Start_Sols_array[gpu_id],  d_Homotopy_Sols_array[gpu_id], \
                            d_Start_Params[gpu_id],      d_Target_Params[gpu_id], d_diffParams[gpu_id], \
                            d_unified_dHdx_dHdt_Index[gpu_id], \
                            d_is_GPU_HC_Sol_Converge[gpu_id], d_is_GPU_HC_Sol_Infinity[gpu_id], d_Debug_Purpose[gpu_id] );
                    }
                    else {
                        gpu_time[gpu_id] = kernel_GPUHC_trifocal_2op1p_30x30_PH_CodeOpt_Volta \
                            (gpu_queues[gpu_id], sub_RANSAC_iters[gpu_id], \
                            GPUHC_Max_Steps, GPUHC_Max_Correction_Steps, GPUHC_delta_t_incremental_steps, \
                            d_Start_Sols_array[gpu_id],  d_Homotopy_Sols_array[gpu_id], \
                            d_Start_Params[gpu_id],      d_Target_Params[gpu_id], d_diffParams[gpu_id], \
                            d_dHdx_Index[gpu_id],        d_dHdt_Index[gpu_id], \
                            d_is_GPU_HC_Sol_Converge[gpu_id], d_is_GPU_HC_Sol_Infinity[gpu_id], d_Debug_Purpose[gpu_id] );
                    }
                    break;
                case 4: //> (#4) PH-(v), CodeOpt-(v), TrunPaths-(v), TrunRANSAC-(x)
                    if (GPU_arch_Ampere_and_above) {
                        gpu_time[gpu_id] = kernel_GPUHC_trifocal_2op1p_30x30_PH_CodeOpt_TrunPaths \
                            (gpu_queues[gpu_id], sub_RANSAC_iters[gpu_id], \
                            GPUHC_Max_Steps, GPUHC_Max_Correction_Steps, GPUHC_delta_t_incremental_steps, \
                            d_Start_Sols_array[gpu_id], d_Homotopy_Sols_array[gpu_id], \
                            d_Start_Params[gpu_id],     d_Target_Params[gpu_id], d_diffParams[gpu_id], 
                            d_unified_dHdx_dHdt_Index[gpu_id], \
                            d_is_GPU_HC_Sol_Converge[gpu_id], d_is_GPU_HC_Sol_Infinity[gpu_id], d_Debug_Purpose[gpu_id] );
                    }
                    else {
                        gpu_time[gpu_id] = kernel_GPUHC_trifocal_2op1p_30x30_PH_CodeOpt_TrunPaths_Volta \
                            (gpu_queues[gpu_id], sub_RANSAC_iters[gpu_id], \
                            GPUHC_Max_Steps, GPUHC_Max_Correction_Steps, GPUHC_delta_t_incremental_steps, \
                            d_Start_Sols_array[gpu_id], d_Homotopy_Sols_array[gpu_id], \
                            d_Start_Params[gpu_id],     d_Target_Params[gpu_id], d_diffParams[gpu_id], 
                            d_dHdx_Index[gpu_id],       d_dHdt_Index[gpu_id], \
                            d_is_GPU_HC_Sol_Converge[gpu_id], d_is_GPU_HC_Sol_Infinity[gpu_id], d_Debug_Purpose[gpu_id] );
                    }
                    break;
                case 5: //> (#5) PH-(v), CodeOpt-(v), TrunPaths-(v), TrunRANSAC-(v)
                    if (GPU_arch_Ampere_and_above) {
                        gpu_time[gpu_id] = kernel_GPUHC_trifocal_2op1p_30x30_PH_CodeOpt_TrunPaths_TrunRANSAC \
                            (gpu_queues[gpu_id], sub_RANSAC_iters[gpu_id], \
                            Num_Of_Triplet_Edgels, GPUHC_Max_Steps, GPUHC_Max_Correction_Steps, GPUHC_delta_t_incremental_steps, \
                            d_Start_Sols_array[gpu_id], d_Homotopy_Sols_array[gpu_id], \
                            d_Start_Params[gpu_id],     d_Target_Params[gpu_id], d_diffParams[gpu_id], 
                            d_unified_dHdx_dHdt_Index[gpu_id], d_Triplet_Edge_Locations[gpu_id], d_Intrinsic_Matrix[gpu_id], \
                            d_is_GPU_HC_Sol_Converge[gpu_id], d_is_GPU_HC_Sol_Infinity[gpu_id], d_Debug_Purpose[gpu_id], \
                            d_Found_Trifocal_Sols[gpu_id], d_Trifocal_Sols_Batch_Index[gpu_id] );
                    }
                    else {
                        gpu_time[gpu_id] = kernel_GPUHC_trifocal_2op1p_30x30_PH_CodeOpt_TrunPaths_TrunRANSAC_Volta \
                            (gpu_queues[gpu_id], sub_RANSAC_iters[gpu_id], \
                            Num_Of_Triplet_Edgels, GPUHC_Max_Steps, GPUHC_Max_Correction_Steps, GPUHC_delta_t_incremental_steps, \
                            d_Start_Sols_array[gpu_id], d_Homotopy_Sols_array[gpu_id], \
                            d_Start_Params[gpu_id],     d_Target_Params[gpu_id], d_diffParams[gpu_id], 
                            d_dHdx_Index[gpu_id],       d_dHdt_Index[gpu_id], \
                            d_Triplet_Edge_Locations[gpu_id], d_Intrinsic_Matrix[gpu_id], \
                            d_is_GPU_HC_Sol_Converge[gpu_id], d_is_GPU_HC_Sol_Infinity[gpu_id], d_Debug_Purpose[gpu_id], \
                            d_Found_Trifocal_Sols[gpu_id], d_Trifocal_Sols_Batch_Index[gpu_id] );
                    }
                    break;
                default:
                    break;
            }   //> switch statement
        }   //> for-loop over gpu_ids
    }   //> if HC_problem statement

    //> Sync across all GPUs before measuring elapsed time
    for(magma_int_t gpu_id = 0; gpu_id < Num_Of_GPUs; gpu_id++) {
        dev_id = (Num_Of_GPUs == 1) ? (SET_GPU_DEVICE_ID) : (gpu_id);
        magma_setdevice(dev_id);
        magma_queue_sync( gpu_queues[gpu_id] );
    }
    //> End the timer
    multi_GPUs_time = magma_wtime() - multi_GPUs_time;

    //> Check returns from the GPU kernel
    for (int gpu_id = 0; gpu_id < Num_Of_GPUs; gpu_id++) {
        dev_id = (Num_Of_GPUs == 1) ? (SET_GPU_DEVICE_ID) : (gpu_id);
        magma_setdevice(dev_id);
        magma_cgetmatrix( (Num_Of_Vars+1), Num_Of_Tracks*sub_RANSAC_iters[gpu_id], d_Homotopy_Sols[gpu_id],  (Num_Of_Vars+1), h_GPU_HC_Track_Sols[gpu_id],    (Num_Of_Vars+1), gpu_queues[gpu_id] );
        cudacheck( cudaMemcpy( h_is_GPU_HC_Sol_Converge[gpu_id], d_is_GPU_HC_Sol_Converge[gpu_id], Num_Of_Tracks*sub_RANSAC_iters[gpu_id]*sizeof(bool), cudaMemcpyDeviceToHost) );
        cudacheck( cudaMemcpy( h_is_GPU_HC_Sol_Infinity[gpu_id], d_is_GPU_HC_Sol_Infinity[gpu_id], Num_Of_Tracks*sub_RANSAC_iters[gpu_id]*sizeof(bool), cudaMemcpyDeviceToHost) );
        
        if (Abort_RANSAC_by_Good_Sol) {
            cudacheck( cudaMemcpy( h_Found_Trifocal_Sols[gpu_id], d_Found_Trifocal_Sols[gpu_id], (1)*sizeof(bool), cudaMemcpyDeviceToHost) );
            cudacheck( cudaMemcpy( h_Trifocal_Sols_Batch_Index[gpu_id], d_Trifocal_Sols_Batch_Index[gpu_id], (Num_Of_Tracks*sub_RANSAC_iters[gpu_id])*sizeof(int), cudaMemcpyDeviceToHost) );
        }
    }

    //> For early aborting RANSAC process debugging...
    if (Abort_RANSAC_by_Good_Sol) {
        for (magma_int_t gpu_id = 0; gpu_id < Num_Of_GPUs; gpu_id++) {
            std::string found_sol_str = (h_Found_Trifocal_Sols[gpu_id][0]) ? "Yes" : "No";
            std::cout << "GPU id " << gpu_id << " found solution? " << found_sol_str;
            if (h_Found_Trifocal_Sols[gpu_id][0]) {
                std::cout << " / Found solution batch index = ";
                for (int i = 0; i < Num_Of_Tracks*sub_RANSAC_iters[gpu_id]; i++) {
                    if ( (h_Trifocal_Sols_Batch_Index[gpu_id])[i] != -1) {
                        candidate_batch_ids.push_back( h_Trifocal_Sols_Batch_Index[gpu_id][i] );
                        std::cout << h_Trifocal_Sols_Batch_Index[gpu_id][i] << ", ";
                    }
                }
            }
            std::cout << std::endl;
        }
    }
    

#if GPU_DEBUG
    magma_cgetmatrix( Num_Of_Tracks, NUM_OF_RANSAC_ITERATIONS, d_Debug_Purpose, Num_Of_Tracks, h_Debug_Purpose, Num_Of_Tracks, gpu_queues[0] );
#endif

    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << "## Solving " << HC_print_problem_name << std::endl;
    std::cout << " - Mode: ";
    print_kernel_mode();
    std::cout << std::endl;

    //> Print out timings
    printf("## Timings:\n");
    printf(" - GPU Computation Time = %7.2f (ms)\n", (multi_GPUs_time)*1000);

    //> First stack all results
    int offset_stack = 0;
    for (int gpu_id = 0; gpu_id < Num_Of_GPUs; gpu_id++) {
        memcpy(h_GPU_HC_Track_Sols_Stack + offset_stack, h_GPU_HC_Track_Sols[gpu_id], Num_Of_Tracks*(Num_Of_Vars+1)*sub_RANSAC_iters[gpu_id]*sizeof(magmaFloatComplex));
        offset_stack += Num_Of_Tracks*(Num_Of_Vars+1)*sub_RANSAC_iters[gpu_id];
    }

    offset_stack = 0;
    for (int gpu_id = 0; gpu_id < Num_Of_GPUs; gpu_id++) {
        memcpy(h_is_GPU_HC_Sol_Converge_Stack + offset_stack, h_is_GPU_HC_Sol_Converge[gpu_id], Num_Of_Tracks*sub_RANSAC_iters[gpu_id]*sizeof(bool));
        memcpy(h_is_GPU_HC_Sol_Infinity_Stack + offset_stack, h_is_GPU_HC_Sol_Infinity[gpu_id], Num_Of_Tracks*sub_RANSAC_iters[gpu_id]*sizeof(bool));
        offset_stack += Num_Of_Tracks*sub_RANSAC_iters[gpu_id];
    }

    //> Object for the Evaluations class
#if WRITE_GPUHC_CONVERGED_SOLS
    Evaluate_GPUHC_Sols->Write_Converged_Sols( h_GPU_HC_Track_Sols[0], h_is_GPU_HC_Sol_Converge[0] );
#endif
    Evaluate_GPUHC_Sols->Evaluate_RANSAC_HC_Sols( h_GPU_HC_Track_Sols_Stack, h_is_GPU_HC_Sol_Converge_Stack, h_is_GPU_HC_Sol_Infinity_Stack );
    // Evaluate_GPUHC_Sols->Find_Unique_Sols( h_GPU_HC_Track_Sols, h_is_GPU_HC_Sol_Converge );

    //> Print out evaluation results
    std::cout << "\n## Evaluation of GPU-HC Solutions: "      << std::endl;
    std::cout << " - Number of Converged Solutions:       " << Evaluate_GPUHC_Sols->Num_Of_Coverged_Sols << std::endl;
    std::cout << " - Number of Real Solutions:            " << Evaluate_GPUHC_Sols->Num_Of_Real_Sols << std::endl;
    std::cout << " - Number of Infinity Failed Solutions: " << Evaluate_GPUHC_Sols->Num_Of_Inf_Sols << std::endl;
    // std::cout << " - Number of Unique Solutions:          " << Evaluate_GPUHC_Sols->Num_Of_Unique_Sols << std::endl;

    Collect_Num_Of_Coverged_Sols.push_back( Evaluate_GPUHC_Sols->Num_Of_Coverged_Sols );
    Collect_Num_Of_Inf_Sols.push_back( Evaluate_GPUHC_Sols->Num_Of_Real_Sols );
    Collect_Num_Of_Real_Sols.push_back( Evaluate_GPUHC_Sols->Num_Of_Inf_Sols );

    Evaluate_GPUHC_Sols->Transform_GPUHC_Sols_to_Trifocal_Relative_Pose( h_GPU_HC_Track_Sols_Stack, h_is_GPU_HC_Sol_Converge_Stack, h_Camera_Intrinsic_Matrix );
    bool find_good_sol = Evaluate_GPUHC_Sols->get_Solution_with_Maximal_Support( Num_Of_Triplet_Edgels, h_Triplet_Edge_Locations, h_Triplet_Edge_Tangents, h_Camera_Intrinsic_Matrix );

    //> Reset the number of solutions
    Evaluate_GPUHC_Sols->Flush_Out_Data();

#if DEBUG_EARLY_RANSAC_ABORT
    if (Abort_RANSAC_by_Good_Sol) {
        Evaluate_GPUHC_Sols->Check_Deviations_of_Veridical_Sol_from_GT( h_GPU_HC_Track_Sols_Stack + (Num_Of_Vars+1)*candidate_batch_ids[0], h_Camera_Pose21, h_Camera_Pose31 );
        candidate_batch_ids.clear();
    }
#endif   
    // if (find_good_sol)
    //     Evaluate_GPUHC_Sols->Measure_Relative_Pose_Error( h_Camera_Pose21, h_Camera_Pose31, h_Debug_Purpose );

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

// template< typename int >
// void GPU_HC_Solver<int>::Export_Data() {
//     Evaluate_GPUHC_Sols->Write_HC_Steps_of_Actual_Solutions( GPUHC_Actual_Sols_Steps_Collections );
// }

void GPU_HC_Solver::Free_Arrays_for_Aborting_RANSAC() {
    if (Abort_RANSAC_by_Good_Sol) {
        for (int gpu_id = 0; gpu_id < Num_Of_GPUs; gpu_id++) {
            dev_id = (Num_Of_GPUs == 1) ? (SET_GPU_DEVICE_ID) : (gpu_id);
            magma_setdevice(dev_id);

            magma_free_cpu( h_Found_Trifocal_Sols[gpu_id] );
            magma_free_cpu( h_Trifocal_Sols_Batch_Index[gpu_id] );
            magma_free( d_Triplet_Edge_Locations[gpu_id] );
            magma_free( d_Found_Trifocal_Sols[gpu_id] );
            magma_free( d_Trifocal_Sols_Batch_Index[gpu_id] );
            magma_free( d_Intrinsic_Matrix[gpu_id] );
        }
    }
}

void GPU_HC_Solver::Free_Triplet_Edgels_Mem() {
    delete [] h_Triplet_Edge_Locations;
    delete [] h_Triplet_Edge_Tangents;
}

GPU_HC_Solver::~GPU_HC_Solver() {

    // delete [] devices;
    // magma_queue_destroy( gpu_queues[0] );
    //> Destroy queues
	for(magma_int_t gpu_id = 0; gpu_id < Num_Of_GPUs; gpu_id++) {
        dev_id = (Num_Of_GPUs == 1) ? (SET_GPU_DEVICE_ID) : (gpu_id);
	    magma_setdevice(dev_id);
	    magma_queue_destroy( gpu_queues[gpu_id] );
	}

    for (int gpu_id = 0; gpu_id < Num_Of_GPUs; gpu_id++) {
        delete [] h_is_GPU_HC_Sol_Converge[gpu_id];
        delete [] h_is_GPU_HC_Sol_Infinity[gpu_id];
        magma_free_cpu( h_Homotopy_Sols[gpu_id] );
        magma_free_cpu( h_Target_Params[gpu_id] );
        magma_free_cpu( h_diffParams[gpu_id] );
        magma_free_cpu( h_GPU_HC_Track_Sols[gpu_id] );
        magma_free_cpu( h_Debug_Purpose[gpu_id] );
        magma_free_cpu( h_dHdx_PHC_Coeffs[gpu_id] );
        magma_free_cpu( h_dHdt_PHC_Coeffs[gpu_id] );
    }
    delete [] h_is_GPU_HC_Sol_Converge_Stack;
    delete [] h_is_GPU_HC_Sol_Infinity_Stack;
    delete [] h_unified_dHdx_dHdt_Index;
    delete [] h_dHdx_Index;
    delete [] h_dHdt_Index;
    delete [] h_Camera_Intrinsic_Matrix;

    magma_free_cpu( h_Start_Sols );
    magma_free_cpu( h_Start_Params );
    magma_free_cpu( h_GPU_HC_Track_Sols_Stack );
    
    for (int gpu_id = 0; gpu_id < Num_Of_GPUs; gpu_id++) {
        dev_id = (Num_Of_GPUs == 1) ? (SET_GPU_DEVICE_ID) : (gpu_id);
        magma_setdevice(dev_id);
        magma_free( d_diffParams[gpu_id] );
        magma_free( d_is_GPU_HC_Sol_Converge[gpu_id] );
        magma_free( d_is_GPU_HC_Sol_Infinity[gpu_id] );
        magma_free( d_Start_Sols[gpu_id] );
        magma_free( d_Homotopy_Sols[gpu_id] );
        magma_free( d_Start_Params[gpu_id] );
        magma_free( d_Target_Params[gpu_id] );
        magma_free( d_dHdx_PHC_Coeffs[gpu_id] );
        magma_free( d_dHdt_PHC_Coeffs[gpu_id] );
        magma_free( d_Debug_Purpose[gpu_id] );
        magma_free( d_Homotopy_Sols_array[gpu_id] );
        magma_free( d_Start_Sols_array[gpu_id] );
        magma_free( d_unified_dHdx_dHdt_Index[gpu_id] );        
        cudacheck( cudaFree( d_dHdx_Index[gpu_id] ) );
        cudacheck( cudaFree( d_dHdt_Index[gpu_id] ) );
    }

    fflush( stdout );
    printf( "\n" );
    magma_finalize();
}

#endif
