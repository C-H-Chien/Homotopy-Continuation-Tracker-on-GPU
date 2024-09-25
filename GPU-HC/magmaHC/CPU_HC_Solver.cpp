#ifndef CPU_HC_Solver_cpp
#define CPU_HC_Solver_cpp
// =============================================================================================================================
//
// ChangLogs
//    24-07-09:   Initially Created (Copied from CPU_HC_Solver.cpp)
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

//> MAGMA
#include "magma_v2.h"
#include "magma_lapack.h"
#include "magma_internal.h"
#undef max
#undef min

#include "CPU_HC_Solver.hpp"
#include "definitions.hpp"

//> Constructor
CPU_HC_Solver::CPU_HC_Solver(YAML::Node Problem_Settings_File)
    : Problem_Setting_YAML_File(Problem_Settings_File)
{
    //> Parse data from the YAML file
    //> (1) Problem Name
    HC_problem                      = Problem_Setting_YAML_File["problem_name"].as<std::string>();
    HC_print_problem_name           = Problem_Setting_YAML_File["problem_print_out_name"].as<std::string>();
    //> (2) Configurations, same as the GPU-HC configurations
    CPUHC_Max_Steps                 = Problem_Setting_YAML_File["GPUHC_Max_Steps"].as<int>();
    CPUHC_Max_Correction_Steps      = Problem_Setting_YAML_File["GPUHC_Max_Correction_Steps"].as<int>();
    CPUHC_delta_t_incremental_steps = Problem_Setting_YAML_File["GPUHC_Num_Of_Steps_to_Increase_Delta_t"].as<int>();
    //> (3) Problem Specifications
    Num_Of_Vars                     = Problem_Setting_YAML_File["Num_Of_Vars"].as<int>();
    Num_Of_Params                   = Problem_Setting_YAML_File["Num_Of_Params"].as<int>();
    Num_Of_Tracks                   = Problem_Setting_YAML_File["Num_Of_Tracks"].as<int>();
    //> (4) RANSAC data
    RANSAC_Dataset_Name             = Problem_Setting_YAML_File["RANSAC_Dataset"].as<std::string>();
    //> (5) CPU-HC Settings
    Num_Of_CPU_Cores                = Problem_Setting_YAML_File["Num_Of_Cores"].as<int>();

    //> Define problem file path for problem data reader and output file path for results evaluations
    Problem_File_Path = std::string("../../problems/") + HC_problem;
    RANSAC_Data_File_Path = std::string("../../RANSAC_Data/") + HC_problem + "/" + RANSAC_Dataset_Name;
    Write_Files_Path = std::string("../../") + WRITE_FILES_FOLDER;

    //> Evaluations
    Evaluate_CPUHC_Sols = std::shared_ptr<Evaluations>(new Evaluations(Write_Files_Path, "CPU-HC", Num_Of_Tracks, Num_Of_Vars));

    //> Initialization
    magma_init();
    Num_Of_Inf_Failed_Sols = 0;
    Num_Of_Successful_Sols = 0;
    RANSAC_Sol_Offset = Num_Of_Vars * Num_Of_Tracks;
}

void CPU_HC_Solver::Allocate_Arrays() {

    //> CPU Allocations
    magma_cmalloc_cpu( &h_Start_Sols,               Num_Of_Tracks*(Num_Of_Vars+1) );
    
    magma_cmalloc_cpu( &h_Start_Params,             (Num_Of_Params+1) );
    
    //> Allocate sizes according to RANSAC iterations
    magma_cmalloc_cpu( &h_Intermediate_Sols,    Num_Of_Tracks*(Num_Of_Vars+1)*NUM_OF_RANSAC_ITERATIONS );
    magma_cmalloc_cpu( &h_CPU_HC_Track_Sols,    Num_Of_Tracks*(Num_Of_Vars+1)*NUM_OF_RANSAC_ITERATIONS );
    magma_cmalloc_cpu( &h_Track_Last_Success,   Num_Of_Tracks*(Num_Of_Vars+1)*NUM_OF_RANSAC_ITERATIONS );
    magma_cmalloc_cpu( &h_Target_Params,        Num_Of_Params*NUM_OF_RANSAC_ITERATIONS );

    //> The following arrays do not neccessary to be allocated with the scale of number of tracks and number of RANSAC iterations,
    //  but for the purpose of OpenMP parallelism where each thread occupies individual values, these arrays need enough space to avoid race condition
    magma_cmalloc_cpu( &h_cgesvA,               Num_Of_Vars*Num_Of_Vars*Num_Of_Tracks*NUM_OF_RANSAC_ITERATIONS );
    magma_cmalloc_cpu( &h_cgesvB,               Num_Of_Vars*Num_Of_Tracks*NUM_OF_RANSAC_ITERATIONS );
    magma_imalloc_cpu( &ipiv,                   Num_Of_Vars*Num_Of_Tracks*NUM_OF_RANSAC_ITERATIONS );

    h_is_Track_Converged    = new bool[ Num_Of_Tracks*NUM_OF_RANSAC_ITERATIONS ];
    h_is_Track_Inf_Failed   = new bool[ Num_Of_Tracks*NUM_OF_RANSAC_ITERATIONS ];
}

bool CPU_HC_Solver::Read_Problem_Data() {

    //> Load problem data to arrays
    bool is_Data_Read_Successfully = false;

    //> Data reader
    Load_Problem_Data = std::shared_ptr<Data_Reader>(new Data_Reader(Problem_File_Path, RANSAC_Data_File_Path, Num_Of_Tracks, Num_Of_Vars, Num_Of_Params, false));

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
    is_Data_Read_Successfully = Load_Problem_Data->Feed_Start_Sols_for_Intermediate_Homotopy(h_Start_Sols, h_Intermediate_Sols, NUM_OF_RANSAC_ITERATIONS);
    if (!is_Data_Read_Successfully) { LOG_DATA_LOAD_ERROR("Start solutions fed to the intermediate HC solutions"); return false; }

    //> (5) For CPU-HC, h_CPU_HC_Track_Sols is initialized as start solutions
    is_Data_Read_Successfully = Load_Problem_Data->Feed_Start_Sols_for_Intermediate_Homotopy(h_Start_Sols, h_CPU_HC_Track_Sols, NUM_OF_RANSAC_ITERATIONS);
    if (!is_Data_Read_Successfully) { LOG_DATA_LOAD_ERROR("Start solutions fed to h_CPU_HC_Track_Sols"); return false; }

    //> (6) This is the same for h_Track_Last_Success
    is_Data_Read_Successfully = Load_Problem_Data->Feed_Start_Sols_for_Intermediate_Homotopy(h_Start_Sols, h_Track_Last_Success, NUM_OF_RANSAC_ITERATIONS);
    if (!is_Data_Read_Successfully) { LOG_DATA_LOAD_ERROR("Start solutions fed to h_Track_Last_Success"); return false; }

    return true;
}

//> This is the same as GPU-HC
bool CPU_HC_Solver::Read_RANSAC_Data( int tp_index ) {
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

//> This is the same as GPU-HC
void CPU_HC_Solver::Prepare_Target_Params( unsigned rand_seed_ ) {
    unsigned Edgel_Indx[3] = {0};
    const int offset_for_tangents = 18;
    std::array<int, 3> picked_samples;
    std::vector< std::array<int, 3> > target_params_match_indices;
#if FEED_RANDOM_SEED 
    srand (time(NULL));
#else
    srand (rand_seed_);
#endif

    for(int ti = 0; ti < NUM_OF_RANSAC_ITERATIONS; ti++) {
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
                (h_Target_Params + ti*(Num_Of_Params))[i*6 + j] = MAGMA_C_MAKE(h_Triplet_Edge_Locations(ei, j), 0.0);
            }
        }
        //> Tangents of the edgels
        for (int i = 0; i < 2; i++) {
            int ei = Edgel_Indx[i];
            for (int j = 0; j < 6; j++) {
                (h_Target_Params + ti*(Num_Of_Params))[i*6 + j + offset_for_tangents] = MAGMA_C_MAKE(h_Triplet_Edge_Tangents(ei, j), 0.0);
            }
        }
        (h_Target_Params + ti*(Num_Of_Params))[30] = MAGMA_C_MAKE(1.0, 0.0);
        (h_Target_Params + ti*(Num_Of_Params))[31] = MAGMA_C_MAKE(0.5, 0.0);
        (h_Target_Params + ti*(Num_Of_Params))[32] = MAGMA_C_MAKE(-1.0, 0.0);
    }
    // magma_cprint(pp->numOfVars+1, 1, (h_startSols + print_set_i*(pp->numOfTracks)*(pp->numOfVars+1) + print_sol_i * (pp->numOfVars+1)), (pp->numOfVars+1));
    // magma_cprint(Num_Of_Params, 1, h_Target_Params, Num_Of_Params);
    // for (int i = 0; i < Num_Of_Params; i++) {
    //     std::cout << std::setprecision(16) << MAGMA_C_REAL(h_Target_Params[i]) << "\t" << MAGMA_C_IMAG(h_Target_Params[i]) << std::endl;
    // }
}

void CPU_HC_Solver::Set_Initial_Array_Vals() {

    //> Set false to the status of all HC solutions
    for (int i = 0; i < Num_Of_Tracks*NUM_OF_RANSAC_ITERATIONS; i++) {
        h_is_Track_Converged[i]     = false;
        h_is_Track_Inf_Failed[i]    = false;
    }
}

void CPU_HC_Solver::Solve_by_CPU_HC() {
    std::cout << "CPU-HC computing ..." << std::endl << std::endl;

    if (HC_problem == "trifocal_2op1p_30x30") {
        CPU_HC_time = CPUHC_Generic_Solver( Num_Of_CPU_Cores, 
                                            &cpu_eval_dHdX_dHdt_trifocal_2op1p_30, 
                                            &cpu_eval_dHdX_H_trifocal_2op1p_30 );
    }

    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << "## Solving " << HC_print_problem_name << " by CPU-HC" << std::endl;

    //> Print out timings
    printf("## Timings:\n");
    printf(" - CPU Computation Time = %7.2f (ms)\n", (CPU_HC_time)*1000);

    //> Object for the Evaluations class
#if WRITE_CPUHC_CONVERGED_SOLS
    Evaluate_CPUHC_Sols->Write_Converged_Sols( h_CPU_HC_Track_Sols, h_is_Track_Converged );
#endif
    Evaluate_CPUHC_Sols->Evaluate_RANSAC_HC_Sols( h_CPU_HC_Track_Sols, h_is_Track_Converged, h_is_Track_Inf_Failed );

    //> Print out evaluation results
    std::cout << "\n## Evaluation of GPU-HC Solutions: "      << std::endl;
    std::cout << " - Number of Converged Solutions:       " << Evaluate_CPUHC_Sols->Num_Of_Coverged_Sols << std::endl;
    std::cout << " - Number of Real Solutions:            " << Evaluate_CPUHC_Sols->Num_Of_Real_Sols << std::endl;
    std::cout << " - Number of Infinity Failed Solutions: " << Evaluate_CPUHC_Sols->Num_Of_Inf_Sols << std::endl;
    
    Evaluate_CPUHC_Sols->Transform_GPUHC_Sols_to_Trifocal_Relative_Pose( h_CPU_HC_Track_Sols, h_is_Track_Converged, h_Camera_Intrinsic_Matrix );
    bool find_good_sol = Evaluate_CPUHC_Sols->get_Solution_with_Maximal_Support( Num_Of_Triplet_Edgels, h_Triplet_Edge_Locations, h_Triplet_Edge_Tangents, h_Camera_Intrinsic_Matrix );

    //> Reset the number of solutions
    Evaluate_CPUHC_Sols->Flush_Out_Data();
}

void CPU_HC_Solver::Free_Triplet_Edgels_Mem() {
    delete [] h_Triplet_Edge_Locations;
    delete [] h_Triplet_Edge_Tangents;
}

CPU_HC_Solver::~CPU_HC_Solver() {

    magma_free_cpu( h_Target_Params );
    magma_free_cpu( h_Track_Last_Success );
    magma_free_cpu( h_CPU_HC_Track_Sols );
    magma_free_cpu( h_Intermediate_Sols );
    magma_free_cpu( h_Start_Sols );
    magma_free_cpu( h_Start_Params );

    magma_free_cpu( h_cgesvA );
    magma_free_cpu( h_cgesvB );
    magma_free_cpu( ipiv );

    delete [] h_is_Track_Converged;
    delete [] h_is_Track_Inf_Failed;

    fflush( stdout );
    printf( "\n" );
    magma_finalize();
}

#endif
