#ifndef magmaHC_kernels_h
#define magmaHC_kernels_h
// ============================================================================
// Header file declaring all kernels
//
// Modifications
//    Chiang-Heng Chien  22-10-31:   Initially Created (Copied from other repos)
//
//> (c) LEMS, Brown University
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// ============================================================================
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <chrono>
#include <memory>

#include "magma_v2.h"
#include "gpu-kernels/magmaHC-kernels.hpp"
#include "Data_Reader.hpp"
#include "Evaluations.hpp"
#include <yaml-cpp/yaml.h>

#define h_Triplet_Edge_Locations(i,j)    h_Triplet_Edge_Locations[(i) * 6 + (j)]
#define h_Triplet_Edge_Tangents(i,j)     h_Triplet_Edge_Tangents[(i) * 6 + (j)]

class Data_Reader;
class Evaluations;

template< typename T_index_mat >
class GPU_HC_Solver {
    
    magma_device_t cdev;                        // variable to indicate current gpu id
    magma_device_t cdevs;                       // variable to indicate gpu ids (for multiple GPUs)
    magma_queue_t  gpu_queues[MAX_NUM_OF_GPUS];  // magma queue variable, internally holds a cuda stream and a cublas handle

    //> Varaibles as sizes of arrays
    magma_int_t             dHdx_Index_Size;
    magma_int_t             dHdt_Index_Size;
    magma_int_t             unified_dHdx_dHdt_Index_Size;
    magma_int_t             dHdx_PHC_Coeffs_Size;
    magma_int_t             dHdt_PHC_Coeffs_Size;

    //> Variables and arrays on the CPU side
    magmaFloatComplex                *h_Homotopy_Sols[MAX_NUM_OF_GPUS] = {NULL};
    magmaFloatComplex                *h_Target_Params[MAX_NUM_OF_GPUS] = {NULL};
    magmaFloatComplex            *h_GPU_HC_Track_Sols[MAX_NUM_OF_GPUS] = {NULL};
    magmaFloatComplex                   *h_diffParams[MAX_NUM_OF_GPUS] = {NULL};
    magmaFloatComplex                *h_Debug_Purpose[MAX_NUM_OF_GPUS] = {NULL};
    bool                    *h_is_GPU_HC_Sol_Converge[MAX_NUM_OF_GPUS] = {NULL};
    bool                    *h_is_GPU_HC_Sol_Infinity[MAX_NUM_OF_GPUS] = {NULL};
    magmaFloatComplex              *h_dHdx_PHC_Coeffs[MAX_NUM_OF_GPUS] = {NULL};
    magmaFloatComplex              *h_dHdt_PHC_Coeffs[MAX_NUM_OF_GPUS] = {NULL};
    magmaFloatComplex       *h_Start_Sols;
    magmaFloatComplex       *h_Start_Params;
    T_index_mat             *h_dHdx_Index;
    T_index_mat             *h_dHdt_Index;
    T_index_mat             *h_unified_dHdx_dHdt_Index;
    float                   *h_Triplet_Edge_Locations;      //> in metrics
    float                   *h_Triplet_Edge_Tangents;       //> in metrics
    float                   h_Camera_Intrinsic_Matrix[9];
    float                   h_Camera_Pose21[12];
    float                   h_Camera_Pose31[12];

    magmaFloatComplex       *h_GPU_HC_Track_Sols_Stack;
    bool                    *h_is_GPU_HC_Sol_Converge_Stack;
    bool                    *h_is_GPU_HC_Sol_Infinity_Stack;

    //> Variables and arrays on the GPU side
    magmaFloatComplex_ptr     d_Homotopy_Sols[MAX_NUM_OF_GPUS] = {NULL};
    magmaFloatComplex_ptr     d_Target_Params[MAX_NUM_OF_GPUS] = {NULL};
    magmaFloatComplex           *d_diffParams[MAX_NUM_OF_GPUS] = {NULL};
    magmaFloatComplex **d_Homotopy_Sols_array[MAX_NUM_OF_GPUS] = {NULL};
    magmaFloatComplex        *d_Debug_Purpose[MAX_NUM_OF_GPUS] = {NULL};
    bool            *d_is_GPU_HC_Sol_Converge[MAX_NUM_OF_GPUS] = {NULL};
    bool            *d_is_GPU_HC_Sol_Infinity[MAX_NUM_OF_GPUS] = {NULL};
    magmaFloatComplex_ptr        d_Start_Sols[MAX_NUM_OF_GPUS] = {NULL};
    magmaFloatComplex_ptr      d_Start_Params[MAX_NUM_OF_GPUS] = {NULL};
    magmaFloatComplex_ptr   d_dHdx_PHC_Coeffs[MAX_NUM_OF_GPUS] = {NULL};
    magmaFloatComplex_ptr   d_dHdt_PHC_Coeffs[MAX_NUM_OF_GPUS] = {NULL};
    T_index_mat                 *d_dHdx_Index[MAX_NUM_OF_GPUS] = {NULL};
    T_index_mat                 *d_dHdt_Index[MAX_NUM_OF_GPUS] = {NULL};
    T_index_mat    *d_unified_dHdx_dHdt_Index[MAX_NUM_OF_GPUS] = {NULL};
    magmaFloatComplex    **d_Start_Sols_array[MAX_NUM_OF_GPUS] = {NULL};
    
public:
    bool                    Use_P2C;

    //> The timers
    real_Double_t           gpu_time[MAX_NUM_OF_GPUS] = {0.0};
    real_Double_t           gpu_max_time_from_multiple_GPUs = 0.0;
    real_Double_t           transfer_h2d_time[MAX_NUM_OF_GPUS] = {0.0};
    real_Double_t           transfer_d2h_time[MAX_NUM_OF_GPUS] = {0.0};
    double                  multi_GPUs_time;

    //> Constructor
    GPU_HC_Solver() {};
    GPU_HC_Solver( YAML::Node, int );
    
    //> Member functions
    bool Read_Problem_Data();
    bool Read_RANSAC_Data( int tp_index );
    void Allocate_Arrays();
    void Prepare_Target_Params( unsigned rand_seed_ );
    void Data_Transfer_From_Host_To_Device();
    void Set_CUDA_Stream_Attributes();
    void Solve_by_GPU_HC();
    void Export_Data();
    void Free_Triplet_Edgels_Mem();

    //> Destructor
    ~GPU_HC_Solver();

private:
    std::shared_ptr<Data_Reader> Load_Problem_Data = nullptr;
    std::shared_ptr<Evaluations> Evaluate_GPUHC_Sols = nullptr;
    YAML::Node Problem_Setting_YAML_File;
    
    std::string Problem_File_Path;
    std::string RANSAC_Data_File_Path;
    std::string Write_Files_Path;

    std::string HC_problem;
    std::string HC_print_problem_name;
    std::string GPUHC_type;
    int GPUHC_Max_Steps;
    int GPUHC_Max_Correction_Steps;
    int GPUHC_delta_t_incremental_steps;
    int Num_Of_Vars;
    int Num_Of_Params;
    int Num_Of_Tracks;
    int dHdx_Max_Terms;
    int dHdx_Max_Parts;
    int dHdt_Max_Terms;
    int dHdt_Max_Parts;
    int Max_Order_Of_T;

    double Num_Of_MBytes_Persistent_Data;
    double Num_Of_MBytes_Persistent_Cache;

    //> GPU kernel settings from YAML file
    bool Use_Runge_Kutta_in_a_Loop;
    int Data_Size_for_Indices;
    std::string Mem_for_Indices;
    bool Inline_Eval_Functions;
    bool Limit_Loop_Unroll;
    bool Use_L2_Persistent_Cache;

    //> Algorithmic settings from YAML file
    bool Truncate_HC_Path_by_Positive_Depths;

    void print_kernel_mode() {
        std::string str_32b         = (Data_Size_for_Indices == 32) ? "v" : "x";
        std::string str_GM          = (Mem_for_Indices == "GM") ? "v" : "x";
        std::string str_RKL         = (Use_Runge_Kutta_in_a_Loop) ? "v" : "x";
        std::string str_inline      = (Inline_Eval_Functions) ? "v" : "x";
        std::string str_lim_unroll  = (Limit_Loop_Unroll) ? "v" : "x";
        std::string str_trunc_path  = (Truncate_HC_Path_by_Positive_Depths) ? "v" : "x";
        std::string str_l2_cache    = (Use_L2_Persistent_Cache) ? "v" : "x";

        if (Use_P2C) printf("PHC-(x) RKL-(x) inline-(x) LimUnroll-(x) TrunPaths-(x)\n");
        else printf("PHC-(v) RKL-(%s) inline-(%s) LimUnroll-(%s) L2Cache-(%s) TrunPaths-(%s)\n", str_RKL.c_str(), str_inline.c_str(), str_lim_unroll.c_str(), str_l2_cache.c_str(), str_trunc_path.c_str());
    }

    unsigned kernel_version;
    unsigned get_kernel_version_number() {
        if (Use_P2C) return 1;
        else if ( !Use_Runge_Kutta_in_a_Loop && !Inline_Eval_Functions && !Limit_Loop_Unroll && !Use_L2_Persistent_Cache && !Truncate_HC_Path_by_Positive_Depths ) return 2;
        else if (  Use_Runge_Kutta_in_a_Loop && !Inline_Eval_Functions && !Limit_Loop_Unroll && !Use_L2_Persistent_Cache && !Truncate_HC_Path_by_Positive_Depths ) return 3;
        else if (  Use_Runge_Kutta_in_a_Loop &&  Inline_Eval_Functions && !Limit_Loop_Unroll && !Use_L2_Persistent_Cache && !Truncate_HC_Path_by_Positive_Depths ) return 4;
        else if (  Use_Runge_Kutta_in_a_Loop &&  Inline_Eval_Functions &&  Limit_Loop_Unroll && !Use_L2_Persistent_Cache && !Truncate_HC_Path_by_Positive_Depths ) return 5;
        else if (  Use_Runge_Kutta_in_a_Loop &&  Inline_Eval_Functions &&  Limit_Loop_Unroll &&  Use_L2_Persistent_Cache && !Truncate_HC_Path_by_Positive_Depths ) return 6;
        else if (  Use_Runge_Kutta_in_a_Loop &&  Inline_Eval_Functions &&  Limit_Loop_Unroll &&  Use_L2_Persistent_Cache &&  Truncate_HC_Path_by_Positive_Depths ) return 7;
        else if (  Use_Runge_Kutta_in_a_Loop && !Inline_Eval_Functions &&  Limit_Loop_Unroll && !Use_L2_Persistent_Cache && !Truncate_HC_Path_by_Positive_Depths ) return 8;
        else if (  Use_Runge_Kutta_in_a_Loop && !Inline_Eval_Functions &&  Limit_Loop_Unroll &&  Use_L2_Persistent_Cache && !Truncate_HC_Path_by_Positive_Depths ) return 9;
        else if (  Use_Runge_Kutta_in_a_Loop && !Inline_Eval_Functions &&  Limit_Loop_Unroll &&  Use_L2_Persistent_Cache &&  Truncate_HC_Path_by_Positive_Depths ) return 10;
        else {
            LOG_ERROR("Invalid configurations on the GPU settings / Algorithmic settings!");
            exit(1);
        }
    }

    //> RANSAC data
    std::string RANSAC_Dataset_Name;
    int Num_Of_Triplet_Edgels;
    int Num_Of_Coeffs_From_Params;
    std::vector<int> GPUHC_Actual_Sols_Steps_Collections;
    int sub_RANSAC_iters[MAX_NUM_OF_GPUS] = {0};
    int RANSAC_Workload_of_Last_GPU;

    //> Number of GPUs in use
    magma_int_t device_count = 0;
    int Num_Of_GPUs;
    // int *devices;
    // int Num_of_Devices_from_magma;
    // bool is_workload_even = true;

    void check_multiGPUs() {
        if( Num_Of_GPUs > MAX_NUM_OF_GPUS) {
            LOG_ERROR("Requested GPUs larger than MAX_NUM_OF_GPUS");
            printf("\033[1;31m[Requested GPUs] %d\t[Max GPUs] %d\033[0m\n", Num_Of_GPUs, MAX_NUM_OF_GPUS);
            exit(1);
        }

        if( Num_Of_GPUs > device_count ) {
            LOG_ERROR("Not enough GPUs");
            printf("\033[1;31m[Requested GPUs] %d\t[Available GPUs] %d\033[0m\n", Num_Of_GPUs, device_count);
            exit(1);
        }
    }
};

#endif
