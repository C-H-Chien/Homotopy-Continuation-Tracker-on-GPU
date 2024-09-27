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
#include <string>
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

class GPU_HC_Solver {
    
    magma_device_t cdev;                        // variable to indicate current gpu id
    magma_device_t cdevs;                       // variable to indicate gpu ids (for multiple GPUs)
    magma_queue_t  gpu_queues[MAX_NUM_OF_GPUS];  // magma queue variable, internally holds a cuda stream and a cublas handle
    magma_int_t    dev_id;

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
    bool                       *h_Found_Trifocal_Sols[MAX_NUM_OF_GPUS] = {NULL};    //> Indication of whether solutions are found for early aborting RANSAC process
    int                  *h_Trifocal_Sols_Batch_Index[MAX_NUM_OF_GPUS] = {NULL};    //> The solution index of which RANSAC is early aborted
    float                   *h_Camera_Intrinsic_Matrix;
    magmaFloatComplex       *h_Start_Sols;
    magmaFloatComplex       *h_Start_Params;
    int                     *h_dHdx_Index;
    int                     *h_dHdt_Index;
    int                     *h_unified_dHdx_dHdt_Index;
    
    float                   *h_Triplet_Edge_Locations;      //> in metrics
    float                   *h_Triplet_Edge_Tangents;       //> in metrics
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
    int                 *d_dHdx_Index[MAX_NUM_OF_GPUS] = {NULL};
    int                 *d_dHdt_Index[MAX_NUM_OF_GPUS] = {NULL};
    int    *d_unified_dHdx_dHdt_Index[MAX_NUM_OF_GPUS] = {NULL};
    magmaFloatComplex    **d_Start_Sols_array[MAX_NUM_OF_GPUS] = {NULL};
    float           *d_Triplet_Edge_Locations[MAX_NUM_OF_GPUS] = {NULL};
    float                 *d_Intrinsic_Matrix[MAX_NUM_OF_GPUS] = {NULL};
    bool               *d_Found_Trifocal_Sols[MAX_NUM_OF_GPUS] = {NULL};
    int          *d_Trifocal_Sols_Batch_Index[MAX_NUM_OF_GPUS] = {NULL};
    
public:
    bool                    Use_P2C;

    //> The timers
    real_Double_t           gpu_time[MAX_NUM_OF_GPUS] = {0.0};
    real_Double_t           transfer_h2d_time[MAX_NUM_OF_GPUS] = {0.0};
    real_Double_t           transfer_d2h_time[MAX_NUM_OF_GPUS] = {0.0};
    double                  multi_GPUs_time;

    //> Constructor
    GPU_HC_Solver() {};
    GPU_HC_Solver( YAML::Node );
    
    //> Member functions
    bool Read_Problem_Data();
    bool Read_RANSAC_Data( int tp_index );
    void Allocate_Arrays();
    void Prepare_Target_Params( unsigned rand_seed_ );
    void Data_Transfer_From_Host_To_Device();
    void Set_CUDA_Stream_Attributes();
    void Set_RANSAC_Abort_Arrays();
    void Solve_by_GPU_HC();
    void Export_Data();
    void Free_Triplet_Edgels_Mem();
    void Free_Arrays_for_Aborting_RANSAC();

    std::vector<unsigned>   Collect_Num_Of_Coverged_Sols;
    std::vector<unsigned>   Collect_Num_Of_Inf_Sols;
    std::vector<unsigned>   Collect_Num_Of_Real_Sols;

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

    //> GPU kernel settings from YAML file (merging all code implementation optimization)
    bool Use_Merge_Code_Optimization;

    //> Algorithmic settings from YAML file
    bool Truncate_HC_Path_by_Positive_Depths;
    bool Abort_RANSAC_by_Good_Sol;

    //> GPU atchitecture enabler
    bool GPU_arch_Ampere_and_above;

    //> Assign -1 to batch index of trifocal solutions. This is used for early aborting RANSAC purpose.
    void initialize_trifocal_sols_batch_index( int Num_Of_Batches, int* &h_Trifocal_Sols_Batch_Index ) {
        for (int i = 0; i < Num_Of_Batches; i++)
            h_Trifocal_Sols_Batch_Index[i] = -1;
    }

    void print_kernel_mode() {
        
        std::string str_code_opt     = (Use_Merge_Code_Optimization) ? "v" : "x";
        std::string str_trunc_path   = (Truncate_HC_Path_by_Positive_Depths) ? "v" : "x";
        std::string str_Abort_RANSAC = (Abort_RANSAC_by_Good_Sol) ? "v" : "x";

        if (Use_P2C) printf("PHC-(x) CodeOpt-(x) TrunPaths-(x) TrunRANSAC-(x)\n");
        else printf("PHC-(v) CodeOpt-(%s) TrunPaths-(%s) TrunRANSAC-(%s)\n", str_code_opt.c_str(), str_trunc_path.c_str(), str_Abort_RANSAC.c_str());
    }

    unsigned kernel_version;
    unsigned get_kernel_version_number() {
        if (Use_P2C) return 1;
        else if ( !Use_Merge_Code_Optimization && !Truncate_HC_Path_by_Positive_Depths && !Abort_RANSAC_by_Good_Sol ) return 2;
        else if (  Use_Merge_Code_Optimization && !Truncate_HC_Path_by_Positive_Depths && !Abort_RANSAC_by_Good_Sol ) return 3;
        else if (  Use_Merge_Code_Optimization &&  Truncate_HC_Path_by_Positive_Depths && !Abort_RANSAC_by_Good_Sol ) return 4;
        else if (  Use_Merge_Code_Optimization &&  Truncate_HC_Path_by_Positive_Depths &&  Abort_RANSAC_by_Good_Sol ) return 5;
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
    std::vector<int> candidate_batch_ids;

    //> Number of GPUs in use
    magma_int_t device_count = 0;
    int Num_Of_GPUs;

    void check_multiGPUs() {
        if (Num_Of_GPUs == 1) {
            std::string out_info = std::string("Only 1 GPU is used. Device ID = ");
            out_info.append( std::to_string(SET_GPU_DEVICE_ID) );
            LOG_INFO_MESG(out_info);
        }

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
