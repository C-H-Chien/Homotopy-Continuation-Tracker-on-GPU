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
    
    magma_device_t cdev;       // variable to indicate current gpu id
    magma_queue_t my_queue;    // magma queue variable, internally holds a cuda stream and a cublas handle

    //> Varaibles as sizes of arrays
    magma_int_t             dHdx_Index_Size;
    magma_int_t             dHdt_Index_Size;
    magma_int_t             dHdx_PHC_Coeffs_Size;
    magma_int_t             dHdt_PHC_Coeffs_Size;
    magma_int_t             ldd_phc_Params_Hx;
    magma_int_t             ldd_phc_Params_Ht;

    //> Variables and arrays on the CPU side
    magmaFloatComplex       *h_GPU_HC_Track_Sols;
    magmaFloatComplex       *h_diffParams;
    magmaFloatComplex       *h_Start_Sols;
    magmaFloatComplex       *h_Homotopy_Sols;
    magmaFloatComplex       *h_Start_Params;
    magmaFloatComplex       *h_Target_Params;
    magmaFloatComplex       *h_dHdx_PHC_Coeffs;
    magmaFloatComplex       *h_dHdt_PHC_Coeffs;
    T_index_mat             *h_dHdx_Index;
    T_index_mat             *h_dHdt_Index;
    float                   *h_Triplet_Edge_Locations;      //> in metrics
    float                   *h_Triplet_Edge_Tangents;       //> in metrics
    magmaFloatComplex       *h_Debug_Purpose;
    bool                    *h_is_GPU_HC_Sol_Converge;
    bool                    *h_is_GPU_HC_Sol_Infinity;
    float                   h_Camera_Intrinsic_Matrix[9];
    float                   h_Camera_Pose21[12];
    float                   h_Camera_Pose31[12];

    //> Variables and arrays on the GPU side
    magmaFloatComplex_ptr   d_Start_Sols, d_Homotopy_Sols;
    magmaFloatComplex_ptr   d_Start_Params, d_Target_Params;
    magmaFloatComplex_ptr   d_dHdx_PHC_Coeffs;
    magmaFloatComplex_ptr   d_dHdt_PHC_Coeffs;
    T_index_mat             *d_dHdx_Index;
    T_index_mat             *d_dHdt_Index;
    magmaFloatComplex       **d_Start_Sols_array;
    magmaFloatComplex       **d_Homotopy_Sols_array;
    magmaFloatComplex       *d_diffParams;
    magmaFloatComplex       *d_Debug_Purpose;
    bool                    *d_is_GPU_HC_Sol_Converge;
    bool                    *d_is_GPU_HC_Sol_Infinity;

public:
    bool                    Use_P2C;

    //> The timers
    real_Double_t           gpu_time;
    real_Double_t           transfer_h2d_time;
    real_Double_t           transfer_d2h_time;

    //> Constructor
    GPU_HC_Solver() {};
    GPU_HC_Solver( YAML::Node, int );
    
    //> Member functions
    bool Read_Problem_Data();
    bool Read_RANSAC_Data( int tp_index );
    void Allocate_Arrays();
    void Prepare_Target_Params( unsigned rand_seed_ );
    void Data_Transfer_From_Host_To_Device();
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

    //> GPU kernel settings from YAML file
    bool Use_Runge_Kutta_in_a_Loop;
    int Data_Size_for_Indices;
    std::string Mem_for_Indices;
    bool Inline_Eval_Functions;
    bool Limit_Loop_Unroll;

    //> Algorithmic settings from YAML file
    bool Truncate_HC_Path_by_Positive_Depths;

    //> RANSAC data
    std::string RANSAC_Dataset_Name;
    int Num_Of_Triplet_Edgels;
    int Num_Of_Coeffs_From_Params;
    std::vector<int> GPUHC_Actual_Sols_Steps_Collections;
};

#endif
