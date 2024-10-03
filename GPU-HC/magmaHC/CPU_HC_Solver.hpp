#ifndef CPU_HC_SOLVER_HPP
#define CPU_HC_SOLVER_HPP
// ============================================================================
// CPU-HC solver class
//
// Modifications
//    Chiang-Heng Chien  24-07-09:  Copied and edited from CPU_HC_Solver.hpp
//
//> (c) LEMS, Brown University
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// ============================================================================
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <functional>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <chrono>
#include <memory>

#include "magma_v2.h"
#include "magma_lapack.h"
#include "Data_Reader.hpp"
#include "Evaluations.hpp"
#include <yaml-cpp/yaml.h>

//> All evaluation functions (TODO: use bash script to automatically add all header files from the cpu-jacobian-evals directory)
#include "./cpu-jacobian-evals/cpu-eval-dHdX_dHdt_trifocal_2op1p_30x30.hpp"
#include "./cpu-jacobian-evals/cpu-eval-dHdX_H_trifocal_2op1p_30x30.hpp"

#include "./cpu-jacobian-evals/cpu-eval-indx_trifocal_2op1p_30x30.hpp"

#define h_Triplet_Edge_Locations(i,j)    h_Triplet_Edge_Locations[(i) * 6 + (j)]
#define h_Triplet_Edge_Tangents(i,j)     h_Triplet_Edge_Tangents[(i) * 6 + (j)]

//> Evaluation functions for CPU-HC
using Eval_dHdX_dHdt = std::function<void(int s, float t, int N, magmaFloatComplex* track, magmaFloatComplex* start_params, magmaFloatComplex* target_params, magmaFloatComplex* cgesvA, magmaFloatComplex* cgesvB)>;
using Eval_dHdX_H    = std::function<void(int s, float t, int N, magmaFloatComplex* track, magmaFloatComplex* start_params, magmaFloatComplex* target_params, magmaFloatComplex* cgesvA, magmaFloatComplex* cgesvB)>;

//> Evaluation method consistent with GPU-HC
using Eval_dHdX = std::function<void(const int Num_Of_Vars, const int dHdx_Max_Terms, const int dHdx_Entry_Offset, const int dHdx_Max_Parts, const int* dHdx_Index, magmaFloatComplex* variable, magmaFloatComplex* param_homotopy, magmaFloatComplex* cgesvA)>;
using Eval_dHdt = std::function<void(const int Num_Of_Vars, const int dHdt_Max_Terms, const int dHdt_Max_Parts,    const int* h_dHdt_Index,  magmaFloatComplex* variable, magmaFloatComplex* param_homotopy,  magmaFloatComplex* cgesvB,   magmaFloatComplex* diff_params)>;
using Eval_H    = std::function<void(const int Num_Of_Vars, const int dHdt_Max_Terms, const int dHdt_Max_Parts,    const int* h_dHdt_Index,  magmaFloatComplex* variable, magmaFloatComplex* param_homotopy,  magmaFloatComplex* cgesvB)>;

class Data_Reader;
class Evaluations;

class CPU_HC_Solver {

    //> Variable as sizes of arrays
    magma_int_t             dHdx_Index_Size;
    magma_int_t             dHdt_Index_Size;
    
    //> Variables and arrays
    magmaFloatComplex       *h_Target_Params;
    magmaFloatComplex       *h_CPU_HC_Track_Sols;
    magmaFloatComplex       *h_Track_Last_Success;
    magmaFloatComplex       *h_Intermediate_Sols;
    magmaFloatComplex       *h_Start_Sols;
    magmaFloatComplex       *h_Start_Params;
    magmaFloatComplex       *h_cgesvA;  //> Matrix A in a linear system Ax=b
    magmaFloatComplex       *h_cgesvB;  //> Vector b in a linear system Ax=b
    magmaFloatComplex       *h_diff_params;
    magma_int_t             *ipiv;
    int                     *h_dHdx_Index;
    int                     *h_dHdt_Index;
    bool                    *h_is_Track_Converged;
    bool                    *h_is_Track_Inf_Failed;

    float                   *h_Triplet_Edge_Locations;      //> in metrics
    float                   *h_Triplet_Edge_Tangents;       //> in metrics
    float                   h_Camera_Intrinsic_Matrix[9];
    float                   h_Camera_Pose21[12];
    float                   h_Camera_Pose31[12];
    
public:

    //> The timers
    double                  CPU_HC_time;

    //> Constructor
    CPU_HC_Solver() {};
    CPU_HC_Solver( YAML::Node );
    
    //> Member functions
    bool Read_Problem_Data();
    bool Read_RANSAC_Data( int tp_index );
    void Allocate_Arrays();
    void Prepare_Target_Params( unsigned rand_seed_ );
    void Set_Initial_Array_Vals();
    void Solve_by_CPU_HC();
    void Free_Triplet_Edgels_Mem();
    real_Double_t CPUHC_Generic_Solver(const int num_of_cores, const Eval_dHdX_dHdt& Eval_dHdX_dHdt_func, const Eval_dHdX_H& Eval_dHdX_H_func );
    real_Double_t CPUHC_Generic_Solver_Eval_by_Indx( const int num_of_cores, const Eval_dHdX& Eval_dHdX_func, const Eval_dHdt& Eval_dHdt_func, const Eval_H& Eval_H_func );
    void CPUHC_get_Runge_Kutta_x_for_k2(magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_Track_cpu, magmaFloatComplex *h_cgesvB, float one_half_delta_t, float delta_t, float &t0 );
    void CPUHC_get_Runge_Kutta_x_for_k3(magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_Track_cpu, magmaFloatComplex *h_cgesvB, magmaFloatComplex *h_Track_last_success, float one_half_delta_t, float delta_t );
    void CPUHC_get_Runge_Kutta_x_for_k4(magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_Track_cpu, magmaFloatComplex *h_cgesvB, magmaFloatComplex *h_Track_last_success, float one_half_delta_t, float &t0, float delta_t );
    void CPUHC_make_Runge_Kutta_prediction(magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_Track_cpu, magmaFloatComplex *h_cgesvB, float delta_t );
    void CPUHC_make_correction( magmaFloatComplex *h_CPU_HC_Track_Sols, magmaFloatComplex *h_cgesvB, bool &is_successful, bool &is_inf_failed );

    void cpu_eval_compute_param_homotopy( float t, magmaFloatComplex* h_param_homotopy, magmaFloatComplex* start_params, magmaFloatComplex* target_params ) {
        for (int i = 0; i <= Num_Of_Params; i++) {
            h_param_homotopy[ i ] = target_params[ i ] * t + start_params[ i ] * (1.0 - t);
        }
    }

    //> Evaluations
    // void cpu_eval_indx_dHdX_trifocal_2op1p_30( magmaFloatComplex* variable, magmaFloatComplex* param_homotopy, magmaFloatComplex* cgesvA, const int* dHdx_indices );
    // void cpu_eval_indx_dHdt_trifocal_2op1p_30( magmaFloatComplex* variable, magmaFloatComplex* param_homotopy, magmaFloatComplex* cgesvB, magmaFloatComplex* diff_params, const int* dHdt_indices );
    // void cpu_eval_indx_H_trifocal_2op1p_30( magmaFloatComplex* variables, magmaFloatComplex* dHdt_indices, magmaFloatComplex* param_homotopy, magmaFloatComplex* cgesvB );

    std::vector<unsigned>   Collect_Num_Of_Coverged_Sols;
    std::vector<unsigned>   Collect_Num_Of_Inf_Sols;
    std::vector<unsigned>   Collect_Num_Of_Real_Sols;

    //> Destructor
    ~CPU_HC_Solver();

private:
    //> OpenMP number of threads
    magma_int_t nthreads;

    std::shared_ptr<Data_Reader> Load_Problem_Data = nullptr;
    std::shared_ptr<Evaluations> Evaluate_CPUHC_Sols = nullptr;
    YAML::Node Problem_Setting_YAML_File;
    
    std::string Problem_File_Path;
    std::string RANSAC_Data_File_Path;
    std::string Write_Files_Path;

    std::string HC_problem;
    std::string HC_print_problem_name;
    int CPUHC_Max_Steps;
    int CPUHC_Max_Correction_Steps;
    int CPUHC_delta_t_incremental_steps;
    int Num_Of_Vars;
    int Num_Of_Params;
    int Num_Of_Tracks;
    int Num_Of_CPU_Cores;
    int RANSAC_Sol_Offset;
    int RANSAC_Sol_Offset_with_Dummy_Variable;
    int dHdx_Max_Terms;
    int dHdx_Max_Parts;
    int dHdt_Max_Terms;
    int dHdt_Max_Parts;
    int dHdx_Entry_Offset;

    int Num_Of_Inf_Failed_Sols;
    int Num_Of_Successful_Sols;

    //> Algorithmic settings from YAML file
    bool Truncate_HC_Path_by_Positive_Depths;

    //> RANSAC data
    std::string RANSAC_Dataset_Name;
    int Num_Of_Triplet_Edgels;
};

#endif
