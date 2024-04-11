#ifndef EVALUATIONS_H
#define EVALUATIONS_H
// ==========================================================================================================
// Header file declaring all kernels
//
// Modifications
//    Chiang-Heng Chien  23-02-05:   Initially created for evaluating the HC solutions, 
//                                   e.g., writing results to files, counting the converged solutions, etc.
//
//> (c) LEMS, Brown University
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// ==========================================================================================================
#include <stdlib.h>
#include <stdio.h>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <set>
#include <memory>

//> MAGMA
#include "magma_v2.h"

#include "definitions.hpp"
#include "util.hpp"

#define K(i,j)      K[(i) * 3 + (j)]

class util;

class Evaluations {
    
public:
    //> Constructor
    Evaluations( std::string, int, int );

    //> Destructor
    ~Evaluations();

    //> Write data to files
    void Write_Converged_Sols( magmaFloatComplex *h_GPU_HC_Track_Sols, bool *h_is_GPU_HC_Sol_Converge );

    //> Evaluate GPU-HC Solutions
    void Evaluate_GPUHC_Sols( magmaFloatComplex *h_GPU_HC_Track_Sols, bool *h_is_GPU_HC_Sol_Converge, bool *h_is_GPU_HC_Sol_Infinity, int ransac_sample_offset );
    void Evaluate_RANSAC_GPUHC_Sols( magmaFloatComplex *h_GPU_HC_Track_Sols, bool *h_is_GPU_HC_Sol_Converge, bool *h_is_GPU_HC_Sol_Infinity );
    void Find_Unique_Sols( magmaFloatComplex *h_GPU_HC_Track_Sols, bool *h_is_GPU_HC_Sol_Converge );

    //> Evaluate RANSAC solutions with the ground-truths
    void Transform_GPUHC_Sols_to_Trifocal_Relative_Pose( magmaFloatComplex *h_GPU_HC_Track_Sols, bool *h_is_GPU_HC_Sol_Converge );
    float get_Rotation_Residual( float* GT_R, std::array<float, 9> Sol_R );
    float get_Translation_Residual( float* GT_Transl, std::array<float, 3> Sol_Transl );
    void Measure_Relative_Pose_Error( float GT_Pose21[12], float GT_Pose31[12] );
    
    //> Some evaluation data
    unsigned Num_Of_Coverged_Sols;
    unsigned Num_Of_Inf_Sols;
    unsigned Num_Of_Real_Sols;
    unsigned Num_Of_Unique_Sols;

    //> Some RANSAC solutions evaluation data
    float Percentage_Of_Convergence;
    float Percentage_Of_Inf_Sols;
    float Percentage_Of_Real_Sols;
    float Percentage_Of_Unique_Sols;
    float Min_Residual_R21;
    float Min_Residual_R31;
    float Min_Residual_t21;
    float Min_Residual_t31;
    bool success_flag;

private:
    //> util
    std::shared_ptr<util> MVG_Utility = nullptr;

    //> output streams for files to be written
    std::ofstream GPUHC_Track_Sols_File;

    float Var_Diff_In_Real_Part( magmaFloatComplex *Sol1, magmaFloatComplex *Sol2, int sol1_offset, int sol2_offset, int var_index ) {
        return std::fabs(MAGMA_C_REAL((Sol1 + sol1_offset*(num_of_variables+1))[var_index]) - MAGMA_C_REAL((Sol2 + sol2_offset*(num_of_variables+1))[var_index]));
    };
    float Var_Diff_In_Imag_Part( magmaFloatComplex *Sol1, magmaFloatComplex *Sol2, int sol1_offset, int sol2_offset, int var_index ) {
        return std::fabs(MAGMA_C_IMAG((Sol1 + sol1_offset*(num_of_variables+1))[var_index]) - MAGMA_C_IMAG((Sol2 + sol2_offset*(num_of_variables+1))[var_index]));
    };

    std::vector< int > Unique_Sols_Index;

    std::string WRITE_FILES_PATH;

    const int num_of_tracks;
    const int num_of_variables;

    //> RANSAC data
    float K[9];

    //> Individual relative pose
    float *Rot21;
    float *Rot31;
    float *Transl21;
    float *Transl31;
    float* Sol_Rotm_21;
    float* Sol_Rotm_31;
    //> Collecting all relative poses from RANSAC
    std::vector<std::array<float, 3>> normalized_t21s;
    std::vector<std::array<float, 3>> normalized_t31s;
    std::vector<std::array<float, 9>> normalized_R21s;
    std::vector<std::array<float, 9>> normalized_R31s;
    // std::vector<float *> F21s;
    // std::vector<float *> F31s;

    std::array<float, 3> normalized_t21;
    std::array<float, 3> normalized_t31;
    std::array<float, 9> normalized_R21;
    std::array<float, 9> normalized_R31;

    float* R_gt_R;
    float* Sols_R_;
    float* Sol_Transl_;
};

#endif
