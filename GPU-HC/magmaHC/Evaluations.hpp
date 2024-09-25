#ifndef EVALUATIONS_H
#define EVALUATIONS_H
// ==========================================================================================================
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

#define h_Triplet_Edge_Locations(i,j)    h_Triplet_Edge_Locations[(i) * 6 + (j)]
#define h_Triplet_Edge_Tangents(i,j)     h_Triplet_Edge_Tangents[(i) * 6 + (j)]

class util;

class Evaluations {
    
public:
    //> Constructor
    Evaluations( std::string, std::string, int, int );

    //> Destructor
    ~Evaluations();

    //> Write data to files
    void Write_Converged_Sols( magmaFloatComplex *h_HC_Track_Sols, bool *h_is_HC_Sol_Converge );

    //> Evaluate GPU-HC Solutions
    void Evaluate_HC_Sols( magmaFloatComplex *h_HC_Track_Sols, bool *h_is_HC_Sol_Converge, bool *h_is_HC_Sol_Infinity, int ransac_sample_offset );
    void Evaluate_RANSAC_HC_Sols( magmaFloatComplex *h_HC_Track_Sols, bool *h_is_HC_Sol_Converge, bool *h_is_HC_Sol_Infinity );
    void Find_Unique_Sols( magmaFloatComplex *h_GPU_HC_Track_Sols, bool *h_is_GPU_HC_Sol_Converge );

    //> Processing rotations and translations from the GPU-HC solutions of trifocal pose estimation problem (30x30)
    void Convert_Trifocal_Translation( magmaFloatComplex *h_GPU_HC_Track_Sols );
    void Convert_Trifocal_Rotation( magmaFloatComplex *h_GPU_HC_Track_Sols );

    //> Evaluate RANSAC solutions with the ground-truths
    void Transform_GPUHC_Sols_to_Trifocal_Relative_Pose( magmaFloatComplex *h_GPU_HC_Track_Sols, bool *h_is_GPU_HC_Sol_Converge, float *IntrinsicMatrix );
    float get_Rotation_Residual( float* GT_R, std::array<float, 9> Sol_R );
    float get_Translation_Residual( float* GT_Transl, std::array<float, 3> Sol_Transl );
    void get_HC_Steps_of_Actual_Sols( magmaFloatComplex *h_Debug_Purpose );
    void Measure_Relative_Pose_Error_from_All_Real_Sols( float GT_Pose21[12], float GT_Pose31[12], magmaFloatComplex *h_Debug_Purpose );
    void Measure_Relative_Pose_Error( float GT_Pose21[12], float GT_Pose31[12] );
    bool get_Solution_with_Maximal_Support( unsigned Num_Of_Triplet_Edgels, float* h_Triplet_Edge_Locations, float* h_Triplet_Edge_Tangents, float *K );

    void Check_Deviations_of_Veridical_Sol_from_GT( magmaFloatComplex *h_GPU_HC_Track_Sols, float GT_Pose21[12], float GT_Pose31[12] );

    //> Others
    void Flush_Out_Data();
    void Write_HC_Steps_of_Actual_Solutions( std::vector<int> GPUHC_Actual_Sols_Steps_Collections ) {
        for (int i = 0; i < GPUHC_Actual_Sols_Steps_Collections.size(); i++) 
            HC_Actual_Sols_Steps_File << GPUHC_Actual_Sols_Steps_Collections[i] << "\n";
    };
    
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
    std::vector<int> HC_steps_of_actual_solutions;
    unsigned Max_Num_Of_Reproj_Inliers_Views21;
    unsigned Max_Num_Of_Reproj_Inliers_Views31;
    std::array<float, 9> R21_w_Max_Supports;
    std::array<float, 9> R31_w_Max_Supports;
    std::array<float, 3> t21_w_Max_Supports;
    std::array<float, 3> t31_w_Max_Supports;
    std::vector<int> Max_Reproj_Inliers_Support_Views21_Index;
    std::vector<int> Max_Reproj_Inliers_Support_Views31_Index;

private:
    //> util
    std::shared_ptr<util> MVG_Utility = nullptr;

    std::string evaluate_GPUHC_or_CPUHC;

    //> output streams for files to be written
    std::ofstream HC_Track_Sols_File;
    std::ofstream HC_Actual_Sols_Steps_File;

    float Var_Diff_In_Real_Part( magmaFloatComplex *Sol1, magmaFloatComplex *Sol2, int sol1_offset, int sol2_offset, int var_index ) {
        return std::fabs(MAGMA_C_REAL((Sol1 + sol1_offset*(num_of_variables+1))[var_index]) - MAGMA_C_REAL((Sol2 + sol2_offset*(num_of_variables+1))[var_index]));
    };
    float Var_Diff_In_Imag_Part( magmaFloatComplex *Sol1, magmaFloatComplex *Sol2, int sol1_offset, int sol2_offset, int var_index ) {
        return std::fabs(MAGMA_C_IMAG((Sol1 + sol1_offset*(num_of_variables+1))[var_index]) - MAGMA_C_IMAG((Sol2 + sol2_offset*(num_of_variables+1))[var_index]));
    };

    void get_GT_Rotation( float GT_Pose[12], float* &GT_Rot )       { std::copy( GT_Pose,      GT_Pose + 9,  GT_Rot   ); }
    void get_GT_Translation( float GT_Pose[12], float* &GT_Transl ) { std::copy( GT_Pose + 9, GT_Pose + 12, GT_Transl ); }

    // template< typename T, int n >
    // void copy_array_values( std::array<T, n> Copy_From, T* &Copy_To ) { for(int i = 0; i < n; i++) Copy_To[i] = Copy_From[i]; }

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
    std::vector<std::array<float, 9>> F21s;
    std::vector<std::array<float, 9>> F31s;

    std::array<float, 3> normalized_t21;
    std::array<float, 3> normalized_t31;
    std::array<float, 9> normalized_R21;
    std::array<float, 9> normalized_R31;
    std::array<float, 9> FundMatrix21;
    std::array<float, 9> FundMatrix31;

    float* R_gt_R;
    float* Sols_R_;
    float* Sol_Transl_;

    float* GT_Rot21;
    float* GT_Rot31;
    float* GT_Transl21;
    float* GT_Transl31;

    std::vector<int> real_track_indices;

    unsigned Num_Of_Reproj_Err_Inliers_Views21;
    unsigned Num_Of_Reproj_Err_Inliers_Views31;
};

#endif
