#ifndef EVALUATIONS_CPP
#define EVALUATIONS_CPP
// =============================================================================================================================
//
// ChangLogs
//    24-02-05:   Initially created for definitions of functions in the ``Evaluations" class.
//    24-03-26:   Replace macro definitions of number of tracks and number of variables by passing them as variables
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
#include <math.h>

#include "definitions.hpp"
#include "Evaluations.hpp"

//> Constructor
Evaluations::Evaluations( std::string Output_Files_Path, std::string GPU_or_CPU, int num_of_tracks, int num_of_vars )
  : WRITE_FILES_PATH(Output_Files_Path), evaluate_GPUHC_or_CPUHC(GPU_or_CPU), num_of_tracks(num_of_tracks), num_of_variables(num_of_vars)
{
  //> Initialize to zeros
  Num_Of_Inf_Sols      = 0;
  Num_Of_Coverged_Sols = 0;
  Num_Of_Real_Sols     = 0;
  Num_Of_Unique_Sols   = 0;

  Percentage_Of_Convergence = 0.0;
  Percentage_Of_Inf_Sols = 0.0;
  Percentage_Of_Real_Sols = 0.0;
  Percentage_Of_Unique_Sols = 0.0;
  success_flag = false;
  Min_Residual_R21 = 100.0;
  Min_Residual_R31 = 100.0;
  Min_Residual_t21 = 100.0;
  Min_Residual_t31 = 100.0;

  //> Write successful HC track solutions to files
  std::string sols_file_name;
  std::string hc_steps_of_actual_sols_file_name;
  if (evaluate_GPUHC_or_CPUHC == "GPU-HC") {
    sols_file_name = "GPU_Converged_HC_tracks.txt";
    hc_steps_of_actual_sols_file_name = "GPUHC_Steps_of_Actual_Solutions.txt";
  }
  else if (evaluate_GPUHC_or_CPUHC == "CPU-HC") {
    sols_file_name = "CPU_Converged_HC_tracks.txt";
    hc_steps_of_actual_sols_file_name = "CPUHC_Steps_of_Actual_Solutions.txt";
  }
  else {
    LOG_ERROR("Invalid GPU_or_CPU input parameter for the Evaluation constructor.");
  }

  std::string write_sols_file_dir = WRITE_FILES_PATH + sols_file_name;
  HC_Track_Sols_File.open(write_sols_file_dir);
  if ( !HC_Track_Sols_File.is_open() ) LOG_FILE_ERROR(write_sols_file_dir);

  //> Write successful HC track solutions to files
  std::string write_actual_sols_HC_steps_file_dir = WRITE_FILES_PATH + hc_steps_of_actual_sols_file_name;
  HC_Actual_Sols_Steps_File.open(write_actual_sols_HC_steps_file_dir);
  if ( !HC_Actual_Sols_Steps_File.is_open() ) LOG_FILE_ERROR(write_actual_sols_HC_steps_file_dir);

  //> util class
  MVG_Utility = std::shared_ptr<util>(new util());

  //> Allocate relative pose arrays
  Rot21     = new float[9];
  Rot31     = new float[9];
  Transl21  = new float[3];
  Transl31  = new float[3];

  GT_Rot21     = new float[9];
  GT_Rot31     = new float[9];
  GT_Transl21  = new float[3];
  GT_Transl31  = new float[3];

  Sol_Rotm_21 = new float[9];
  Sol_Rotm_31 = new float[9];

  Sol_Transl_ = new float[3];
  R_gt_R      = new float[9];
  Sols_R_     = new float[9];
}

void Evaluations::Flush_Out_Data() {
  Num_Of_Inf_Sols = 0;
  Num_Of_Coverged_Sols = 0;
  Num_Of_Real_Sols = 0;
  Num_Of_Unique_Sols = 0;
  Percentage_Of_Convergence = 0.0;
  Percentage_Of_Inf_Sols = 0.0;
  Percentage_Of_Real_Sols = 0.0;
  Percentage_Of_Unique_Sols = 0.0;
  success_flag = false;
  Min_Residual_R21 = 100.0;
  Min_Residual_R31 = 100.0;
  Min_Residual_t21 = 100.0;
  Min_Residual_t31 = 100.0;
  real_track_indices.clear();
  HC_steps_of_actual_solutions.clear();

  normalized_t21s.clear();
  normalized_t31s.clear();
  normalized_R21s.clear();
  normalized_R31s.clear();
  F21s.clear();
  F31s.clear();

  Max_Reproj_Inliers_Support_Views21_Index.clear();
  Max_Reproj_Inliers_Support_Views31_Index.clear();
}

void Evaluations::Write_Converged_Sols( \
    magmaFloatComplex *h_HC_Track_Sols, \
    bool *h_is_HC_Sol_Converge ) 
{
  LOG_INFO_MESG("Writing HC converged solutions to a file ...");
  int counter = 0;
  for (int ri = 0; ri < NUM_OF_RANSAC_ITERATIONS; ri++) {
    HC_Track_Sols_File << "-------------------- RANSAC Iteration " << ri+1 << " --------------------\n\n";
    for (int bs = 0; bs < num_of_tracks; bs++) {
      HC_Track_Sols_File << std::setprecision(10);

      if ((h_is_HC_Sol_Converge + ri * num_of_tracks)[ bs ] == 1) {
        HC_Track_Sols_File << counter << "\n";
        for (int vs = 0; vs < num_of_variables; vs++) {
          HC_Track_Sols_File << std::setprecision(20) << MAGMA_C_REAL((h_HC_Track_Sols + ri * num_of_tracks * (num_of_variables+1) + bs * (num_of_variables+1))[vs]) << "\t" \
                             << std::setprecision(20) << MAGMA_C_IMAG((h_HC_Track_Sols + ri * num_of_tracks * (num_of_variables+1) + bs * (num_of_variables+1))[vs]) << "\n";
        }
        HC_Track_Sols_File << "\n";
      }
      counter++;
    }
    HC_Track_Sols_File << "\n";
  }
}

void Evaluations::Evaluate_HC_Sols( \
    magmaFloatComplex *h_HC_Track_Sols, \
    bool *h_is_HC_Sol_Converge, \
    bool *h_is_HC_Sol_Infinity, \
    int ransac_sample_offset ) 
{
  //> Count the number of converged solutions, the number of infinity failed solutions, and the number of real solutions
  for (int bs = 0; bs < num_of_tracks; bs++) {
    if ( (h_is_HC_Sol_Converge + num_of_tracks * ransac_sample_offset)[ bs ] ) Num_Of_Coverged_Sols++;
    if ( (h_is_HC_Sol_Infinity + num_of_tracks * ransac_sample_offset)[ bs ] ) Num_Of_Inf_Sols++;

    int Num_Of_Real_Vars = 0;
    if ((h_is_HC_Sol_Converge + num_of_tracks * ransac_sample_offset)[ bs ] == 1) {
      for (int vs = 0; vs < num_of_variables; vs++) {
        if (fabs(MAGMA_C_IMAG((h_HC_Track_Sols + num_of_tracks * (num_of_variables+1) * ransac_sample_offset + bs * (num_of_variables+1))[vs])) <= ZERO_IMAG_PART_TOL_FOR_SP) {
            Num_Of_Real_Vars++;
        }
      }
    }

    if (Num_Of_Real_Vars == num_of_variables) Num_Of_Real_Sols++;
  }
}

void Evaluations::Evaluate_RANSAC_HC_Sols( \
    magmaFloatComplex *h_HC_Track_Sols, \
    bool *h_is_HC_Sol_Converge, \
    bool *h_is_HC_Sol_Infinity )
{
  //> Loop over all RANSAC iterations
  for (int ri = 0; ri < NUM_OF_RANSAC_ITERATIONS; ri++) {
    Evaluate_HC_Sols( h_HC_Track_Sols, h_is_HC_Sol_Converge, h_is_HC_Sol_Infinity, ri );
  }

  Percentage_Of_Convergence = (float)Num_Of_Coverged_Sols / (float)(num_of_tracks * NUM_OF_RANSAC_ITERATIONS);
  Percentage_Of_Inf_Sols    = (float)Num_Of_Inf_Sols / (float)(num_of_tracks * NUM_OF_RANSAC_ITERATIONS);
  Percentage_Of_Real_Sols   = (float)Num_Of_Real_Sols / (float)(num_of_tracks * NUM_OF_RANSAC_ITERATIONS);
}

void Evaluations::Find_Unique_Sols( magmaFloatComplex *h_GPU_HC_Track_Sols, bool *h_is_GPU_HC_Sol_Converge ) {

  std::set< int > Duplicate_Sol_Index;
  std::set< int > Skip_Sol_Index;
  for (int bs = 0; bs < num_of_tracks; bs++) {

    //> Make sure that the solution is converged before evaluating its uniqueness
    if ( h_is_GPU_HC_Sol_Converge[ bs ] ) {
      
      if (!Skip_Sol_Index.empty()) {

        //if(std::find(Skip_Sol_Index.begin(), Skip_Sol_Index.end(), bs) != Skip_Sol_Index.end())
        if (Skip_Sol_Index.find(bs) != Skip_Sol_Index.end())
          continue;

        //> Flush out the set
        Duplicate_Sol_Index.clear();
      }

      for (int ds = bs+1; ds < num_of_tracks; ds++) {

        int Num_Of_Duplicate_Vars = 0;
        for (int vs = 0; vs < num_of_variables; vs++) {

          //> If both the real and imaginery parts are very close
          if ( Var_Diff_In_Real_Part(h_GPU_HC_Track_Sols, h_GPU_HC_Track_Sols, bs, ds, vs) < DUPLICATE_SOL_DIFF_TOL && \
               Var_Diff_In_Imag_Part(h_GPU_HC_Track_Sols, h_GPU_HC_Track_Sols, bs, ds, vs) < DUPLICATE_SOL_DIFF_TOL ) {
              Num_Of_Duplicate_Vars++;
          }
        }
        if (Num_Of_Duplicate_Vars == num_of_variables) Duplicate_Sol_Index.insert(ds);
      }

      //> If the duplicate solution index vector is empty
      if (Duplicate_Sol_Index.empty()) {
        Num_Of_Unique_Sols++;
        Unique_Sols_Index.push_back(bs);
      }
      else Skip_Sol_Index = Duplicate_Sol_Index;
    }
  }

#if DEBUG_EVALUATOR
  std::cout << "Indices of unique solutions: ";
  for (int i = 0; i < Unique_Sols_Index.size(); i++) {
    std::cout << Unique_Sols_Index[i] << "\t";
  }
  std::cout << std::endl;
#endif
}

void Evaluations::Convert_Trifocal_Translation( magmaFloatComplex *h_GPU_HC_Track_Sols ) {
  //> \transl_{21}
  Transl21[0] = MAGMA_C_REAL(h_GPU_HC_Track_Sols[18]);
  Transl21[1] = MAGMA_C_REAL(h_GPU_HC_Track_Sols[19]);
  Transl21[2] = MAGMA_C_REAL(h_GPU_HC_Track_Sols[20]);
  MVG_Utility->Normalize_Translation_Vector( Transl21 );
  std::copy(Transl21, Transl21 + 3, begin(normalized_t21));
  
  //> \transl_{31}
  Transl31[0] = MAGMA_C_REAL(h_GPU_HC_Track_Sols[21]);
  Transl31[1] = MAGMA_C_REAL(h_GPU_HC_Track_Sols[22]);
  Transl31[2] = MAGMA_C_REAL(h_GPU_HC_Track_Sols[23]);
  MVG_Utility->Normalize_Translation_Vector( Transl31 );
  std::copy(Transl31, Transl31 + 3, begin(normalized_t31));
}

void Evaluations::Convert_Trifocal_Rotation( magmaFloatComplex *h_GPU_HC_Track_Sols ) {
  //> \rot_{21}
  Rot21[0] = MAGMA_C_REAL(h_GPU_HC_Track_Sols[24]);
  Rot21[1] = MAGMA_C_REAL(h_GPU_HC_Track_Sols[25]);
  Rot21[2] = MAGMA_C_REAL(h_GPU_HC_Track_Sols[26]);
  MVG_Utility->Cayley_To_Rotation_Matrix( Rot21, Sol_Rotm_21 );
  std::copy(Sol_Rotm_21, Sol_Rotm_21 + 9, begin(normalized_R21));
  
  //> \rot_{31}
  Rot31[0] = MAGMA_C_REAL(h_GPU_HC_Track_Sols[27]);
  Rot31[1] = MAGMA_C_REAL(h_GPU_HC_Track_Sols[28]);
  Rot31[2] = MAGMA_C_REAL(h_GPU_HC_Track_Sols[29]);
  MVG_Utility->Cayley_To_Rotation_Matrix( Rot31, Sol_Rotm_31 );
  std::copy(Sol_Rotm_31, Sol_Rotm_31 + 9, begin(normalized_R31));
}

void Evaluations::Check_Deviations_of_Veridical_Sol_from_GT( magmaFloatComplex *h_GPU_HC_Track_Sols, float GT_Pose21[12], float GT_Pose31[12] ) {
  //> Retrieve translation. The results are normalized_t21 and normalized_t31.
  Convert_Trifocal_Translation( h_GPU_HC_Track_Sols );

  //> Retrieve rotation. The results are normalized_R21 and normalized_R31.
  Convert_Trifocal_Rotation( h_GPU_HC_Track_Sols );

  //> Decompose the GT pose into rotation and translation
  get_GT_Rotation( GT_Pose21, GT_Rot21 );
  get_GT_Rotation( GT_Pose31, GT_Rot31 );
  get_GT_Translation( GT_Pose21, GT_Transl21 );
  get_GT_Translation( GT_Pose31, GT_Transl31 );

  //> Normalize the GT translations
  MVG_Utility->Normalize_Translation_Vector( GT_Transl21 );
  MVG_Utility->Normalize_Translation_Vector( GT_Transl31 );
  std::cout << "GT translation_21 = (" << GT_Transl21[0] << ", " << GT_Transl21[1] << ", " << GT_Transl21[2] << ")" << std::endl;
  std::cout << "GT translation_31 = (" << GT_Transl31[0] << ", " << GT_Transl31[1] << ", " << GT_Transl31[2] << ")" << std::endl;
  std::cout << "Sol translation_21 = (" << normalized_t21[0] << ", " << normalized_t21[1] << ", " << normalized_t21[2] << ")" << std::endl;
  std::cout << "Sol translation_31 = (" << normalized_t31[0] << ", " << normalized_t31[1] << ", " << normalized_t31[2] << ")" << std::endl;

  //> Calculate the residuals 
  float residual_R21 = get_Rotation_Residual( GT_Rot21, normalized_R21 );
  float residual_R31 = get_Rotation_Residual( GT_Rot31, normalized_R31 );
  float residual_t21 = get_Translation_Residual( GT_Transl21, normalized_t21 );
  float residual_t31 = get_Translation_Residual( GT_Transl31, normalized_t31 );

  std::cout << "Residuals in Rotations:    (R21) " << residual_R21 << " (R31) " << residual_R31 << std::endl;
  std::cout << "Residuals in Translations: (t21) " << residual_t21 << " (t31) " << residual_t31 << std::endl;
}

void Evaluations::Transform_GPUHC_Sols_to_Trifocal_Relative_Pose( \
  magmaFloatComplex *h_GPU_HC_Track_Sols, bool *h_is_GPU_HC_Sol_Converge, float *IntrinsicMatrix ) {
//> ----------------------------------------------------------------------------------
//> !!Note!! This function is specifically desgined for trifocal 2op1p 30x30 problem
//> ----------------------------------------------------------------------------------

  for (int i = 0; i < 9; i++) K[i] = IntrinsicMatrix[i];

  unsigned RANSAC_loop_index;
  unsigned track_index_in_one_RANSAC;

  //> Loop over all paths
  for (int bs = 0; bs < num_of_tracks*NUM_OF_RANSAC_ITERATIONS; bs++) {

    //> Current RANSAC loop index
    RANSAC_loop_index = std::floor(bs / num_of_tracks);
    track_index_in_one_RANSAC = bs - RANSAC_loop_index*num_of_tracks;
    
    //> if the solution converges
    if ( (h_is_GPU_HC_Sol_Converge + num_of_tracks * RANSAC_loop_index)[ bs ] ) {

      //> solution index offset under a RANSAC loop
      unsigned RANSAC_Sol_Offset = num_of_tracks * (num_of_variables+1) * RANSAC_loop_index + track_index_in_one_RANSAC * (num_of_variables+1);
      int small_imag_part_counter = 0;
      int positive_depth_counter = 0;

      //> Check the imaginary part of the two relative rotations
      for (int vi = 0; vi < 6; vi++) {
          if ( fabs(MAGMA_C_IMAG((h_GPU_HC_Track_Sols + RANSAC_Sol_Offset)[24 + vi])) < IMAG_PART_TOL ) 
            small_imag_part_counter++;
      }
      if ( small_imag_part_counter < 6 ) continue;

      //> Check whether the depths are all positive
      for (int di = 0; di < 8; di++) {
          if ( MAGMA_C_REAL((h_GPU_HC_Track_Sols + RANSAC_Sol_Offset)[di]) >= 0 ) positive_depth_counter++;
      }
      if (positive_depth_counter < 8) continue;

      Convert_Trifocal_Translation( h_GPU_HC_Track_Sols );
      normalized_t21s.push_back( normalized_t21 );
      normalized_t31s.push_back( normalized_t31 );

      Convert_Trifocal_Rotation( h_GPU_HC_Track_Sols );
      normalized_R21s.push_back( normalized_R21 );
      normalized_R31s.push_back( normalized_R31 );

      //> Compute the fundamental matrix used to find the maximal inliers support
      //> F21 and F31
      MVG_Utility->get_Fundamental_Matrix( K, Rot21, Transl21 );
      std::copy( MVG_Utility->F, MVG_Utility->F + 9, begin(FundMatrix21));
      F21s.push_back( FundMatrix21 );
      MVG_Utility->get_Fundamental_Matrix( K, Rot31, Transl31 );
      std::copy( MVG_Utility->F, MVG_Utility->F + 9, begin(FundMatrix31));
      F31s.push_back( FundMatrix31 );

      //> Push back track index with real solutions. Use for GPU DEBUG
      real_track_indices.push_back(bs);
    }
  }
}

float Evaluations::get_Rotation_Residual(float* GT_R, std::array<float, 9> Sol_R) {
  for (int i = 0; i < 9; i++) {
    R_gt_R[i] = GT_R[i];
    Sols_R_[i] = Sol_R[i];
  }

  MVG_Utility->get_Matrix_Transpose< 3 >(GT_R);                         //> R_{gt}'
  MVG_Utility->get_Matrix_Matrix_Product< 3 >(GT_R, Sols_R_, R_gt_R);    //> R_{gt}' * R
  float trace_RgR = MVG_Utility->get_Matrix_Trace< 3 >(R_gt_R);         //> trace( R_{gt}' * R )

  //> Remember to transpose the GT rotation back
  MVG_Utility->get_Matrix_Transpose< 3 >(GT_R);                         //> R_{gt}

  return acos(0.5 * (trace_RgR - 1.0));
}

float Evaluations::get_Translation_Residual(float* GT_Transl, std::array<float, 3> Sol_Transl) {
  for (int i = 0; i < 3; i++) Sol_Transl_[i] = Sol_Transl[i];
  float Transl_Dot_Prod = MVG_Utility->get_Vector_Dot_Product< 3 >(GT_Transl, Sol_Transl_);
  return std::fabs(Transl_Dot_Prod - 1.0);
}

bool Evaluations::get_Solution_with_Maximal_Support( unsigned Num_Of_Triplet_Edgels, float* h_Triplet_Edge_Locations, float* h_Triplet_Edge_Tangents, float *K ) {

  std::array<float, 3> Rel_t21;
  std::array<float, 3> Rel_t31;
  std::array<float, 9> Rel_R21;
  std::array<float, 9> Rel_R31;

  float *gamma1   = new float[3];
  float *gamma2   = new float[3];
  float *gamma3   = new float[3];
  float *tangent1 = new float[3];
  float *tangent2 = new float[3];
  float *tangent3 = new float[3];
  float *hypo_R21 = new float[9];
  float *hypo_R31 = new float[9];
  float *hypo_t21 = new float[3];
  float *hypo_t31 = new float[3];
  float rho1_21, rho1_31;
  float reproj_error_21, reproj_error_31;

  gamma1[2]   = 1.0; gamma2[2]   = 1.0; gamma3[2]   = 1.0;
  tangent1[2] = 0.0; tangent2[2] = 0.0; tangent3[2] = 0.0;

  //> Initialize
  Max_Num_Of_Reproj_Inliers_Views21 = 0;
  Max_Num_Of_Reproj_Inliers_Views31 = 0;

  //> Loop over all hypothesis pose
  for ( int cp_i = 0; cp_i < normalized_t21s.size(); cp_i++ ) {

    //> Reset inlier counter
    Num_Of_Reproj_Err_Inliers_Views21 = 0;
    Num_Of_Reproj_Err_Inliers_Views31 = 0;

    //> Fetch hypothesis trifocal relative pose can convert data type
    Rel_R21 = normalized_R21s[cp_i];
    Rel_R31 = normalized_R31s[cp_i];
    Rel_t21 = normalized_t21s[cp_i];
    Rel_t31 = normalized_t31s[cp_i];
    std::copy(Rel_R21.begin(), Rel_R21.end(), hypo_R21);
    std::copy(Rel_R31.begin(), Rel_R31.end(), hypo_R31);
    std::copy(Rel_t21.begin(), Rel_t21.end(), hypo_t21);
    std::copy(Rel_t31.begin(), Rel_t31.end(), hypo_t31);

    //> Loop over all triplet edgels
    for ( int ei = 0; ei < Num_Of_Triplet_Edgels; ei++ ) {
      //> Assigne triplet edgels
      gamma1[0] = h_Triplet_Edge_Locations(ei, 0);
      gamma1[1] = h_Triplet_Edge_Locations(ei, 1);
      gamma2[0] = h_Triplet_Edge_Locations(ei, 2);
      gamma2[1] = h_Triplet_Edge_Locations(ei, 3);
      gamma3[0] = h_Triplet_Edge_Locations(ei, 4);
      gamma3[1] = h_Triplet_Edge_Locations(ei, 5);

      tangent1[0] = h_Triplet_Edge_Tangents(ei, 0);
      tangent1[1] = h_Triplet_Edge_Tangents(ei, 1);
      tangent2[0] = h_Triplet_Edge_Tangents(ei, 2);
      tangent2[1] = h_Triplet_Edge_Tangents(ei, 3);
      tangent3[0] = h_Triplet_Edge_Tangents(ei, 4);
      tangent3[1] = h_Triplet_Edge_Tangents(ei, 5);

      //> Use Sampson error / reprojection error / tangent transfer error ?
      //> Reprojection error of view 1&2
      rho1_21 = MVG_Utility->get_depth_rho( gamma1, gamma2, hypo_R21, hypo_t21 );
      reproj_error_21 = MVG_Utility->get_Reprojection_Pixels_Error( gamma1, gamma2, hypo_R21, hypo_t21, K, rho1_21 );

      //> Reprojection error of view 1&3
      rho1_31 = MVG_Utility->get_depth_rho( gamma1, gamma3, hypo_R31, hypo_t31 );
      reproj_error_31 = MVG_Utility->get_Reprojection_Pixels_Error( gamma1, gamma3, hypo_R31, hypo_t31, K, rho1_31 );

      if ( reproj_error_21 < REPROJ_ERROR_INLIER_THRESH ) Num_Of_Reproj_Err_Inliers_Views21++;
      if ( reproj_error_31 < REPROJ_ERROR_INLIER_THRESH ) Num_Of_Reproj_Err_Inliers_Views31++;
    }

    //> Record maximal inlier support for view 1&2 and view 1&3 individually
    if (Num_Of_Reproj_Err_Inliers_Views21 >= Max_Num_Of_Reproj_Inliers_Views21) {
    // if (Num_Of_Reproj_Err_Inliers_Views21 == Num_Of_Triplet_Edgels) {
      Max_Num_Of_Reproj_Inliers_Views21 = Num_Of_Reproj_Err_Inliers_Views21;
      // Max_Reproj_Inliers_Support_Views21_Index = cp_i;
      Max_Reproj_Inliers_Support_Views21_Index.push_back( cp_i );
    }
    if (Num_Of_Reproj_Err_Inliers_Views31 >= Max_Num_Of_Reproj_Inliers_Views31) {
    // if (Num_Of_Reproj_Err_Inliers_Views31 == Num_Of_Triplet_Edgels) {
      Max_Num_Of_Reproj_Inliers_Views31 = Num_Of_Reproj_Err_Inliers_Views31;
      // Max_Reproj_Inliers_Support_Views31_Index = cp_i;
      Max_Reproj_Inliers_Support_Views31_Index.push_back( cp_i );
    }
  }

  //> Delete intermediate arrays
  delete [] gamma1;
  delete [] gamma2;
  delete [] gamma3;
  delete [] tangent1;
  delete [] tangent2;
  delete [] tangent3;

  delete [] hypo_R21;
  delete [] hypo_R31;
  delete [] hypo_t21;
  delete [] hypo_t31;

  if ( Max_Reproj_Inliers_Support_Views21_Index.size() == 0 || Max_Reproj_Inliers_Support_Views31_Index.size() == 0 ) {
    return false;
  }
  else {
    // std::cout << "Index of poses with max number of inliers views 1&2: ";
    // for (int i = 0; i < Max_Reproj_Inliers_Support_Views21_Index.size(); i++) std::cout << Max_Reproj_Inliers_Support_Views21_Index[i] << ", ";
    // std::cout << std::endl;
    // std::cout << "Index of poses with max number of inliers views 1&3: ";
    // for (int i = 0; i < Max_Reproj_Inliers_Support_Views31_Index.size(); i++) std::cout << Max_Reproj_Inliers_Support_Views31_Index[i] << ", ";
    // std::cout << std::endl;
    if (Max_Reproj_Inliers_Support_Views21_Index.size() > 0 && Max_Reproj_Inliers_Support_Views31_Index.size() > 0)
      std::cout << "### Found GT pose!" << std::endl;

    //> Now with solution index supported by maximal inliers, fetch the corresponding solution(s)
    R21_w_Max_Supports = normalized_R21s[ Max_Reproj_Inliers_Support_Views21_Index[0] ];
    t21_w_Max_Supports = normalized_t21s[ Max_Reproj_Inliers_Support_Views21_Index[0] ];
    R31_w_Max_Supports = normalized_R31s[ Max_Reproj_Inliers_Support_Views31_Index[0] ];
    t31_w_Max_Supports = normalized_t31s[ Max_Reproj_Inliers_Support_Views31_Index[0] ];
    return true;
  }
}

void Evaluations::get_HC_Steps_of_Actual_Sols( magmaFloatComplex *h_Debug_Purpose ) {

  int fetch_HC_steps;

  //> Find the union of HC steps from views 1&2 and 1&3
  std::vector<int>::iterator it, st;
  std::vector<int> union_HC_steps( Max_Reproj_Inliers_Support_Views21_Index.size() + Max_Reproj_Inliers_Support_Views31_Index.size());
  it = std::set_union( Max_Reproj_Inliers_Support_Views21_Index.begin(), Max_Reproj_Inliers_Support_Views21_Index.end(), \
                       Max_Reproj_Inliers_Support_Views31_Index.begin(), Max_Reproj_Inliers_Support_Views31_Index.end(), \
                       union_HC_steps.begin() );
  
  for (st = union_HC_steps.begin(); st != it; ++st) {
    fetch_HC_steps = MAGMA_C_REAL( h_Debug_Purpose[ real_track_indices[ *st ] ] );
    HC_steps_of_actual_solutions.push_back( fetch_HC_steps );
  } 
}

void Evaluations::Measure_Relative_Pose_Error( float GT_Pose21[12], float GT_Pose31[12] ) {
  //> Decompose the GT pose into rotation and translation
  get_GT_Rotation( GT_Pose21, GT_Rot21 );
  get_GT_Rotation( GT_Pose31, GT_Rot31 );
  get_GT_Translation( GT_Pose21, GT_Transl21 );
  get_GT_Translation( GT_Pose31, GT_Transl31 );

  //> Normalize the GT translations
  MVG_Utility->Normalize_Translation_Vector( GT_Transl21 );
  MVG_Utility->Normalize_Translation_Vector( GT_Transl31 );

  Min_Residual_R21 = get_Rotation_Residual( GT_Rot21, R21_w_Max_Supports );
  Min_Residual_R31 = get_Rotation_Residual( GT_Rot31, R31_w_Max_Supports );
  Min_Residual_t21 = get_Translation_Residual( GT_Transl21, t21_w_Max_Supports );
  Min_Residual_t31 = get_Translation_Residual( GT_Transl31, t31_w_Max_Supports );

  if (Min_Residual_t21 < TRANSL_RESIDUAL_TOL && Min_Residual_t31 < TRANSL_RESIDUAL_TOL && \
      Min_Residual_R21 < ROT_RESIDUAL_TOL && Min_Residual_R31 < ROT_RESIDUAL_TOL) {
    success_flag = true;
  }
}

void Evaluations::Measure_Relative_Pose_Error_from_All_Real_Sols( float GT_Pose21[12], float GT_Pose31[12], magmaFloatComplex *h_Debug_Purpose ) {

  //> Decompose the GT pose into rotation and translation
  get_GT_Rotation( GT_Pose21, GT_Rot21 );
  get_GT_Rotation( GT_Pose31, GT_Rot31 );
  get_GT_Translation( GT_Pose21, GT_Transl21 );
  get_GT_Translation( GT_Pose31, GT_Transl31 );

  //> Normalize the GT translations
  MVG_Utility->Normalize_Translation_Vector( GT_Transl21 );
  MVG_Utility->Normalize_Translation_Vector( GT_Transl31 );

  std::cout << " ====================================================== " << std::endl;

  //> Measure the relative pose error only when there is a valid solution
  if (!normalized_R21s.empty()) {
    
    float Residual_R21, Residual_R31;
    float Residual_t21, Residual_t31;
    for (int si = 0; si < normalized_R21s.size(); si++) {

      Residual_R21 = get_Rotation_Residual(GT_Rot21, normalized_R21s[si]);
      Residual_R31 = get_Rotation_Residual(GT_Rot31, normalized_R31s[si]);
      Residual_t21 = get_Translation_Residual(GT_Transl21, normalized_t21s[si]);
      Residual_t31 = get_Translation_Residual(GT_Transl31, normalized_t31s[si]);

      if (Residual_R21 < Min_Residual_R21) Min_Residual_R21 = Residual_R21;
      if (Residual_R31 < Min_Residual_R31) Min_Residual_R31 = Residual_R31;
      if (Residual_t21 < Min_Residual_t21) Min_Residual_t21 = Residual_t21;
      if (Residual_t31 < Min_Residual_t31) Min_Residual_t31 = Residual_t31;

      if (Residual_t21 < TRANSL_RESIDUAL_TOL && Residual_t31 < TRANSL_RESIDUAL_TOL && \
          Residual_R21 < ROT_RESIDUAL_TOL && Residual_R31 < ROT_RESIDUAL_TOL) {
            
          success_flag = true;
            
          //> HC steps, if GPU_DEBUG is activated
// #if GPU_DEBUG
          // std::cout << si << ", ";
          // int fetch_HC_steps = MAGMA_C_REAL( h_Debug_Purpose[ real_track_indices[si] ] );
          // HC_steps_of_actual_solutions.push_back( fetch_HC_steps );          
// #endif
      }
    }
  }
}

Evaluations::~Evaluations() {
  //> Close all files
  HC_Track_Sols_File.close();
  HC_Actual_Sols_Steps_File.close();

  //> Free memory
  delete [] Rot21;
  delete [] Rot31;
  delete [] Transl21;
  delete [] Transl31;

  delete [] GT_Rot21;
  delete [] GT_Rot31;
  delete [] GT_Transl21;
  delete [] GT_Transl31;

  delete [] Sol_Rotm_21;
  delete [] Sol_Rotm_31;

  delete [] R_gt_R;
  delete [] Sols_R_;
  delete [] Sol_Transl_;
}


#endif
