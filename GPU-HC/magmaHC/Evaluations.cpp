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
Evaluations::Evaluations( std::string Output_Files_Path, int num_of_tracks, int num_of_vars )
  : WRITE_FILES_PATH(Output_Files_Path), num_of_tracks(num_of_tracks), num_of_variables(num_of_vars)
{
  //> Initialize to zeros
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

  //> Write successful HC track solutions to files
  std::string write_sols_file_dir = WRITE_FILES_PATH.append("GPU_Converged_HC_tracks.txt");
  GPUHC_Track_Sols_File.open(write_sols_file_dir);
  if ( !GPUHC_Track_Sols_File.is_open() ) LOG_FILE_ERROR("write_sols_file_dir");

  //> util class
  MVG_Utility = std::shared_ptr<util>(new util());

  //> Allocate relative pose arrays
  Rot21     = new float[9];
  Rot31     = new float[9];
  Transl21  = new float[3];
  Transl31  = new float[3];

  Sol_Rotm_21 = new float[9];
  Sol_Rotm_31 = new float[9];

  Sol_Transl_ = new float[3];
  R_gt_R      = new float[9];
  Sols_R_     = new float[9];
}

void Evaluations::Write_Converged_Sols( \
    magmaFloatComplex *h_GPU_HC_Track_Sols, \
    bool *h_is_GPU_HC_Sol_Converge ) 
{
  for (int ri = 0; ri < NUM_OF_RANSAC_ITERATIONS; ri++) {
    GPUHC_Track_Sols_File << "-------------------- RANSAC Iteration " << ri+1 << " --------------------\n\n";
    for (int bs = 0; bs < num_of_tracks; bs++) {
      GPUHC_Track_Sols_File << std::setprecision(10);

      if (h_is_GPU_HC_Sol_Converge[ bs ] == 1) {
        //GPUHC_Track_Sols_File << h_is_GPU_HC_Sol_Converge[ bs ] << "\n";
        for (int vs = 0; vs < num_of_variables; vs++) {
          GPUHC_Track_Sols_File << std::setprecision(20) << MAGMA_C_REAL((h_GPU_HC_Track_Sols + ri * num_of_tracks * (num_of_variables+1) + bs * (num_of_variables+1))[vs]) << "\t" \
                                << std::setprecision(20) << MAGMA_C_IMAG((h_GPU_HC_Track_Sols + ri * num_of_tracks * (num_of_variables+1) + bs * (num_of_variables+1))[vs]) << "\n";
        }
        GPUHC_Track_Sols_File << "\n";
      }
    }
    GPUHC_Track_Sols_File << "\n";
  }
}

void Evaluations::Evaluate_GPUHC_Sols( \
    magmaFloatComplex *h_GPU_HC_Track_Sols, \
    bool *h_is_GPU_HC_Sol_Converge, \
    bool *h_is_GPU_HC_Sol_Infinity, \
    int ransac_sample_offset ) 
{
  //> Count the number of converged solutions, the number of infinity failed solutions, and the number of real solutions
  for (int bs = 0; bs < num_of_tracks; bs++) {
    if ( (h_is_GPU_HC_Sol_Converge + num_of_tracks * ransac_sample_offset)[ bs ] ) Num_Of_Coverged_Sols++;
    if ( (h_is_GPU_HC_Sol_Infinity + num_of_tracks * ransac_sample_offset)[ bs ] ) Num_Of_Inf_Sols++;

    int Num_Of_Real_Vars = 0;
    if ((h_is_GPU_HC_Sol_Converge + num_of_tracks * ransac_sample_offset)[ bs ] == 1) {
      for (int vs = 0; vs < num_of_variables; vs++) {
        if (fabs(MAGMA_C_IMAG((h_GPU_HC_Track_Sols + num_of_tracks * (num_of_variables+1) * ransac_sample_offset + bs * (num_of_variables+1))[vs])) <= ZERO_IMAG_PART_TOL_FOR_SP) {
            Num_Of_Real_Vars++;
        }
      }
    }

    if (Num_Of_Real_Vars == num_of_variables) Num_Of_Real_Sols++;
  }
}

void Evaluations::Evaluate_RANSAC_GPUHC_Sols( \
    magmaFloatComplex *h_GPU_HC_Track_Sols, \
    bool *h_is_GPU_HC_Sol_Converge, \
    bool *h_is_GPU_HC_Sol_Infinity )
{
  //> Loop over all RANSAC iterations
  for (int ri = 0; ri < NUM_OF_RANSAC_ITERATIONS; ri++) {
    Evaluate_GPUHC_Sols( h_GPU_HC_Track_Sols, h_is_GPU_HC_Sol_Converge, h_is_GPU_HC_Sol_Infinity, ri );
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

void Evaluations::Transform_GPUHC_Sols_to_Trifocal_Relative_Pose( \
  magmaFloatComplex *h_GPU_HC_Track_Sols, bool *h_is_GPU_HC_Sol_Converge, float IntrinsicMatrix[9] ) {
//> ----------------------------------------------------------------------------------
//> !!Note!! This function is specifically desgined for trifocal 2op1p 30x30 problem
//> ----------------------------------------------------------------------------------

  for (int i = 0; i < 9; i++) K[i] = IntrinsicMatrix[i];

  float norm_t21, norm_t31;
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

      //> \transl_{21}
      Transl21[0] = MAGMA_C_REAL((h_GPU_HC_Track_Sols + RANSAC_Sol_Offset)[18]);
      Transl21[1] = MAGMA_C_REAL((h_GPU_HC_Track_Sols + RANSAC_Sol_Offset)[19]);
      Transl21[2] = MAGMA_C_REAL((h_GPU_HC_Track_Sols + RANSAC_Sol_Offset)[20]);
      MVG_Utility->Normalize_Translation_Vector( Transl21 );
      std::copy(Transl21, Transl21 + 3, begin(normalized_t21));
      normalized_t21s.push_back( normalized_t21 );

      //> \transl_{31}
      Transl31[0] = MAGMA_C_REAL((h_GPU_HC_Track_Sols + RANSAC_Sol_Offset)[21]);
      Transl31[1] = MAGMA_C_REAL((h_GPU_HC_Track_Sols + RANSAC_Sol_Offset)[22]);
      Transl31[2] = MAGMA_C_REAL((h_GPU_HC_Track_Sols + RANSAC_Sol_Offset)[23]);
      MVG_Utility->Normalize_Translation_Vector( Transl31 );
      std::copy(Transl31, Transl31 + 3, begin(normalized_t31));
      normalized_t31s.push_back( normalized_t31 );

      //> \rot_{21}
      Rot21[0] = MAGMA_C_REAL((h_GPU_HC_Track_Sols + RANSAC_Sol_Offset)[24]);
      Rot21[1] = MAGMA_C_REAL((h_GPU_HC_Track_Sols + RANSAC_Sol_Offset)[25]);
      Rot21[2] = MAGMA_C_REAL((h_GPU_HC_Track_Sols + RANSAC_Sol_Offset)[26]);
      MVG_Utility->Cayley_To_Rotation_Matrix( Rot21, Sol_Rotm_21 );
      std::copy(Sol_Rotm_21, Sol_Rotm_21 + 9, begin(normalized_R21));
      normalized_R21s.push_back( normalized_R21 );

      //> \rot_{31}
      Rot31[0] = MAGMA_C_REAL((h_GPU_HC_Track_Sols + RANSAC_Sol_Offset)[27]);
      Rot31[1] = MAGMA_C_REAL((h_GPU_HC_Track_Sols + RANSAC_Sol_Offset)[28]);
      Rot31[2] = MAGMA_C_REAL((h_GPU_HC_Track_Sols + RANSAC_Sol_Offset)[29]);
      MVG_Utility->Cayley_To_Rotation_Matrix( Rot31, Sol_Rotm_31 );
      std::copy(Sol_Rotm_31, Sol_Rotm_31 + 9, begin(normalized_R31));
      normalized_R31s.push_back( normalized_R31 );

      //> Compute the fundamental matrix used to find the maximal inliers support
      //> F21 and F31
      MVG_Utility->get_Fundamental_Matrix( K, Rot21, Transl21 );
      std::copy( MVG_Utility->F, MVG_Utility->F + 9, begin(FundMatrix21));
      F21s.push_back( FundMatrix21 );
      MVG_Utility->get_Fundamental_Matrix( K, Rot31, Transl31 );
      std::copy( MVG_Utility->F, MVG_Utility->F + 9, begin(FundMatrix31));
      F31s.push_back( FundMatrix31 );

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
  return acos(0.5 * (trace_RgR - 1.0));
}

float Evaluations::get_Translation_Residual(float* GT_Transl, std::array<float, 3> Sol_Transl) {
  for (int i = 0; i < 3; i++) Sol_Transl_[i] = Sol_Transl[i];
  float Transl_Dot_Prod = MVG_Utility->get_Vector_Dot_Product< 3 >(GT_Transl, Sol_Transl_);
  return std::fabs(Transl_Dot_Prod - 1.0);
}

void Evaluations::Measure_Relative_Pose_Error( float GT_Pose21[12], float GT_Pose31[12] ) {
  //> Decompose the GT pose into rotation and translation
  float* GT_Rot21     = new float[9];
  float* GT_Rot31     = new float[9];
  float* GT_Transl21  = new float[3];
  float* GT_Transl31  = new float[3];

  for (int i = 0; i < 9; i++) {
    GT_Rot21[i] = GT_Pose21[i];
    GT_Rot31[i] = GT_Pose31[i];
  }

  int t_index = 0;
  for (int i = 9; i < 12; i++) {
    GT_Transl21[t_index] = GT_Pose21[i];
    GT_Transl31[t_index] = GT_Pose31[i];
    t_index++;
  }

  //> Normalize the GT translations
  MVG_Utility->Normalize_Translation_Vector( GT_Transl21 );
  MVG_Utility->Normalize_Translation_Vector( GT_Transl31 );

  //> Measure the relative pose error only when there is a valid solution
  if (!normalized_R21s.empty()) {
    
    float Residual_R21, Residual_R31;
    float Residual_t21, Residual_t31;
    for (int si = 0; si < normalized_R21s.size(); si++) {

      Residual_R21 = get_Rotation_Residual(GT_Rot21, normalized_R21s[si]);
      Residual_R31 = get_Rotation_Residual(GT_Rot31, normalized_R31s[si]);
      Residual_t21 = get_Translation_Residual(GT_Transl21, normalized_t21s[si]);
      Residual_t31 = get_Translation_Residual(GT_Transl31, normalized_t31s[si]);

      if (Residual_t21 < TRANSL_RESIDUAL_TOL && Residual_t31 < TRANSL_RESIDUAL_TOL && \
          Residual_R21 < ROT_RESIDUAL_TOL && Residual_R31 < ROT_RESIDUAL_TOL) {
            if (Residual_R21 < Min_Residual_R21) Min_Residual_R21 = Residual_R21;
            if (Residual_R31 < Min_Residual_R31) Min_Residual_R31 = Residual_R31;
            if (Residual_t21 < Min_Residual_t21) Min_Residual_t21 = Residual_t21;
            if (Residual_t31 < Min_Residual_t31) Min_Residual_t31 = Residual_t31;
            success_flag = true;
      }
    }
  }
  
  delete [] GT_Rot21;
  delete [] GT_Rot31;
  delete [] GT_Transl21;
  delete [] GT_Transl31;
}

Evaluations::~Evaluations() {
  //> Close all files
  GPUHC_Track_Sols_File.close();

  //> Free memory
  delete [] Rot21;
  delete [] Rot31;
  delete [] Transl21;
  delete [] Transl31;

  delete [] Sol_Rotm_21;
  delete [] Sol_Rotm_31;

  delete [] R_gt_R;
  delete [] Sols_R_;
  delete [] Sol_Transl_;
}


#endif
