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

  //> Write successful HC track solutions to files
  std::string write_sols_file_dir = WRITE_FILES_PATH + "GPU_Converged_HC_tracks.txt";
  GPUHC_Track_Sols_File.open(write_sols_file_dir);
  if ( !GPUHC_Track_Sols_File.is_open() ) LOG_FILE_ERROR(write_sols_file_dir);

  //> Write successful HC track solutions to files
  std::string write_actual_sols_HC_steps_file_dir = WRITE_FILES_PATH + "HC_Steps_of_Actual_Solutions.txt";
  GPUHC_Actual_Sols_Steps_File.open(write_actual_sols_HC_steps_file_dir);
  if ( !GPUHC_Actual_Sols_Steps_File.is_open() ) LOG_FILE_ERROR(write_actual_sols_HC_steps_file_dir);

  //> util class
  MVG_Utility = std::shared_ptr<util>(new util());
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

Evaluations::~Evaluations() {
  //> Close all files
  GPUHC_Track_Sols_File.close();
  GPUHC_Actual_Sols_Steps_File.close();
}


#endif
