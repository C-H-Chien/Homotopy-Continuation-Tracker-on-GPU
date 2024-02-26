#ifndef EVALUATIONS_CPP
#define EVALUATIONS_CPP
// =============================================================================================================================
//
// ChangLogs
//    24-02-05:   Initially created for definitions of functions in the ``Evaluations" class.
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

#include "definitions.hpp"
#include "Evaluations.hpp"

//> Constructor
Evaluations::Evaluations( ) {
  
  //> Initialize to zeros
  Num_Of_Inf_Sols = 0;
  Num_Of_Coverged_Sols = 0;
  Num_Of_Real_Sols = 0;
  Num_Of_Unique_Sols = 0;

  //> Prepare files for write:

  //> 1) Write successful HC track solutions to files
  std::string write_sols_file_dir = WRITE_FILES_PATH.append("GPU_Converged_HC_tracks.txt");
  GPUHC_Track_Sols_File.open(write_sols_file_dir);
  if ( !GPUHC_Track_Sols_File.is_open() ) LOG_FILE_ERROR("write_sols_file_dir");

  //> 2) ...
}

void Evaluations::Write_Converged_Sols( \
    magmaFloatComplex *h_GPU_HC_Track_Sols, \
    bool *h_is_GPU_HC_Sol_Converge ) 
{
  for (int bs = 0; bs < NUM_OF_TRACKS; bs++) {
    GPUHC_Track_Sols_File << std::setprecision(10);

    GPUHC_Track_Sols_File << h_is_GPU_HC_Sol_Converge[ bs ] << "\n";
    for (int vs = 0; vs < NUM_OF_VARS; vs++) {
      GPUHC_Track_Sols_File << std::setprecision(20) << MAGMA_C_REAL((h_GPU_HC_Track_Sols + bs * (NUM_OF_VARS+1))[vs]) << "\t" \
                            << std::setprecision(20) << MAGMA_C_IMAG((h_GPU_HC_Track_Sols + bs * (NUM_OF_VARS+1))[vs]) << "\n";
    }
    GPUHC_Track_Sols_File << "\n";
  }
}

void Evaluations::Evaluate_Sols( \
    magmaFloatComplex *h_GPU_HC_Track_Sols, \
    bool *h_is_GPU_HC_Sol_Converge, \
    bool *h_is_GPU_HC_Sol_Infinity ) 
{
  //> Count the number of converged solutions, the number of infinity failed solutions, and the number of real solutions
  for (int bs = 0; bs < NUM_OF_TRACKS; bs++) {
    if ( h_is_GPU_HC_Sol_Converge[ bs ] ) Num_Of_Coverged_Sols++;
    if ( h_is_GPU_HC_Sol_Infinity[ bs ] ) Num_Of_Inf_Sols++;

    int Num_Of_Real_Vars = 0;
    for (int vs = 0; vs < NUM_OF_VARS; vs++) {
      if (MAGMA_C_IMAG((h_GPU_HC_Track_Sols + bs * (NUM_OF_VARS+1))[vs]) <= ZERO_IMAG_PART_TOL_FOR_SP) {
          Num_Of_Real_Vars++;
      }
    }

    if (Num_Of_Real_Vars == NUM_OF_VARS) Num_Of_Real_Sols++;
  }
}

void Evaluations::Find_Unique_Sols( magmaFloatComplex *h_GPU_HC_Track_Sols, bool *h_is_GPU_HC_Sol_Converge ) {

  std::set< int > Duplicate_Sol_Index;
  std::set< int > Skip_Sol_Index;
  for (int bs = 0; bs < NUM_OF_TRACKS; bs++) {

    //> Make sure that the solution is converged before evaluating its uniqueness
    if ( h_is_GPU_HC_Sol_Converge[ bs ] ) {
      
      if (!Skip_Sol_Index.empty()) {

        //if(std::find(Skip_Sol_Index.begin(), Skip_Sol_Index.end(), bs) != Skip_Sol_Index.end())
        if (Skip_Sol_Index.find(bs) != Skip_Sol_Index.end())
          continue;

        //> Flush out the set
        Duplicate_Sol_Index.clear();
      }

      for (int ds = bs+1; ds < NUM_OF_TRACKS; ds++) {

        int Num_Of_Duplicate_Vars = 0;
        for (int vs = 0; vs < NUM_OF_VARS; vs++) {

          //> If both the real and imaginery parts are very close
          if ( Var_Diff_In_Real_Part(h_GPU_HC_Track_Sols, h_GPU_HC_Track_Sols, bs, ds, vs) < DUPLICATE_SOL_DIFF_TOL && \
               Var_Diff_In_Imag_Part(h_GPU_HC_Track_Sols, h_GPU_HC_Track_Sols, bs, ds, vs) < DUPLICATE_SOL_DIFF_TOL ) {
              Num_Of_Duplicate_Vars++;
          }
        }

        if (Num_Of_Duplicate_Vars == NUM_OF_VARS) Duplicate_Sol_Index.insert(ds);
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
}


#endif
