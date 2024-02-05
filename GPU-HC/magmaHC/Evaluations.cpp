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

#include <cuda.h>
#include <cuda_runtime.h>

//> MAGMA
#include "magma_v2.h"
#include "magma_lapack.h"
#include "magma_internal.h"

#include "definitions.hpp"
#include "Evaluations.hpp"


//> Constructor
Evaluations::Evaluations( ) {
    
    Num_Of_Inf_Sols = 0;
    Num_Of_Coverged_Sols = 0;
    Num_Of_Real_Sols = 0;

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


Evaluations::~Evaluations() {

    //> Close all files
    GPUHC_Track_Sols_File.close();
}


#endif
