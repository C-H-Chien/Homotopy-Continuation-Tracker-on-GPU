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
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <set>


//> MAGMA
#include "magma_v2.h"

#include "definitions.hpp"

class Evaluations {
    
public:
    //> Constructor
    Evaluations( std::string );

    //> Destructor
    ~Evaluations();

    //> Write data to files
    void Write_Converged_Sols( magmaFloatComplex *h_GPU_HC_Track_Sols, bool *h_is_GPU_HC_Sol_Converge );

    //> Evaluate GPU-HC Solutions
    void Evaluate_Sols( magmaFloatComplex *h_GPU_HC_Track_Sols, bool *h_is_GPU_HC_Sol_Converge, bool *h_is_GPU_HC_Sol_Infinity );
    void Find_Unique_Sols( magmaFloatComplex *h_GPU_HC_Track_Sols, bool *h_is_GPU_HC_Sol_Converge );
    
    //> Some evaluation data
    unsigned Num_Of_Coverged_Sols;
    unsigned Num_Of_Inf_Sols;
    unsigned Num_Of_Real_Sols;
    unsigned Num_Of_Unique_Sols;

private:
    //> output streams for files to be written
    std::ofstream GPUHC_Track_Sols_File;

    float Var_Diff_In_Real_Part( magmaFloatComplex *Sol1, magmaFloatComplex *Sol2, int sol1_offset, int sol2_offset, int var_index ) {
        return std::fabs(MAGMA_C_REAL((Sol1 + sol1_offset*(NUM_OF_VARS+1))[var_index]) - MAGMA_C_REAL((Sol2 + sol2_offset*(NUM_OF_VARS+1))[var_index]));
    };
    float Var_Diff_In_Imag_Part( magmaFloatComplex *Sol1, magmaFloatComplex *Sol2, int sol1_offset, int sol2_offset, int var_index ) {
        return std::fabs(MAGMA_C_IMAG((Sol1 + sol1_offset*(NUM_OF_VARS+1))[var_index]) - MAGMA_C_IMAG((Sol2 + sol2_offset*(NUM_OF_VARS+1))[var_index]));
    };

    std::vector< int > Unique_Sols_Index;

    std::string WRITE_FILES_PATH;
};

#endif
