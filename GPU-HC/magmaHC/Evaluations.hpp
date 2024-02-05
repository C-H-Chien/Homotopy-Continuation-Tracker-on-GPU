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
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <chrono>

//> MAGMA
#include "magma_v2.h"

class Evaluations {
    
public:
    //> Constructor
    Evaluations();

    //> Destructor
    ~Evaluations();

    //> Write data to files
    void Write_Converged_Sols( magmaFloatComplex *h_GPU_HC_Track_Sols, bool *h_is_GPU_HC_Sol_Converge );

    //> Evaluate GPU-HC Solutions
    void Evaluate_Sols( magmaFloatComplex *h_GPU_HC_Track_Sols, bool *h_is_GPU_HC_Sol_Converge, bool *h_is_GPU_HC_Sol_Infinity );
    

    //> Some evaluation data
    unsigned int Num_Of_Coverged_Sols;
    unsigned int Num_Of_Inf_Sols;
    unsigned int Num_Of_Real_Sols;

private:
    //> output streams for files to be written
    std::ofstream GPUHC_Track_Sols_File;
};

#endif
