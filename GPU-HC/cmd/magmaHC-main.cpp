#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <string>
#include <assert.h>
// =======================================================================================================
// main function
//
// Modifications
//    Chien  21-12-29    Initially created a parametric HC for a General Computer Vision Problems
//    Chien  22-11-14    Add geometric form of 5pt relative pose problem
//    Chien  22-11-16    Add algebraic form of 5pt relative pose problem
//    Chien  23-10-20    Add gamma trick and trifocal relative pose from lines at points problem
//    Chien  23-12-27    Add definitions.hpp placing all Macros
//    Chien  24-02-26    Shift most of the code to GPU_HC_Solver class. Make main code clean.
//
//> (c) LEMS, Brown University
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// =======================================================================================================

//> magma
#include "magma_v2.h"

//> Macros
#include "magmaHC/definitions.hpp"
#include "magmaHC/GPU_HC_Solver.hpp"

int main(int argc, char **argv) {

  //> Assertion failed if HC problem is undefined
  assert(UNDEFINE_HC_PROBLEM == true);

  //> Initialization from GPU-HC constructor
  GPU_HC_Solver GPU_HC_;

  //> (1) Allocate CPU and GPU arrays
  GPU_HC_.Allocate_Arrays();

  //> (2) Read Problem-Specific Data
  bool pass_Data_Read_Test = GPU_HC_.Read_Problem_Data();

  if (pass_Data_Read_Test) {
    
    //> (4) Compute and assign parameter homotopy coefficients
    GPU_HC_.Construct_Coeffs_From_Params();

    //> (5) Transfer data from CPU to GPU
    GPU_HC_.Data_Transfer_From_Host_To_Device();

    //> (6) Solve the problem by GPU-HC
    GPU_HC_.Solve_by_GPU_HC();
  }

  

  return 0;
}
