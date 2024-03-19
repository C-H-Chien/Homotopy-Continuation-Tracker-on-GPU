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

//> Macros
#include "magmaHC/definitions.hpp"
#include "magmaHC/GPU_HC_Solver.hpp"

int main(int argc, char **argv) {

  --argc; ++argv;
  std::string arg;
  int argIndx = 0, argTotal = 4;
  std::string REPO_PATH;

  if (argc) {
    arg = std::string(*argv);
    if (arg == "-h" || arg == "--help") {
      LOG_PRINT_HELP_MESSAGE;
      return 0;
    }
    else if (argc <= argTotal) {
      while(argIndx <= argTotal-1) {
        if (arg == "-d" || arg == "--directory") {
          argv++;
          arg = std::string(*argv);
          REPO_PATH = arg;
          argIndx+=2;
          break;
        }
        else {
          LOG_ERROR("Invalid input arguments!");
          LOG_PRINT_HELP_MESSAGE;
          return 0;
        }
        argv++;
      }
    }
    else if (argc > argTotal) {
      LOG_ERROR("Too many input arguments!");
      LOG_PRINT_HELP_MESSAGE;
      return 0;
    }
  }
  else {
    LOG_PRINT_HELP_MESSAGE;
    return 0;
  }

  //> Assertion failed if HC problem is undefined
  assert(UNDEFINE_HC_PROBLEM == true);

  //> Initialization from GPU-HC constructor
  GPU_HC_Solver GPU_HC_( REPO_PATH );

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
