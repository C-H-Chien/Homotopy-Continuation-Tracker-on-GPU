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
//    Chien  24-03-26    Use yaml-cpp to parse data from problem yaml files so that no recompilation is needed when switching problems
//
//> (c) LEMS, Brown University
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// =======================================================================================================

#include "magmaHC/definitions.hpp"
#include "magmaHC/GPU_HC_Solver.hpp"

#include <yaml-cpp/yaml.h>

int main(int argc, char **argv) {
  //> Get input argument
  --argc; ++argv;
  std::string arg;
  int argIndx = 0, argTotal = 4;
  std::string PROBLEM_NAME;

  if (argc) {
    arg = std::string(*argv);
    if (arg == "-h" || arg == "--help") {
      LOG_PRINT_HELP_MESSAGE;
      return 0;
    }
    else if (argc <= argTotal) {
      while(argIndx <= argTotal-1) {
        if (arg == "-p" || arg == "--problem") {
          argv++;
          arg = std::string(*argv);
          PROBLEM_NAME = arg;
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

  //> Path to GPU-HC problem settings yaml file
  std::string Problem_Path = "../../problems/" + PROBLEM_NAME + "/gpuhc_settings.yaml";

  YAML::Node Problem_Settings_Map;
  try {
		Problem_Settings_Map = YAML::LoadFile(Problem_Path);
#if SHOW_PROBLEM_SETTINGS
		std::cout << std::endl << Problem_Settings_Map << std::endl << std::endl;
#endif
	}
	catch (const std::exception& e) {
		std::cerr << "Exception: " << e.what() << std::endl;
    return 0;
	}

  //> Initialization from GPU-HC constructor
  GPU_HC_Solver GPU_HC_( Problem_Settings_Map );

  //> (1) Allocate CPU and GPU arrays
  GPU_HC_.Allocate_Arrays();

  //> Loop over many times of RANSACs
  for (int ti = 0; ti < TEST_RANSAC_TIMES; ti++) {
    //> (2) Read Problem-Specific Data
    bool pass_Data_Read_Test = GPU_HC_.Read_Problem_Data();
    if (!pass_Data_Read_Test) return 0;

    // //> (3) Read RANSAC Data, if required
    // if (GPU_HC_.RANSAC_Dataset_Name != "None") {
    //   pass_Data_Read_Test = GPU_HC_.Read_RANSAC_Data( ti );
    //   if (!pass_Data_Read_Test) return 0;
    // }

    // //> (4) Convert from triplet edgels to target parameters
    // GPU_HC_.Prepare_Target_Params();

    //> (5) Transfer data from CPU to GPU
    GPU_HC_.Data_Transfer_From_Host_To_Device();

    //> (6) Solve the problem by GPU-HC
    GPU_HC_.Solve_by_GPU_HC();
  }

  return 0;
}
