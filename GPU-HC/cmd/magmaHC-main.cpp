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

template<typename T>
bool run_GPU_HC_Solver( YAML::Node Problem_Settings_Map, int Data_Size_for_Indices) {

  double GPUHC_time_ms = 0.0;
  double avg_gpu_runtime = 0.0;
  double max_gpu_runtime = 0.0;
  double min_gpu_runtime = 10000.0;
  double all_gpu_runtime[TEST_RANSAC_TIMES];

  GPU_HC_Solver<T> GPU_HC_( Problem_Settings_Map, Data_Size_for_Indices );

  //> (1) Allocate CPU and GPU arrays
  GPU_HC_.Allocate_Arrays();

  //> Loop over many times of RANSACs
  for (int ti = 0; ti < TEST_RANSAC_TIMES; ti++) {
    //> (2) Read Problem-Specific Data
    bool pass_Data_Read_Test = GPU_HC_.Read_Problem_Data();
    if (!pass_Data_Read_Test) return false;

    //> (3) Read RANSAC Data
    pass_Data_Read_Test = GPU_HC_.Read_RANSAC_Data( ti );
    if (!pass_Data_Read_Test) return false;

    //> (4) Convert from triplet edgels to target parameters
    GPU_HC_.Prepare_Target_Params( ti );

    //> (5) Transfer data from CPU to GPU
    GPU_HC_.Data_Transfer_From_Host_To_Device();

    //> (6) Solve the problem by GPU-HC
    GPU_HC_.Solve_by_GPU_HC();

    //> (7) Free triplet edgels memory
    GPU_HC_.Free_Triplet_Edgels_Mem();

    GPUHC_time_ms = (GPU_HC_.gpu_max_time_from_multiple_GPUs)*1000;

    all_gpu_runtime[ti] = GPUHC_time_ms;
    avg_gpu_runtime += GPUHC_time_ms;
    max_gpu_runtime = (GPUHC_time_ms > max_gpu_runtime) ? (GPUHC_time_ms) : (max_gpu_runtime);
    min_gpu_runtime = (GPUHC_time_ms < min_gpu_runtime) ? (GPUHC_time_ms) : (min_gpu_runtime);

    //> Reset GPU max elapse time across all GPUs 
    GPU_HC_.gpu_max_time_from_multiple_GPUs = 0.0;
  }

  avg_gpu_runtime /= TEST_RANSAC_TIMES;
  std::cout << std::endl;
  printf("## Running %d rounds of %d RANSAC iterations:\n", TEST_RANSAC_TIMES, NUM_OF_RANSAC_ITERATIONS);
  printf(" - [Average GPU Computation Time] %7.2f (ms)\n", avg_gpu_runtime);
  printf(" - [Maximal GPU Computation Time] %7.2f (ms)\n", max_gpu_runtime);
  printf(" - [Minimal GPU Computation Time] %7.2f (ms)\n", min_gpu_runtime);

  //> Calculate the standard deviation
  double sigma = 0.0;
  for (int i = 0; i < TEST_RANSAC_TIMES; i++) {
    sigma += (all_gpu_runtime[i] - avg_gpu_runtime) * (all_gpu_runtime[i] - avg_gpu_runtime);
  }
  sigma /= TEST_RANSAC_TIMES;
  sigma = sqrt(sigma);
  printf(" - [Std dev GPU Computation Time] %7.2f (ms)\n", sigma);

  //> Export all data (HC steps for now) to files
  // GPU_HC_.Export_Data();

  return true;
}

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

  //> Get data dize (32-bit or 8-bit) for indices in evaluations
  int Data_Size_for_Indices = Problem_Settings_Map["Data_Size_for_Indices"].as<int>();

  //> Initialization from GPU-HC constructor
  bool should_continue = false;
  if (Data_Size_for_Indices == 8) {
    should_continue = run_GPU_HC_Solver<char>( Problem_Settings_Map, Data_Size_for_Indices );
  }
  else if (Data_Size_for_Indices == 32) {
    should_continue = run_GPU_HC_Solver<int>( Problem_Settings_Map, Data_Size_for_Indices );
  }
  else
    LOG_ERROR("Data size for indices defined in the YAML file is incorrect!");

  if (!should_continue)
    LOG_ERROR("Something's wrong with run_GPU_HC_Solver!");

  //> Maybe we can do something with should_continue

  return 0;
}
