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
#include "magmaHC/CPU_HC_Solver.hpp"
#include <yaml-cpp/yaml.h>

//> ======================================
//> GPU-HC
//> ======================================
bool run_GPU_HC_Solver( YAML::Node Problem_Settings_Map) {

  double GPUHC_time_ms = 0.0;
  double avg_gpu_runtime = 0.0;
  double max_gpu_runtime = 0.0;
  double min_gpu_runtime = 10000.0;
  double all_gpu_runtime[TEST_RANSAC_TIMES];

  GPU_HC_Solver GPU_HC_( Problem_Settings_Map );

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

    //> (5) Set arrays used for early aborting RANSAC process
    GPU_HC_.Set_RANSAC_Abort_Arrays();

    //> (6) Transfer data from CPU to GPU
    GPU_HC_.Data_Transfer_From_Host_To_Device();

    //> (7) Set CUDA stream attribute values, if necessary
    GPU_HC_.Set_CUDA_Stream_Attributes();

    //> (8) Solve the problem by GPU-HC
    GPU_HC_.Solve_by_GPU_HC();

    //> (9) Free triplet edgels memory
    GPU_HC_.Free_Triplet_Edgels_Mem();

    //> (10) Free allocated arrays for early aborting RANSAC
    GPU_HC_.Free_Arrays_for_Aborting_RANSAC();

    // GPUHC_time_ms = (GPU_HC_.gpu_max_time_from_multiple_GPUs)*1000;
    GPUHC_time_ms = (GPU_HC_.multi_GPUs_time)*1000;

    all_gpu_runtime[ti] = GPUHC_time_ms;
    avg_gpu_runtime += GPUHC_time_ms;
    max_gpu_runtime = (GPUHC_time_ms > max_gpu_runtime) ? (GPUHC_time_ms) : (max_gpu_runtime);
    min_gpu_runtime = (GPUHC_time_ms < min_gpu_runtime) ? (GPUHC_time_ms) : (min_gpu_runtime);
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
  
  //> Write GPU-HC Timings
  std::string write_timings_file_dir = std::string("../../") + WRITE_FILES_FOLDER + std::string("GPU_Timings.txt");
  std::ofstream GPUHC_Timings_File;
  GPUHC_Timings_File.open(write_timings_file_dir);
  if ( !GPUHC_Timings_File.is_open() ) LOG_FILE_ERROR(write_timings_file_dir);
  for (int i = 0; i < TEST_RANSAC_TIMES; i++) {
    GPUHC_Timings_File << all_gpu_runtime[i] << "\n";
  }
  GPUHC_Timings_File.close();

  //> Write GPU-HC number of converged solutions, real solutions, and infinity-failed solutions
  std::string write_solution_statistics_dir = std::string("../../") + WRITE_FILES_FOLDER + std::string("GPU_Sols_Statistics.txt");
  std::ofstream GPUHC_Sols_Statistics_File;
  GPUHC_Sols_Statistics_File.open(write_solution_statistics_dir);
  if ( !GPUHC_Sols_Statistics_File.is_open() ) LOG_FILE_ERROR(write_solution_statistics_dir);
  for (int i = 0; i < TEST_RANSAC_TIMES; i++) {
    GPUHC_Sols_Statistics_File << GPU_HC_.Collect_Num_Of_Coverged_Sols[i] << "\t" \
                               << GPU_HC_.Collect_Num_Of_Inf_Sols[i]      << "\t" \
                               << GPU_HC_.Collect_Num_Of_Real_Sols[i]     << "\n";
  }
  GPUHC_Sols_Statistics_File.close();

  return true;
}

//> ======================================
//> CPU-HC
//> ======================================
bool run_CPU_HC_Solver( YAML::Node Problem_Settings_Map) {

  double CPUHC_time_ms = 0.0;
  double avg_cpu_runtime = 0.0;
  double max_cpu_runtime = 0.0;
  double min_cpu_runtime = 100000000000000.0;
  double all_cpu_runtime[TEST_RANSAC_TIMES];

  CPU_HC_Solver CPU_HC_( Problem_Settings_Map );

  //> (1) Allocate CPU and GPU arrays
  CPU_HC_.Allocate_Arrays();

  //> Loop over many times of RANSACs
  for (int ti = 0; ti < TEST_RANSAC_TIMES; ti++) {
    //> (2) Read Problem-Specific Data
    bool pass_Data_Read_Test = CPU_HC_.Read_Problem_Data();
    if (!pass_Data_Read_Test) return false;

    //> (3) Read RANSAC Data
    pass_Data_Read_Test = CPU_HC_.Read_RANSAC_Data( ti );
    if (!pass_Data_Read_Test) return false;

    //> (4) Convert from triplet edgels to target parameters
    CPU_HC_.Prepare_Target_Params( ti );

    CPU_HC_.Set_Initial_Array_Vals();

    //> (8) Solve the problem by GPU-HC
    CPU_HC_.Solve_by_CPU_HC();

    //> (9) Free triplet edgels memory
    CPU_HC_.Free_Triplet_Edgels_Mem();

    CPUHC_time_ms = (CPU_HC_.CPU_HC_time)*1000;

    all_cpu_runtime[ti] = CPUHC_time_ms;
    avg_cpu_runtime += CPUHC_time_ms;
    max_cpu_runtime = (CPUHC_time_ms > max_cpu_runtime) ? (CPUHC_time_ms) : (max_cpu_runtime);
    min_cpu_runtime = (CPUHC_time_ms < min_cpu_runtime) ? (CPUHC_time_ms) : (min_cpu_runtime);
  }

  avg_cpu_runtime /= TEST_RANSAC_TIMES;
  std::cout << std::endl;
  printf("## Running %d rounds of %d RANSAC iterations:\n", TEST_RANSAC_TIMES, NUM_OF_RANSAC_ITERATIONS);
  printf(" - [Average CPU Computation Time] %7.2f (ms)\n", avg_cpu_runtime);
  printf(" - [Maximal CPU Computation Time] %7.2f (ms)\n", max_cpu_runtime);
  printf(" - [Minimal CPU Computation Time] %7.2f (ms)\n", min_cpu_runtime);

  //> Calculate the standard deviation
  double sigma = 0.0;
  for (int i = 0; i < TEST_RANSAC_TIMES; i++) {
    sigma += (all_cpu_runtime[i] - avg_cpu_runtime) * (all_cpu_runtime[i] - avg_cpu_runtime);
  }
  sigma /= TEST_RANSAC_TIMES;
  sigma = sqrt(sigma);
  printf(" - [Std dev CPU Computation Time] %7.2f (ms)\n", sigma);

  //> Write CPU-HC number of converged solutions, real solutions, and infinity-failed solutions
  std::string write_solution_statistics_dir = std::string("../../") + WRITE_FILES_FOLDER + std::string("CPU_Sols_Statistics.txt");
  std::ofstream CPUHC_Sols_Statistics_File;
  CPUHC_Sols_Statistics_File.open(write_solution_statistics_dir);
  if ( !CPUHC_Sols_Statistics_File.is_open() ) LOG_FILE_ERROR(write_solution_statistics_dir);
  for (int i = 0; i < TEST_RANSAC_TIMES; i++) {
    CPUHC_Sols_Statistics_File << CPU_HC_.Collect_Num_Of_Coverged_Sols[i] << "\t" \
                               << CPU_HC_.Collect_Num_Of_Inf_Sols[i]      << "\t" \
                               << CPU_HC_.Collect_Num_Of_Real_Sols[i]     << "\n";
  }
  CPUHC_Sols_Statistics_File.close();

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

  //> GPU-HC
  run_GPU_HC_Solver( Problem_Settings_Map );

  //> CPU-HC
  run_CPU_HC_Solver( Problem_Settings_Map );

  return 0;
}
