#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <string>
// ============================================================================
// main function
//
// Modifications
//    Chien  21-06-10    Revised from the original main code
//    Chien  21-11-16    Switch from CVPR paper to ISSAC paper
//
// ============================================================================
// -- nvidia cuda --
#include <cuda.h>
#include <cuda_runtime.h>

// -- magma --
#include "magma_v2.h"

// -- magma --
#include "magmaHC/magmaHC-problems.cuh"

// -- global repo directory --
std::string repo_dir = "/users/cchien3/data/cchien3/MyBitBucket/issac-benchmark-hc/";

int main(int argc, char **argv) {
  --argc; ++argv;
  std::string arg;
  int argIndx = 0;
  int argTotal = 4;
  std::string HC_problem = "default";

  // -- declare class objects (put the long lasting object in dynamic memory) --
  magmaHCWrapper::problem_params* pp = new magmaHCWrapper::problem_params;
  magmaHCWrapper::const_mats* cm = new magmaHCWrapper::const_mats;

  if (argc) {
    arg = std::string(*argv);
    if (arg == "-h" || arg == "--help") {
      pp->print_usage();
      exit(1);
    }
    else if (argc <= argTotal) {
      while(argIndx <= argTotal-1) {
        if (arg == "-p" || arg == "--problem") {
          argv++;
          arg = std::string(*argv);
          HC_problem = arg;
          argIndx+=2;
          break;
        }
        else {
          std::cerr<<"invalid input arguments! See examples: \n";
          pp->print_usage();
          exit(1);
        }
        argv++;
      }
    }
    else if (argc > argTotal) {
      std::cerr<<"too many arguments!\n";
      pp->print_usage();
      exit(1);
    }
  }
  else {
    pp->print_usage();
    exit(1);
  }

  // =============================================================================
  // -- read files: start solutions, start coefficients, and target parameters --
  // =============================================================================
  magmaFloatComplex *h_startSols;
  magmaFloatComplex *h_Track;
  magmaFloatComplex *h_startCoefs;
  magmaFloatComplex *h_targetCoefs;
  magma_int_t *h_Hx_idx;
  magma_int_t *h_Ht_idx;
  
  pp->define_problem_params(HC_problem);

  // -- files to be read and define problem filename --
  std::string problem_filename = repo_dir;
  problem_filename.append("problems/");
  problem_filename.append(HC_problem);

  // -- allocate tracks, coeffs, and Hx, Ht arrays in cpu --
  magma_cmalloc_cpu( &h_startSols, pp->numOfTracks*(pp->numOfVars+1) );
  magma_cmalloc_cpu( &h_Track, pp->numOfTracks*(pp->numOfVars+1) );
  magma_cmalloc_cpu( &h_startCoefs, pp->numOfCoeffs );
  magma_cmalloc_cpu( &h_targetCoefs, pp->numOfCoeffs );
  magma_imalloc_cpu( &h_Hx_idx, pp->numOfVars*pp->numOfVars*pp->Hx_maximal_terms*(pp->Hx_maximal_parts) );
  magma_imalloc_cpu( &h_Ht_idx, pp->numOfVars*pp->Ht_maximal_terms*(pp->Ht_maximal_parts) );
  
  std::string targetCoef_filename = problem_filename;
  std::string startCoef_filename = problem_filename;
  std::string startSols_filename = problem_filename;
  std::string filename_Hx = problem_filename;
  std::string filename_Ht = problem_filename;  
  startSols_filename.append("/start_sols.txt");
  targetCoef_filename.append("/target_coefs.txt");
  startCoef_filename.append("/start_coefs.txt");
  filename_Hx.append("/Hx_idx.txt");
  filename_Ht.append("/Ht_idx.txt");
  std::fstream startCoef_file;
  std::fstream targetParams_file;
  std::fstream startSols_file;
  std::fstream Hx_idx_file;
  std::fstream Ht_idx_file;
  bool read_success = 0;
  bool start_sols_read_success = 0;
  bool start_coeffs_read_success = 0;
  bool targetParams_read_success = 0;
  bool Hx_file_read_success = false;
  bool Ht_file_read_success = false;

  float s_real, s_imag;
  int d = 0, i = 0;
  // -- 1) read start system solutions --
  startSols_file.open(startSols_filename, std::ios_base::in);
  std::cout<<"============= start solutions h_startSols ============"<<std::endl;
  if (!startSols_file) {
    std::cerr << "problem start solutions file not existed!\n";
  }
  else {
    for (int read_sols = 0; read_sols < pp->numOfTracks*pp->numOfVars; read_sols++) {
      startSols_file >> s_real >> s_imag;
      (h_startSols + i * (pp->numOfVars+1))[d] = MAGMA_C_MAKE(s_real, s_imag);
      (h_Track + i * (pp->numOfVars+1))[d] = MAGMA_C_MAKE(s_real, s_imag);
      if (d < pp->numOfVars-1) {
        d++;
      }
      else {
        d = 0;
        i++;
      }
    }
    /*while (startSols_file >> s_real >> s_imag) {
      (h_startSols + i * (pp->numOfVars+1))[d] = MAGMA_C_MAKE(s_real, s_imag);
      (h_Track + i * (pp->numOfVars+1))[d] = MAGMA_C_MAKE(s_real, s_imag);
      if (d < pp->numOfVars-1) {
        d++;
      }
      else {
        d = 0;
        i++;
      }
    }*/
    for(int k = 0; k < pp->numOfTracks; k++) {
      (h_startSols + k * (pp->numOfVars+1))[pp->numOfVars] = MAGMA_C_MAKE(1.0, 0.0);
      (h_Track + k * (pp->numOfVars+1))[pp->numOfVars] = MAGMA_C_MAKE(1.0, 0.0);
    }
    start_sols_read_success = 1;
  }

  //magma_cprint((pp->numOfVars+1), 1, h_startSols + 28*(pp->numOfVars+1), (pp->numOfVars+1));

  d = 0;
  // -- 2) read start system coefficients --
  startCoef_file.open(startCoef_filename, std::ios_base::in);
  if (!startCoef_file) {
    std::cerr << "problem start coefficients file not existed!\n";
  }
  else {
    while (startCoef_file >> s_real >> s_imag) {
      (h_startCoefs)[d] = MAGMA_C_MAKE(s_real, s_imag);
      d++;
    }
    start_coeffs_read_success = 1;
  }

  //magma_cprint(pp->numOfCoeffs, 1, h_startCoefs, pp->numOfCoeffs);
  d = 0;
  // -- 3) read start system coefficients --
  targetParams_file.open(targetCoef_filename, std::ios_base::in);
  if (!targetParams_file) {
    std::cerr << "problem target parameters file not existed!\n";
  }
  else {
    while (targetParams_file >> s_real >> s_imag) {
      (h_targetCoefs)[d] = MAGMA_C_MAKE(s_real, s_imag);
      d++;
    }
    targetParams_read_success = 1;
  }

  // -- 4) read Hx index matrix --
  int index;
  d = 0;
  Hx_idx_file.open(filename_Hx, std::ios_base::in);
  if (!Hx_idx_file) {
    std::cerr << "problem Hx index matrix file not existed!\n";
  }
  else {
    while (Hx_idx_file >> index) {
      (h_Hx_idx)[d] = index;
      d++;
    }
    Hx_file_read_success = 1;
  }

  // -- 5) read Ht index matrix --
  d = 0;
  Ht_idx_file.open(filename_Ht, std::ios_base::in);
  if (!Ht_idx_file) {
    std::cerr << "problem Ht index matrix file not existed!\n";
  }
  else {
    while (Ht_idx_file >> index) {
      (h_Ht_idx)[d] = index;
      d++;
    }
    Ht_file_read_success = 1;
  }
  // =============================================================================

  // =============================================================================
  // -- Evaluation Files --
  // =============================================================================
  // -- files to be written --
  std::ofstream track_sols_file;
  std::ofstream target_coeffs_file;
  std::ofstream track_success_file;
  std::string write_sols_file_dir = repo_dir;
  std::string write_coeffs_file_dir = repo_dir;
  std::string write_success_file_dir = repo_dir;
  write_sols_file_dir.append("converged_HC_tracks.txt");
  write_coeffs_file_dir.append("target_coefficients.txt");
  write_success_file_dir.append("HC_number_of_successes.txt");
  track_sols_file.open(write_sols_file_dir);
  target_coeffs_file.open(write_coeffs_file_dir);
  track_success_file.open(write_success_file_dir);
  if ( !track_sols_file.is_open() || !target_coeffs_file.is_open() || !track_success_file.is_open() )
    std::cout<<"evaluation files cannot be opened!"<<std::endl;

  // -- write target coefficients to the file --
  for (int c = 0; c < pp->numOfCoeffs; c++) {
    target_coeffs_file << MAGMA_C_REAL(h_targetCoefs[c]) << "\n";
  }
  // =============================================================================

  // =============================================================================
  // -- call homotopy continuation solver --
  // =============================================================================
  // read_success = (start_sols_read_success && start_coeffs_read_success && targetParams_read_success && Hx_file_read_success && Ht_file_read_success);
  read_success = (start_sols_read_success && start_coeffs_read_success && targetParams_read_success);
  if (read_success) {
    magmaHCWrapper::homotopy_continuation_solver(h_startSols, h_Track, h_startCoefs, h_targetCoefs, h_Hx_idx, h_Ht_idx, pp, cm, HC_problem, track_sols_file, track_success_file);
  }
  else {
    std::cout<<"read files failed!"<<std::endl;
    exit(1);
  }
  // =============================================================================

  track_sols_file.close();
  target_coeffs_file.close();
  track_success_file.close();

  delete pp;
  delete cm;
  magma_free_cpu( h_startSols );
  magma_free_cpu( h_Track );
  magma_free_cpu( h_startCoefs );
  magma_free_cpu( h_targetCoefs );
  magma_free_cpu( h_Hx_idx );
  magma_free_cpu( h_Ht_idx );

  return 0;
}
