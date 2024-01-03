#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <string>
// =======================================================================================================
// main function
//
// Modifications
//    Chiang-Heng Chien  21-12-29    Initially created a parametric HC for a General Computer Vision Problems
//    Chiang-Heng Chien  22-11-14    Add geometric form of 5pt relative pose problem
//    Chiang-Heng Chien  22-11-16    Add algebraic form of 5pt relative pose problem
//    Chiang-Heng Chien  23-10-20    Add gamma trick and trifocal relative pose from lines at points problem
//    Chiang-Heng Chien  23-12-27    Add definitions.hpp placing all Macros
//
//> (c) LEMS, Brown University
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// =======================================================================================================

//> magma
#include "magma_v2.h"

//> Macros
#include "magmaHC/definitions.hpp"
#include "magmaHC/GPU_HC_Solver.hpp"

//> p2c
#include "magmaHC/const-matrices/p2c-trifocal_2op1p_30x30.h"
#include "magmaHC/const-matrices/p2c-5pt_rel_pos_alg_form_quat.h"
#include "magmaHC/const-matrices/p2c-5pt_rel_pos_geo_form_quat.h"

int main(int argc, char **argv) {

  magmaFloatComplex *h_startSols;
  magmaFloatComplex *h_Track;
  magmaFloatComplex *h_startParams;
  magmaFloatComplex *h_targetParams;
  magmaFloatComplex *h_phc_coeffs_Hx;
  magmaFloatComplex *h_phc_coeffs_Ht;
  magma_int_t *h_Hx_idx;
  magma_int_t *h_Ht_idx;

  //> files to be read
  std::string repo_root_dir = REPO_PATH.append("problems/");
  std::string problem_filename = repo_root_dir.append(HC_PROBLEM);

  //> allocate tracks and coeffs arrays in cpu
  magma_cmalloc_cpu( &h_startSols,    NUM_OF_TRACKS*(NUM_OF_VARS+1) );
  magma_cmalloc_cpu( &h_Track,        NUM_OF_TRACKS*(NUM_OF_VARS+1) );
  magma_cmalloc_cpu( &h_startParams,  NUM_OF_PARAMS );
  magma_cmalloc_cpu( &h_targetParams, NUM_OF_PARAMS );

  magma_cmalloc_cpu( &h_phc_coeffs_Hx, (NUM_OF_COEFFS_FROM_PARAMS+1)*(MAX_ORDER_OF_T+1) );
  magma_cmalloc_cpu( &h_phc_coeffs_Ht, (NUM_OF_COEFFS_FROM_PARAMS+1)*(MAX_ORDER_OF_T) );
  magma_imalloc_cpu( &h_Hx_idx,        NUM_OF_VARS * NUM_OF_VARS * HX_MAXIMAL_TERMS * HX_MAXIMAL_PARTS );
  magma_imalloc_cpu( &h_Ht_idx,        NUM_OF_VARS * HT_MAXIMAL_TERMS * HT_MAXIMAL_PARTS );

  // =============================================================================
  //> read files: start solutions, start coefficients, and target parameters
  // =============================================================================
  std::string targetParam_filename_test = problem_filename;
  std::string startParams_filename_test = problem_filename;
  std::string startSols_filename_test = problem_filename;
  startSols_filename_test.append("/start_sols.txt");
  targetParam_filename_test.append("/target_params.txt");
  startParams_filename_test.append("/start_params.txt");
  std::fstream startCoef_file;
  std::fstream targetParams_file;
  std::fstream startSols_file;
  bool read_success = 0;
  bool start_sols_read_success = 0;
  bool start_coeffs_read_success = 0;
  bool targetParams_read_success = 0;
  
  float s_real, s_imag;
  int d = 0, i = 0; 
  startSols_file.open(startSols_filename_test, std::ios_base::in);
  if (!startSols_file) { std::cerr << "problem start solutions file " << startSols_filename_test << " not existed!\n"; exit(1); }
  else {
    while (startSols_file >> s_real >> s_imag) {
      (h_startSols + i * (NUM_OF_VARS+1))[d] = MAGMA_C_MAKE(s_real, s_imag);
      (h_Track + i * (NUM_OF_VARS+1))[d] = MAGMA_C_MAKE(s_real, s_imag);
      if (d < NUM_OF_VARS-1) {
        d++;
      }
      else {
        d = 0;
        i++;
      }
    }
    for(int k = 0; k < NUM_OF_TRACKS; k++) {
      (h_startSols + k * (NUM_OF_VARS+1))[NUM_OF_VARS] = MAGMA_C_MAKE(1.0, 0.0);
      (h_Track + k * (NUM_OF_VARS+1))[NUM_OF_VARS] = MAGMA_C_MAKE(1.0, 0.0);
    }
    start_sols_read_success = 1;
  }

  //> read start system coefficients
  d = 0;
  startCoef_file.open(startParams_filename_test, std::ios_base::in);
  if (!startCoef_file) {
    std::cerr << "problem start coefficients file not existed!\n";
  }
  else {
    while (startCoef_file >> s_real >> s_imag) {
      (h_startParams)[d] = MAGMA_C_MAKE(s_real, s_imag);
      d++;
    }
    start_coeffs_read_success = 1;
  }

  d = 0;
  targetParams_file.open(targetParam_filename_test, std::ios_base::in);
  if (!targetParams_file) {
    std::cerr << "problem target parameters file not existed!\n";
  }
  else {
    while (targetParams_file >> s_real >> s_imag) {
      (h_targetParams)[d] = MAGMA_C_MAKE(s_real, s_imag);
      d++;
    }
    targetParams_read_success = 1;
  }

  //>-------------------------------------------------------------------------------------------------
  bool Hx_file_read_success = false;
  bool Ht_file_read_success = false;
  
  std::string filename_Hx = problem_filename;
  std::string filename_Ht = problem_filename;
  filename_Hx.append("/Hx_idx.txt");
  filename_Ht.append("/Ht_idx.txt");
  std::fstream Hx_idx_file;
  std::fstream Ht_idx_file;
  
  //> 4) read Hx index matrix, if required
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
  //> 5) read Ht index matrix
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

  read_success = (start_sols_read_success && start_coeffs_read_success && targetParams_read_success \
                  && Hx_file_read_success && Ht_file_read_success);
  
  if (read_success) {

    std::cout << "Solving Problem " << HC_PROBLEM << std::endl;
#if REL_POSE_5PT_GEO_FORM_QUAT
    magmaHCWrapper::p2c_5pt_rel_pos_geo_form_quat(h_targetParams, h_startParams, h_phc_coeffs_Hx, h_phc_coeffs_Ht);
#elif REL_POSE_5PT_ALG_FORM_QUAT
    magmaHCWrapper::p2c_5pt_rel_pos_alg_form_quat(h_targetParams, h_startParams, h_phc_coeffs_Hx, h_phc_coeffs_Ht);
#endif

    GPU_HC_Solver GPU_HC_( h_startSols, h_Track, h_startParams, h_targetParams, h_Hx_idx, h_Ht_idx, h_phc_coeffs_Hx, h_phc_coeffs_Ht );
    
    //> Member functions
    GPU_HC_.Prepare_Files_for_Write();
    GPU_HC_.Allocate_Arrays();
    GPU_HC_.Data_Transfer_From_Host_To_Device();
    GPU_HC_.Solve_by_GPU_HC();
  }
  else {
    std::cout<<"read files failed!"<<std::endl;
    exit(1);
  }

  magma_free_cpu( h_startSols );
  magma_free_cpu( h_Track );
  magma_free_cpu( h_startParams );
  magma_free_cpu( h_targetParams );
  magma_free_cpu( h_phc_coeffs_Hx );
  magma_free_cpu( h_phc_coeffs_Ht );
  magma_free_cpu( h_Hx_idx );
  magma_free_cpu( h_Ht_idx );

  return 0;
}
