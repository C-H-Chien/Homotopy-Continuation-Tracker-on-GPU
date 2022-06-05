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
//    Chien  21-12-29    initial developments
//
// ============================================================================
// -- nvidia cuda --
#include <cuda.h>
#include <cuda_runtime.h>

// -- magma --
#include "magma_v2.h"

// -- input cmd --
//#include "magmaHC/input-info.h"

// -- magma --
#include "magmaHC/magmaHC-problems.cuh"

// -- p2c --
#include "magmaHC/const-matrices/p2c-3view_unknownf_pHC.h"
#include "magmaHC/const-matrices/p2c-4vTrg.h"
#include "magmaHC/const-matrices/p2c-3vTrg.h"
#include "magmaHC/const-matrices/p2c-3vTrg_relax.h"
#include "magmaHC/const-matrices/p2c-5pt_rel_pose_w_depth_recon.h"
#include "magmaHC/const-matrices/p2c-optimalPnP_w_quaternion.h"
#include "magmaHC/const-matrices/p2c-3pt_rel_pose_w_homo_constraint.h"
#include "magmaHC/const-matrices/p2c-r6p.h"
#include "magmaHC/const-matrices/p2c-refractive_p5p.h"
#include "magmaHC/const-matrices/p2c-refractive_p6p.h"

// -- global repo directory --
std::string repo_dir = "/users/cchien3/data/cchien3/MyBitBucket/issac-parametric-hc/";

int main(int argc, char **argv) {
  --argc; ++argv;
  std::string arg;
  int argIndx = 0;
  int argTotal = 4;
  std::string HC_problem = "default";

  if (argc) {
    arg = std::string(*argv);
    if (arg == "-h" || arg == "--help") {
      magmaHCWrapper::print_usage();
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
          magmaHCWrapper::print_usage();
          exit(1);
        }
        argv++;
      }
    }
    else if (argc > argTotal) {
      std::cerr<<"too many arguments!\n";
      magmaHCWrapper::print_usage();
      exit(1);
    }
  }
  else {
    magmaHCWrapper::print_usage();
    exit(1);
  }

  magmaFloatComplex *h_startSols;
  magmaFloatComplex *h_Track;
  magmaFloatComplex *h_startParams;
  magmaFloatComplex *h_targetParams;
  magmaFloatComplex *h_phc_coeffs_Hx;
  magmaFloatComplex *h_phc_coeffs_Ht;
  magma_int_t *h_Hx_idx;
  magma_int_t *h_Ht_idx;

  // -- files to be read --
  std::string repo_root_dir = repo_dir;
  repo_dir.append("problems/");
  std::string problem_filename = repo_dir.append(HC_problem);

  // -- declare class objects (put the long lasting object in dynamic memory) --
  magmaHCWrapper::problem_params* pp = new magmaHCWrapper::problem_params;
  magmaHCWrapper::const_mats* cm = new magmaHCWrapper::const_mats;

  pp->define_problem_params(problem_filename, HC_problem);

  // -- allocate tracks and coeffs arrays in cpu --
  magma_cmalloc_cpu( &h_startSols, pp->numOfTracks*(pp->numOfVars+1) );
  magma_cmalloc_cpu( &h_Track, pp->numOfTracks*(pp->numOfVars+1) );
  magma_cmalloc_cpu( &h_startParams, pp->numOfParams );
  magma_cmalloc_cpu( &h_targetParams, pp->numOfParams );

  magma_cmalloc_cpu( &h_phc_coeffs_Hx, pp->numOfCoeffsFromParams*(pp->max_orderOf_t+1) );
  magma_cmalloc_cpu( &h_phc_coeffs_Ht, pp->numOfCoeffsFromParams*(pp->max_orderOf_t) );
  magma_imalloc_cpu( &h_Hx_idx, pp->numOfVars*pp->numOfVars*pp->Hx_maximal_terms*pp->Hx_maximal_parts );
  magma_imalloc_cpu( &h_Ht_idx, pp->numOfVars*pp->Ht_maximal_terms*pp->Ht_maximal_parts );

  // =============================================================================
  // -- read files: start solutions, start coefficients, and target parameters --
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
  std::cout<<"============= start solutions h_startSols ============"<<std::endl;
  if (!startSols_file) {
    std::cerr << "problem start solutions file not existed!\n";
  }
  else {
    while (startSols_file >> s_real >> s_imag) {
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
    for(int k = 0; k < pp->numOfTracks; k++) {
      (h_startSols + k * (pp->numOfVars+1))[pp->numOfVars] = MAGMA_C_MAKE(1.0, 0.0);
      (h_Track + k * (pp->numOfVars+1))[pp->numOfVars] = MAGMA_C_MAKE(1.0, 0.0);
    }
    start_sols_read_success = 1;
  }

  // -- read start system coefficients --
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

  // ---------------------------------------------------------------------------------------------------
  bool Hx_file_read_success = false;
  bool Ht_file_read_success = false;
  
  std::string filename_Hx = problem_filename;
  std::string filename_Ht = problem_filename;
  filename_Hx.append("/Hx_idx.txt");
  filename_Ht.append("/Ht_idx.txt");
  std::fstream Hx_idx_file;
  std::fstream Ht_idx_file;
  
  // -- 4) read Hx index matrix, if required --
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

  magma_cprint(pp->numOfParams, 1, h_targetParams, pp->numOfParams);

  // -- params2coeffs --
  // -- coefficients of function t for H, used in gpu --
  if(HC_problem == "3view_unknownf_pHC")
    magmaHCWrapper::p2c_3view_unknwonf_pHC(h_targetParams, h_startParams, h_phc_coeffs_Hx, h_phc_coeffs_Ht);
  else if (HC_problem == "4vTrg")
    magmaHCWrapper::p2c_4vTrg(h_targetParams, h_startParams, h_phc_coeffs_Hx, h_phc_coeffs_Ht);
  else if (HC_problem == "3vTrg")
    magmaHCWrapper::p2c_3vTrg(h_targetParams, h_startParams, h_phc_coeffs_Hx, h_phc_coeffs_Ht);
  else if (HC_problem == "3vTrg_relax")
    magmaHCWrapper::p2c_3vTrg_relax(h_targetParams, h_startParams, h_phc_coeffs_Hx, h_phc_coeffs_Ht);
  else if (HC_problem == "5pt_rel_pose_w_depth_recon")
    magmaHCWrapper::p2c_5pt_rel_pose_w_depth_recon(h_targetParams, h_startParams, h_phc_coeffs_Hx, h_phc_coeffs_Ht);
  else if (HC_problem == "optimalPnP_w_quaternion")
    magmaHCWrapper::p2c_optimalPnP_w_quaternion(h_targetParams, h_startParams, h_phc_coeffs_Hx, h_phc_coeffs_Ht);
  else if (HC_problem == "3pt_rel_pose_w_homo_constraint")
    magmaHCWrapper::p2c_3pt_rel_pose_w_homo_constraint(h_targetParams, h_startParams, h_phc_coeffs_Hx, h_phc_coeffs_Ht);
  else if (HC_problem == "r6p")
    magmaHCWrapper::p2c_r6p(h_targetParams, h_startParams, h_phc_coeffs_Hx, h_phc_coeffs_Ht);
  else if (HC_problem == "refractive_p5p")
    magmaHCWrapper::p2c_refractive_p5p(h_targetParams, h_startParams, h_phc_coeffs_Hx, h_phc_coeffs_Ht);
  else if (HC_problem == "refractive_p6p")
    magmaHCWrapper::p2c_refractive_p6p(h_targetParams, h_startParams, h_phc_coeffs_Hx, h_phc_coeffs_Ht);

  read_success = (start_sols_read_success && start_coeffs_read_success && targetParams_read_success && Hx_file_read_success && Ht_file_read_success);

  // -- write successful track solutions to files --
  std::ofstream track_sols_file;
  std::string write_sols_file_dir = repo_root_dir;
  write_sols_file_dir.append("converged_HC_tracks.txt");
  track_sols_file.open(write_sols_file_dir);
  if ( !track_sols_file.is_open() )
    std::cout<<"successful solutions for evaluation write files cannot be opened!"<<std::endl;

  // -- call homotopy continuation solver --
  //read_success = (start_sols_read_success && start_coeffs_read_success && targetParams_read_success);
  if (read_success) {
    magmaHCWrapper::homotopy_continuation_solver(h_startSols, h_Track, h_startParams, h_targetParams, h_Hx_idx, h_Ht_idx, h_phc_coeffs_Hx, h_phc_coeffs_Ht, pp, cm, HC_problem, track_sols_file);
  }
  else {
    std::cout<<"read files failed!"<<std::endl;
    exit(1);
  }

  delete pp;
  delete cm;
  magma_free_cpu( h_startSols );
  magma_free_cpu( h_Track );
  magma_free_cpu( h_startParams );
  magma_free_cpu( h_targetParams );
  magma_free_cpu( h_phc_coeffs_Hx );
  magma_free_cpu( h_phc_coeffs_Ht );

  magma_free_cpu( h_Hx_idx );
  magma_free_cpu( h_Ht_idx );

  track_sols_file.close();

  return 0;
}
