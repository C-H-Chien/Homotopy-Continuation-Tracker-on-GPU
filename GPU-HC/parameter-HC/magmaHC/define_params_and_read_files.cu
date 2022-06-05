#ifndef define_params_and_read_files_cu
#define define_params_and_read_files_cu
// =======================================================================
//
// Modifications
//    Chien  21-09-13:   Originally Created
//
// =======================================================================
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <chrono>

#include "magmaHC-problems.cuh"
//#include "magmaHC/input-info.h"

namespace magmaHCWrapper {

  void print_usage()
  {
    std::cerr << "===============================================================================================================\n";
    std::cerr << "Usage: ./magmaHC-main <input-argument> <command>\n\n";
    std::cerr << "Choices of input arguments and commands\n";
    std::cerr << "<input-argument>      <command>\n"
                 "       -p             <problem>         # (or --problem)  : minimal problem name \n"
                 "       -h                               # (or --help)     : print this help message\n\n";
    std::cerr << "----------------------- NOTICE -----------------------\n";
    std::cerr << "1. Order matters.\n";
    std::cerr << "2. If <input-argument> and <command> are not specified, the help message will be automatically shown.\n\n";
    std::cerr << "----------------------- Examples -----------------------\n";
    std::cerr << "./magmaHC-main -p 3vTrg                 # solve 3-view triangulation problem\n";
    std::cerr << "===============================================================================================================\n";
  }

  void problem_params::define_problem_params(std::string problem_filename, std::string HC_problem)
  {
    if (HC_problem == "3view_unknownf_pHC") {
      // -- problem specifications --
      numOfParams = 36;
      numOfTracks = 1784;
      numOfVars = 18;
      numOfCoeffsFromParams = 153;

      // -- constant matrices parameters --
      Hx_maximal_terms = 12;
      Hx_maximal_parts = 6;
      Ht_maximal_terms = 14;
      Ht_maximal_parts = 7;

      max_orderOf_t = 2;
    }
    else if (HC_problem == "3vTrg") {
      // -- problem specifications --
      numOfParams = 33;
      numOfTracks = 94;
      numOfVars = 9;
      numOfCoeffsFromParams = 34;

      // -- constant matrices parameters --
      Hx_maximal_terms = 3;
      Hx_maximal_parts = 3;
      Ht_maximal_terms = 9;
      Ht_maximal_parts = 4;

      max_orderOf_t = 1;
    }
    else if (HC_problem == "3vTrg_relax") {
      // -- problem specifications --
      numOfParams = 24;
      numOfTracks = 31;
      numOfVars = 8;
      numOfCoeffsFromParams = 25;

      // -- constant matrices parameters --
      Hx_maximal_terms = 3;
      Hx_maximal_parts = 3;
      Ht_maximal_terms = 9;
      Ht_maximal_parts = 4;

      max_orderOf_t = 1;
    }
    else if (HC_problem == "4vTrg") {
      // -- problem specifications --
      numOfParams = 35;
      numOfTracks = 142;
      numOfVars = 11;
      numOfCoeffsFromParams = 36;

      // -- constant matrices parameters --
      Hx_maximal_terms = 3;
      Hx_maximal_parts = 3;
      Ht_maximal_terms = 9;
      Ht_maximal_parts = 4;

      max_orderOf_t = 1;
    }
    else if (HC_problem == "5pt_rel_pose_w_depth_recon") {
      // -- problem specifications --
      numOfParams = 20;
      numOfTracks = 40;
      numOfVars = 16;
      numOfCoeffsFromParams = 44;

      // -- constant matrices parameters --
      Hx_maximal_terms = 7;
      Hx_maximal_parts = 4;
      Ht_maximal_terms = 9;
      Ht_maximal_parts = 5;

      max_orderOf_t = 1;
    }
    else if (HC_problem == "optimalPnP_w_quaternion") {
      // -- problem specifications --
      numOfParams = 81;
      numOfTracks = 80;
      numOfVars = 4;
      numOfCoeffsFromParams = 107;

      // -- constant matrices parameters --
      Hx_maximal_terms = 20;
      Hx_maximal_parts = 5;
      Ht_maximal_terms = 35;
      Ht_maximal_parts = 6;

      max_orderOf_t = 1;
    }
    else if (HC_problem == "3pt_rel_pose_w_homo_constraint") {
      // -- problem specifications --
      numOfParams = 18;
      numOfTracks = 8;
      numOfVars = 8;
      numOfCoeffsFromParams = 44;

      // -- constant matrices parameters --
      Hx_maximal_terms = 6;
      Hx_maximal_parts = 4;
      Ht_maximal_terms = 9;
      Ht_maximal_parts = 5;

      max_orderOf_t = 2;
    }
    else if (HC_problem == "r6p") {
      // -- problem specifications --
      numOfParams = 111;
      numOfTracks = 20;
      numOfVars = 6;
      numOfCoeffsFromParams = 96;

      // -- constant matrices parameters --
      Hx_maximal_terms = 4;
      Hx_maximal_parts = 3;
      Ht_maximal_terms = 16;
      Ht_maximal_parts = 4;

      max_orderOf_t = 3;
    }
    else if (HC_problem == "refractive_p5p") {
      // -- problem specifications --
      numOfParams = 29;
      numOfTracks = 16;
      numOfVars = 5;
      numOfCoeffsFromParams = 90;

      // -- constant matrices parameters --
      Hx_maximal_terms = 10;
      Hx_maximal_parts = 4;
      Ht_maximal_terms = 26;
      Ht_maximal_parts = 5;

      max_orderOf_t = 3;
    }
    else if (HC_problem == "refractive_p6p") {
      // -- problem specifications --
      numOfParams = 35;
      numOfTracks = 36;
      numOfVars = 6;
      numOfCoeffsFromParams = 108;

      // -- constant matrices parameters --
      Hx_maximal_terms = 10;
      Hx_maximal_parts = 5;
      Ht_maximal_terms = 26;
      Ht_maximal_parts = 6;

      max_orderOf_t = 3;
    }
    else {
      std::cout<<"You are entering invalid HC problem in your input argument!"<<std::endl;
      print_usage();
      exit(1);
    }
  }
} // end of namespace

#endif
