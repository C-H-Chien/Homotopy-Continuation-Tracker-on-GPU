#ifndef define_params_and_read_files_cu
#define define_params_and_read_files_cu
// ==============================================================================
//
// Modifications
//    Chiang-Heng Chien  22-10-31:   Initially Created (Copied from other repos)
//
//> (c) LEMS, Brown University
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// ==============================================================================
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <chrono>

#include "magmaHC-problems.cuh"

namespace magmaHCWrapper {

  void print_usage()
  {
    std::cerr << "===================================================================================================================\n";
    std::cerr << "Usage: ./magmaHC-main <input-argument> <command>\n\n";
    std::cerr << "Choices of input arguments and commands\n";
    std::cerr << "<input-argument>      <command>\n"
                 "       -p             <problem>                  # (or --problem)  : minimal problem name \n"
                 "       -h                                        # (or --help)     : print this help message\n\n";
    std::cerr << "----------------------- NOTICE -----------------------\n";
    std::cerr << "1. Order matters.\n";
    std::cerr << "2. If <input-argument> and <command> are not specified, the help message will be automatically shown.\n\n";
    std::cerr << "----------------------- Examples -----------------------\n";
    std::cerr << "./magmaHC-main -p six_lines_6x6      # solve six lines minimal problem under a generalized camera model\n";
    std::cerr << "===================================================================================================================\n";
  }

  void problem_params::define_problem_params(std::string problem_filename, std::string HC_problem)
  {
    if (HC_problem == "5pt_rel_pos_geo_form_quat") {
      //> problem specifications
      numOfParams = 20;
      numOfTracks = 40;
      numOfVars = 6;
      numOfCoeffsFromParams = 75;

      //> constant matrices parameters
      Hx_maximal_terms = 12;
      Hx_maximal_parts = 4;
      Ht_maximal_terms = 30;
      Ht_maximal_parts = 5;

      max_orderOf_t = 2;
    }
    else if (HC_problem == "5pt_rel_pos_alg_form_quat") {
      //> problem specifications
      numOfParams = 20;
      numOfTracks = 40;
      numOfVars = 6;
      numOfCoeffsFromParams = 75;

      //> constant matrices parameters
      Hx_maximal_terms = 12;
      Hx_maximal_parts = 4;
      Ht_maximal_terms = 30;
      Ht_maximal_parts = 5;

      max_orderOf_t = 2;
    }
    else {
      std::cout<<"You are entering invalid HC problem in your input argument!"<<std::endl;
      print_usage();
      exit(1);
    }
  }
} // end of namespace

#endif
