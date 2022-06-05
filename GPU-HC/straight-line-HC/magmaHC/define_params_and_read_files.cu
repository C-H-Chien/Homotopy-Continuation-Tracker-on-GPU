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

namespace magmaHCWrapper {

  void problem_params::define_problem_params(std::string HC_problem)
  {
    if (HC_problem == "alea6") {
      numOfCoeffs = 29;
      numOfVars = 6;
      numOfTracks = 387;
      
      // -- constant matrices parameters --
      Hx_maximal_terms = 4;
      Hx_maximal_parts = 4;
      Ht_maximal_terms = 6;
      Ht_maximal_parts = 5;
    }
    else if (HC_problem == "eco12") {
      numOfCoeffs = 12;
      numOfVars = 12;
      numOfTracks = 1024;

      // -- constant matrices parameters --
      Hx_maximal_terms = 11;
      Hx_maximal_parts = 4;
      Ht_maximal_terms = 12;
      Ht_maximal_parts = 5;        
    }
    else if (HC_problem == "d1") {
      numOfCoeffs = 10;
      numOfVars = 12;
      numOfTracks = 48;

      // -- constant matrices parameters --
      Hx_maximal_terms = 3;
      Hx_maximal_parts = 4;
      Ht_maximal_terms = 7;
      Ht_maximal_parts = 5;
    }
    else if (HC_problem == "cyclic7") {
      numOfCoeffs = 2;
      numOfVars = 7;
      numOfTracks = 924;

      // -- constant matrices parameters --
      Hx_maximal_terms = 6;
      Hx_maximal_parts = 8;
      Ht_maximal_terms = 7;
      Ht_maximal_parts = 9;
    }
    else if (HC_problem == "cyclic8") {
      numOfCoeffs = 2;
      numOfVars = 8;
      numOfTracks = 1152;

      // -- constant matrices parameters --
      Hx_maximal_terms = 7;
      Hx_maximal_parts = 9;
      Ht_maximal_terms = 8;
      Ht_maximal_parts = 10;
    }
    else if (HC_problem == "cyclic9") {
      numOfCoeffs = 2;
      numOfVars = 9;
      numOfTracks = 5994;

      // -- constant matrices parameters --
      Hx_maximal_terms = 8;
      Hx_maximal_parts = 10;
      Ht_maximal_terms = 9;
      Ht_maximal_parts = 11;
    }
    else if (HC_problem == "katsura6") {
      numOfCoeffs = 3;
      numOfVars = 7;
      numOfTracks = 64;

      // -- constant matrices parameters --
      Hx_maximal_terms = 3;
      Hx_maximal_parts = 3;
      Ht_maximal_terms = 8;
      Ht_maximal_parts = 4;
    }
    else if (HC_problem == "katsura7") {
      numOfCoeffs = 3;
      numOfVars = 8;
      numOfTracks = 128;

      // -- constant matrices parameters --
      Hx_maximal_terms = 3;
      Hx_maximal_parts = 3;
      Ht_maximal_terms = 9;
      Ht_maximal_parts = 4;
    }
    else if (HC_problem == "katsura8") {
      numOfCoeffs = 3;
      numOfVars = 9;
      numOfTracks = 256;

      // -- constant matrices parameters --
      Hx_maximal_terms = 3;
      Hx_maximal_parts = 3;
      Ht_maximal_terms = 10;
      Ht_maximal_parts = 4;
    }
    else if (HC_problem == "katsura9") {
      numOfCoeffs = 3;
      numOfVars = 10;
      numOfTracks = 512;

      // -- constant matrices parameters --
      Hx_maximal_terms = 3;
      Hx_maximal_parts = 3;
      Ht_maximal_terms = 11;
      Ht_maximal_parts = 4;
    }
    else if (HC_problem == "katsura10") {
      numOfCoeffs = 3;
      numOfVars = 11;
      numOfTracks = 1024;

      // -- constant matrices parameters --
      Hx_maximal_terms = 3;
      Hx_maximal_parts = 3;
      Ht_maximal_terms = 12;
      Ht_maximal_parts = 4;
    }
    else if (HC_problem == "katsura11") {
      numOfCoeffs = 3;
      numOfVars = 12;
      numOfTracks = 2048;

      // -- constant matrices parameters --
      Hx_maximal_terms = 3;
      Hx_maximal_parts = 3;
      Ht_maximal_terms = 13;
      Ht_maximal_parts = 4;
    }
    else if (HC_problem == "katsura12") {
      numOfCoeffs = 3;
      numOfVars = 13;
      numOfTracks = 4096;

      // -- constant matrices parameters --
      Hx_maximal_terms = 3;
      Hx_maximal_parts = 3;
      Ht_maximal_terms = 14;
      Ht_maximal_parts = 4;
    }
    else if (HC_problem == "katsura13") {
      numOfCoeffs = 3;
      numOfVars = 14;
      numOfTracks = 8192;

      // -- constant matrices parameters --
      Hx_maximal_terms = 3;
      Hx_maximal_parts = 3;
      Ht_maximal_terms = 15;
      Ht_maximal_parts = 4;
    }
    else if (HC_problem == "katsura14") {
      numOfCoeffs = 3;
      numOfVars = 15;
      numOfTracks = 16384;

      // -- constant matrices parameters --
      Hx_maximal_terms = 3;
      Hx_maximal_parts = 3;
      Ht_maximal_terms = 16;
      Ht_maximal_parts = 4;
    }
    else if (HC_problem == "katsura15") {
      numOfCoeffs = 3;
      numOfVars = 16;
      numOfTracks = 32768;

      // -- constant matrices parameters --
      Hx_maximal_terms = 3;
      Hx_maximal_parts = 3;
      Ht_maximal_terms = 17;
      Ht_maximal_parts = 4;
    }
    else if (HC_problem == "katsura20") {
      numOfCoeffs = 3;
      numOfVars = 21;
      //numOfTracks = 1046288;
      numOfTracks = 262144;

      // -- constant matrices parameters --
      Hx_maximal_terms = 3;
      Hx_maximal_parts = 3;
      Ht_maximal_terms = 22;
      Ht_maximal_parts = 4;
    }
    else if (HC_problem == "katsura21") {
      numOfCoeffs = 3;
      numOfVars = 22;
      numOfTracks = 262144;
      //numOfTracks = 2097152;

      // -- constant matrices parameters --
      Hx_maximal_terms = 3;
      Hx_maximal_parts = 3;
      Ht_maximal_terms = 23;
      Ht_maximal_parts = 4;
    }
    else if (HC_problem == "game6two") {
      numOfCoeffs = 191;
      numOfVars = 6;
      numOfTracks = 259;

      // -- constant matrices parameters --
      Hx_maximal_terms = 16;
      Hx_maximal_parts = 6;
      Ht_maximal_terms = 32;
      Ht_maximal_parts = 7;
    }
    else if (HC_problem == "game7two") {
      numOfCoeffs = 448;
      numOfVars = 7;
      numOfTracks = 1854;

      // -- constant matrices parameters --
      Hx_maximal_terms = 32;
      Hx_maximal_parts = 7;
      Ht_maximal_terms = 64;
      Ht_maximal_parts = 8;
    }
    else if (HC_problem == "pole28sys") {
      numOfCoeffs = 1162;
      numOfVars = 16;
      numOfTracks = 12862;

      // -- constant matrices parameters --
      Hx_maximal_terms = 8;
      Hx_maximal_parts = 3;
      Ht_maximal_terms = 73;
      Ht_maximal_parts = 4;
    }
    else {
      std::cout<<"You are entering invalid HC problem in your input argument!"<<std::endl;
      print_usage();
      exit(1);
    }
  }

  void problem_params::print_usage() {
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
    std::cerr << "./magmaHC-main -p eco12                 # solve eco12 benchmark polynomial problem\n";
    std::cerr << "===============================================================================================================\n";
  }


} // end of namespace

#endif
