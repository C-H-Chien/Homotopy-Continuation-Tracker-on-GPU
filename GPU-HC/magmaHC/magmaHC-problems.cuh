#ifndef magmaHC_problems_hpp
#define magmaHC_problems_hpp
// =======================================================================
//
// Modifications
//    Chiang-Heng Chien  21-10-12:   Intiailly Created
//
// Notes
//    Chiang-Heng Chien  22-11-12:   (TODO) This script has to be reorganized.
//
//> (c) LEMS, Brown University
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// =======================================================================
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>
#include "magma_v2.h"

namespace magmaHCWrapper {

    // -- informaiton of the benchmark problem --
    class problem_params {
    public:
      int numOfTracks;
      int numOfParams;
      int numOfVars;
      int numOfCoeffsFromParams;

      int Hx_maximal_terms;
      int Hx_maximal_parts;
      int Ht_maximal_terms;
      int Ht_maximal_parts;

      int max_orderOf_t;

      void define_problem_params(std::string problem_filename, std::string HC_problem);

    private:
      std::string startSols_filename;
      std::string startCoef_filename;
      std::string targetParam_filename;
    };

    void print_usage();

    void homotopy_continuation_solver(
      magmaFloatComplex *h_startSols, magmaFloatComplex *h_Track,
      magmaFloatComplex *h_startParams, magmaFloatComplex *h_targetParams,
      magma_int_t *h_Hx_idx, magma_int_t *h_Ht_idx,
      magmaFloatComplex *h_phc_coeffs_H, magmaFloatComplex *h_phc_coeffs_Ht,
      problem_params* pp, std::string hc_problem, std::ofstream &track_sols_file
    );
}

#endif
