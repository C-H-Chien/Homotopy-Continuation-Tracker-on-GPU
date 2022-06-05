#ifndef magmaHC_problems_cuh
#define magmaHC_problems_cuh
// =======================================================================
//
// Modifications
//    Chien  21-10-12:   Intiailly Created
//
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

    // -- structs for constant matrices --
    class const_mats {
    public:
      magmaFloatComplex *h_const_cd;

      // -- Hx --
      magma_int_t *h_Hx;
      magma_int_t *h_Ht;

      magmaFloatComplex_ptr d_const_cd;

      // -- member functions --
      void const_matrices_allocations(
        problem_params *pp, 
        magmaFloatComplex *h_startParams, magmaFloatComplex *h_targetParams,
        magma_queue_t my_queue
      );

      void define_const_matrices(
        magmaFloatComplex *h_startParams, magmaFloatComplex *h_targetParams, int coefsCount
      );

      void free_const_matrices();
    };

    void print_usage();

    void homotopy_continuation_solver(
      magmaFloatComplex *h_startSols, magmaFloatComplex *h_Track,
      magmaFloatComplex *h_startParams, magmaFloatComplex *h_targetParams,
      magma_int_t *h_Hx_idx, magma_int_t *h_Ht_idx,
      magmaFloatComplex *h_phc_coeffs_H, magmaFloatComplex *h_phc_coeffs_Ht,
      problem_params* pp, const_mats *cm, std::string hc_problem, std::ofstream &track_sols_file
    );
}

#endif
