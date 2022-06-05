#ifndef magmaHC_problems_cuh
#define magmaHC_problems_cuh
// =======================================================================
//
// Modifications
//    Chien  21-05-23:   Originally created
//    Chien  21-09-13:   Add struct of problem_params
//    Chien  21-09-14:   Add struct of
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
      int numOfCoeffs;
      int numOfVars;

      int Hx_maximal_terms;
      int Hx_maximal_parts;
      int Ht_maximal_terms;
      int Ht_maximal_parts;

      void define_problem_params(std::string HC_problem);
      void print_usage();

    private:
      std::string startSols_filename;
      std::string startCoef_filename;
      std::string targetCoef_filename;
    };

    // -- structs for constant matrices --
    class const_mats {
    public:

      // -- variables arrays --
      // -- Hx has to declare maximal possible Xs --
      // -- 1) CPU --
      magmaFloatComplex *h_const_cd;

      // -- 2) GPU --
      magmaFloatComplex_ptr d_const_cd;

      // -- member functions --
      void const_matrices_allocations(
        problem_params *pp, magma_int_t ldd_const_matrices_Hx_collection, magma_int_t ldd_const_matrices_Ht_collection,
        magma_int_t ldd_const_matrices_Hx, magma_int_t ldd_const_matrices_Ht,
        magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs, magma_int_t ldd_coefs,
        magma_queue_t my_queue, std::string hc_problem
      );

      void define_const_matrices(
        magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs, int coefsCount, std::string hc_problem
      );

      void free_const_matrices();
    };

    void homotopy_continuation_solver(
      magmaFloatComplex *h_startSols, magmaFloatComplex *h_Track,
      magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
      magma_int_t *h_Hx_idx, magma_int_t *h_Ht_idx,
      problem_params* pp, const_mats *cm, std::string hc_problem,
      std::ofstream &track_sols_file, std::ofstream &tracks_success_file
    );
}

#endif
