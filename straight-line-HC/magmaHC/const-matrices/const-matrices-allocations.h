#ifndef const_matrices_allocations_h
#define const_matrices_allocations_h
// ============================================================================
// Header file declaring all kernels
//
// Modifications
//    Chien  21-05-04:   Add kernel_cgesv_batched_four_rand_fuse_reg,
//                           kernel_cgesv_batched_four_rand_mulk,
//                           kernel_get_new_data_four_rand,
//                           CPU_four_rand,
//                           compute_residual
//
// ============================================================================
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <chrono>

#include "../magmaHC-problems.cuh"

extern "C" {
namespace magmaHCWrapper {

  void const_mats::const_matrices_allocations(
    problem_params *pp, magma_int_t ldd_const_matrices_Hx_collection, magma_int_t ldd_const_matrices_Ht_collection,
    magma_int_t ldd_const_matrices_Hx, magma_int_t ldd_const_matrices_Ht,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs, magma_int_t ldd_coefs,
    magma_queue_t my_queue, std::string hc_problem
  )
  {
    int N = pp->numOfVars;
    // -------------------- CONSTANT MATRICES ---------------------
    magma_cmalloc_cpu( &h_const_cd, pp->numOfCoeffs );

    // -- constant matrices allocation in GPU --
    magma_cmalloc( &d_const_cd, ldd_coefs );

    // -- define the values of elements in the constant matrices --
    define_const_matrices( h_startCoefs, h_targetCoefs, pp->numOfCoeffs, hc_problem);

    magma_csetmatrix( pp->numOfCoeffs, 1, h_const_cd, pp->numOfCoeffs, d_const_cd, ldd_coefs, my_queue );
  }

  void const_mats::free_const_matrices()
  {
    // -- free constant matrices in CPU memory --
    magma_free_cpu( h_const_cd );

    // -- free constant matrices in GPU memory --
    magma_free( d_const_cd );
  }
}
}

#endif
