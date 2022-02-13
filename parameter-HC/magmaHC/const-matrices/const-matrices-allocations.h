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
    problem_params *pp, 
    magmaFloatComplex *h_startParams, magmaFloatComplex *h_targetParams,
    magma_queue_t my_queue
  )
  {
    int N = pp->numOfVars;
    // -------------------- CONSTANT MATRICES ---------------------
    magma_cmalloc_cpu( &h_const_cd, pp->numOfParams+1 );

    magma_int_t ldd_padded_coeffs = magma_roundup( pp->numOfParams+1, 32 );
    magma_cmalloc( &d_const_cd, ldd_padded_coeffs );

    // -- define the values of elements in the constant matrices --
    define_const_matrices( h_startParams, h_targetParams, pp->numOfParams);

    magma_csetmatrix( pp->numOfParams+1, 1, h_const_cd, pp->numOfParams+1, d_const_cd, ldd_padded_coeffs, my_queue );
  }

  void const_mats::free_const_matrices()
  {
    magma_free_cpu( h_const_cd );
    magma_free( d_const_cd );
  }
}
}

#endif
