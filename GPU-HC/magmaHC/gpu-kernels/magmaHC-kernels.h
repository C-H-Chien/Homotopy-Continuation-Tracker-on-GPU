#ifndef magmaHC_kernels_h
#define magmaHC_kernels_h
// ============================================================================
// Header file declaring all kernels
//
// Modifications
//    Chiang-Heng Chien  22-10-31:   Initially Created (Copied from other repos)
//
//> (c) LEMS, Brown University
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// ============================================================================
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <chrono>

#include "../magmaHC-problems.cuh"

// -- magma --
#include "flops.h"
#include "magma_v2.h"

extern "C" {
namespace magmaHCWrapper {

  //> Geometric Form
  real_Double_t kernel_HC_Solver_5pt_rel_pos_geo_form_quat(
    magma_int_t N, magma_int_t batchCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t* d_Hx_idx_array, magma_int_t* d_Ht_idx_array,
    magmaFloatComplex_ptr d_phc_coeffs_Hx, magmaFloatComplex_ptr d_phc_coeffs_Ht,
    magma_int_t numOf_phc_coeffs
  );

  //> Algebraic Form
  real_Double_t kernel_HC_Solver_5pt_rel_pos_alg_form_quat(
    magma_int_t N, magma_int_t batchCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t* d_Hx_idx_array, magma_int_t* d_Ht_idx_array,
    magmaFloatComplex_ptr d_phc_coeffs_Hx, magmaFloatComplex_ptr d_phc_coeffs_Ht,
    magma_int_t numOf_phc_coeffs
  );

}
}

#endif
