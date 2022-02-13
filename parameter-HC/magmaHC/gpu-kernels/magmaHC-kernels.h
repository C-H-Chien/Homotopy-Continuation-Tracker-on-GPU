#ifndef magmaHC_kernels_h
#define magmaHC_kernels_h
// ============================================================================
// Header file declaring all kernels
//
// Modifications
//    Chien  21-05-04:
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

// -- magma --
#include "flops.h"
#include "magma_v2.h"

extern "C" {
namespace magmaHCWrapper {
  
  // -- GPU kernels of parametric HC --

  real_Double_t kernel_HC_Solver_3view_unknownf_pHC(
    magma_int_t N, magma_int_t batchCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_startCoefs_array, magmaFloatComplex** d_targetCoefs_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array,
    magmaFloatComplex_ptr d_phc_coeffs_H, magmaFloatComplex_ptr d_phc_coeffs_Ht,
    magma_int_t numOf_phc_coeffs
  );

  real_Double_t kernel_HC_Solver_4vTrg(
    magma_int_t N, magma_int_t batchCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array,
    magmaFloatComplex_ptr d_phc_coeffs_Hx, magmaFloatComplex_ptr d_phc_coeffs_Ht,
    magma_int_t numOf_phc_coeffs
  );

  real_Double_t kernel_HC_Solver_3vTrg(
    magma_int_t N, magma_int_t batchCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array,
    magmaFloatComplex_ptr d_phc_coeffs_Hx, magmaFloatComplex_ptr d_phc_coeffs_Ht,
    magma_int_t numOf_phc_coeffs
  );

  real_Double_t kernel_HC_Solver_3vTrg_relax(
    magma_int_t N, magma_int_t batchCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array,
    magmaFloatComplex_ptr d_phc_coeffs_Hx, magmaFloatComplex_ptr d_phc_coeffs_Ht,
    magma_int_t numOf_phc_coeffs
  );

  real_Double_t kernel_HC_Solver_5pt_rel_pose_w_depth_recon(
    magma_int_t N, magma_int_t batchCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array,
    magmaFloatComplex_ptr d_phc_coeffs_Hx, magmaFloatComplex_ptr d_phc_coeffs_Ht,
    magma_int_t numOf_phc_coeffs
  );

  real_Double_t kernel_HC_Solver_PnP_wo_principal_point(
    magma_int_t N, magma_int_t batchCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array,
    magmaFloatComplex_ptr d_phc_coeffs_Hx, magmaFloatComplex_ptr d_phc_coeffs_Ht,
    magma_int_t numOf_phc_coeffs
  );

  real_Double_t kernel_HC_Solver_optimalPnP_w_quaternion(
    magma_int_t N, magma_int_t batchCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array,
    magmaFloatComplex_ptr d_phc_coeffs_Hx, magmaFloatComplex_ptr d_phc_coeffs_Ht,
    magma_int_t numOf_phc_coeffs
  );
  
  real_Double_t kernel_HC_Solver_3pt_rel_pose_w_homo_constraint(
    magma_int_t N, magma_int_t batchCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array,
    magmaFloatComplex_ptr d_phc_coeffs_Hx, magmaFloatComplex_ptr d_phc_coeffs_Ht,
    magma_int_t numOf_phc_coeffs
  );

  real_Double_t kernel_HC_Solver_r6p(
    magma_int_t N, magma_int_t batchCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array,
    magmaFloatComplex_ptr d_phc_coeffs_Hx, magmaFloatComplex_ptr d_phc_coeffs_Ht,
    magma_int_t numOf_phc_coeffs
  );

  real_Double_t kernel_HC_Solver_refractive_p5p(
    magma_int_t N, magma_int_t batchCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array,
    magmaFloatComplex_ptr d_phc_coeffs_Hx, magmaFloatComplex_ptr d_phc_coeffs_Ht,
    magma_int_t numOf_phc_coeffs
  );

  real_Double_t kernel_HC_Solver_refractive_p6p(
    magma_int_t N, magma_int_t batchCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array,
    magmaFloatComplex_ptr d_phc_coeffs_Hx, magmaFloatComplex_ptr d_phc_coeffs_Ht,
    magma_int_t numOf_phc_coeffs
  );
}
}

#endif
