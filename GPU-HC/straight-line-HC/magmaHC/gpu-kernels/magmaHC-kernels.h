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

  // -- homotopy continuation solver - alea6 kernel --
  real_Double_t kernel_HC_Solver_alea6(
    magma_int_t N, magma_int_t batchCount, magma_int_t coefsCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_startCoefs_array, magmaFloatComplex** d_targetCoefs_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array
  );

  // -- homotopy continuation solver - echo12 kernel --
  real_Double_t kernel_HC_Solver_eco12(
    magma_int_t N, magma_int_t batchCount, magma_int_t coefsCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_startCoefs_array, magmaFloatComplex** d_targetCoefs_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array
  );

  // -- homotopy continuation solver - d1 kernel --
  real_Double_t kernel_HC_Solver_d1(
    magma_int_t N, magma_int_t batchCount, magma_int_t coefsCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_startCoefs_array, magmaFloatComplex** d_targetCoefs_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array
  );

  // -- homotopy continuation solver - cyclic7 kernel --
  real_Double_t kernel_HC_Solver_cyclic7(
    magma_int_t N, magma_int_t batchCount, magma_int_t coefsCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_startCoefs_array, magmaFloatComplex** d_targetCoefs_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array
  );

  // -- homotopy continuation solver - cyclic8 kernel --
  real_Double_t kernel_HC_Solver_cyclic8(
    magma_int_t N, magma_int_t batchCount, magma_int_t coefsCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_startCoefs_array, magmaFloatComplex** d_targetCoefs_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array
  );

  // -- homotopy continuation solver - cyclic9 kernel --
  real_Double_t kernel_HC_Solver_cyclic9(
    magma_int_t N, magma_int_t batchCount, magma_int_t coefsCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_startCoefs_array, magmaFloatComplex** d_targetCoefs_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array
  );

  real_Double_t kernel_HC_Solver_cyclic8_extractClkData(
    magma_int_t N, magma_int_t batchCount, magma_int_t coefsCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_startCoefs_array, magmaFloatComplex** d_targetCoefs_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array, long long *clocks
  );

  // -- homotopy continuation solver - katsura6 kernel --
  real_Double_t kernel_HC_Solver_katsura6(
    magma_int_t N, magma_int_t batchCount, magma_int_t coefsCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_startCoefs_array, magmaFloatComplex** d_targetCoefs_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array
  );

  // -- homotopy continuation solver - katsura9 kernel --
  real_Double_t kernel_HC_Solver_katsura7(
    magma_int_t N, magma_int_t batchCount, magma_int_t coefsCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_startCoefs_array, magmaFloatComplex** d_targetCoefs_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array
  );

  // -- homotopy continuation solver - katsura9 kernel --
  real_Double_t kernel_HC_Solver_katsura8(
    magma_int_t N, magma_int_t batchCount, magma_int_t coefsCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_startCoefs_array, magmaFloatComplex** d_targetCoefs_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array
  );

  // -- homotopy continuation solver - katsura9 kernel --
  real_Double_t kernel_HC_Solver_katsura9(
    magma_int_t N, magma_int_t batchCount, magma_int_t coefsCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_startCoefs_array, magmaFloatComplex** d_targetCoefs_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array
  );

  // -- homotopy continuation solver - katsura10 kernel --
  real_Double_t kernel_HC_Solver_katsura10(
    magma_int_t N, magma_int_t batchCount, magma_int_t coefsCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_startCoefs_array, magmaFloatComplex** d_targetCoefs_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array
  );

  // -- homotopy continuation solver - katsura11 kernel --
  real_Double_t kernel_HC_Solver_katsura11(
    magma_int_t N, magma_int_t batchCount, magma_int_t coefsCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_startCoefs_array, magmaFloatComplex** d_targetCoefs_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array
  );

  // -- homotopy continuation solver - katsura12 kernel --
  real_Double_t kernel_HC_Solver_katsura12(
    magma_int_t N, magma_int_t batchCount, magma_int_t coefsCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_startCoefs_array, magmaFloatComplex** d_targetCoefs_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array
  );

  // -- homotopy continuation solver - katsura13 kernel --
  real_Double_t kernel_HC_Solver_katsura13(
    magma_int_t N, magma_int_t batchCount, magma_int_t coefsCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_startCoefs_array, magmaFloatComplex** d_targetCoefs_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array
  );

  // -- homotopy continuation solver - katsura14 kernel --
  real_Double_t kernel_HC_Solver_katsura14(
    magma_int_t N, magma_int_t batchCount, magma_int_t coefsCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_startCoefs_array, magmaFloatComplex** d_targetCoefs_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array
  );

  // -- homotopy continuation solver - katsura15 kernel --
  real_Double_t kernel_HC_Solver_katsura15(
    magma_int_t N, magma_int_t batchCount, magma_int_t coefsCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_startCoefs_array, magmaFloatComplex** d_targetCoefs_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array
  );

  // -- homotopy continuation solver - katsura20 kernel --
  real_Double_t kernel_HC_Solver_katsura20(
    magma_int_t N, magma_int_t batchCount, magma_int_t coefsCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_startCoefs_array, magmaFloatComplex** d_targetCoefs_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array
  );

  // -- homotopy continuation solver - katsura20 kernel --
  real_Double_t kernel_HC_Solver_katsura21(
    magma_int_t N, magma_int_t batchCount, magma_int_t coefsCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_startCoefs_array, magmaFloatComplex** d_targetCoefs_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array
  );

  // -- extracting clocks within a kernel to measure partial time --
  // -- homotopy continuation solver - katsura6 kernel --
  real_Double_t kernel_HC_Solver_katsura6_extractClkData(
    magma_int_t N, magma_int_t batchCount, magma_int_t coefsCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_startCoefs_array, magmaFloatComplex** d_targetCoefs_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array, long long *clocks
  );

  // -- homotopy continuation solver - game6two kernel --
  real_Double_t kernel_HC_Solver_game6two(
    magma_int_t N, magma_int_t batchCount, magma_int_t coefsCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_startCoefs_array, magmaFloatComplex** d_targetCoefs_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array
  );

  real_Double_t kernel_HC_Solver_game7two_extractClkData(
    magma_int_t N, magma_int_t batchCount, magma_int_t coefsCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_startCoefs_array, magmaFloatComplex** d_targetCoefs_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array, long long *clocks
  );

  // -- homotopy continuation solver - game7two kernel --
  real_Double_t kernel_HC_Solver_game7two(
    magma_int_t N, magma_int_t batchCount, magma_int_t coefsCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_startCoefs_array, magmaFloatComplex** d_targetCoefs_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array
  );

  // -- homotopy continuation solver - pole28sys kernel --
  real_Double_t kernel_HC_Solver_pole28sys(
    magma_int_t N, magma_int_t batchCount, magma_int_t coefsCount, magma_int_t ldda,
    magma_queue_t my_queue,
    magmaFloatComplex** d_startSols_array, magmaFloatComplex** d_Track_array,
    magmaFloatComplex** d_startCoefs_array, magmaFloatComplex** d_targetCoefs_array,
    magmaFloatComplex** d_cgesvA_array, magmaFloatComplex** d_cgesvB_array, const_mats *cm,
    magma_int_t** d_Hx_idx_array, magma_int_t** d_Ht_idx_array
  );
}
}

#endif
