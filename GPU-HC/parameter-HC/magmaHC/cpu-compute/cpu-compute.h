#ifndef cpu_compute_h
#define cpu_compute_h
// ============================================================================
// Header file declaring all kernels
//
// Modifications
//    Chien  21-6-16:   Originally created
//
// ============================================================================
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <chrono>

#include "/gpfs/runtime/opt/openblas/0.3.7/include/cblas.h"

// -- magma --
#include "flops.h"
#include "magma_v2.h"

extern "C" {
namespace magmaHCWrapper {

  // -- cpu computation of homotopy continuation for trifocal pose estimation problem with unknown focal length --
  real_Double_t cpu_hc_solver_3view_unknownf_pHC(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_sols, magmaFloatComplex *h_Track,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    bool continue_from_gpu, magma_queue_t my_queue, int max_steps,
    magmaFloatComplex *h_track_numOfPred_s
  );

  // -- cpu computation of homotopy continuation for 4-view triangulation problem --
  real_Double_t cpu_hc_solver_4vTrg(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_sols, magmaFloatComplex *h_Track,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    bool continue_from_gpu, magma_queue_t my_queue, int max_steps,
    magmaFloatComplex *h_track_numOfPred_s
  );

  // -- cpu computation of homotopy continuation for 3-view triangulation problem --
  real_Double_t cpu_hc_solver_3vTrg(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_sols, magmaFloatComplex *h_Track,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    bool continue_from_gpu, magma_queue_t my_queue, int max_steps,
    magmaFloatComplex *h_track_numOfPred_s
  );

  // -- cpu computation of homotopy continuation for relaxed 3-view triangulation problem --
  real_Double_t cpu_hc_solver_3vTrg_relax(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_sols, magmaFloatComplex *h_Track,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    bool continue_from_gpu, magma_queue_t my_queue, int max_steps,
    magmaFloatComplex *h_track_numOfPred_s
  );

  // -- cpu computation of homotopy continuation for 5 point relative pose with depth reconstruction problem --
  real_Double_t cpu_hc_solver_5pt_rel_pose_w_depth_recon(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_sols, magmaFloatComplex *h_Track,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    bool continue_from_gpu, magma_queue_t my_queue, int max_steps,
    magmaFloatComplex *h_track_numOfPred_s
  );

  // -- cpu computation of homotopy continuation for optimal PnP with quaternion problem --
  real_Double_t cpu_hc_solver_optimalPnP_w_quaternion(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_sols, magmaFloatComplex *h_Track,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    bool continue_from_gpu, magma_queue_t my_queue, int max_steps,
    magmaFloatComplex *h_track_numOfPred_s
  );

  // -- cpu computation of homotopy continuation for 3-point relative pose with homography constraint problem --
  real_Double_t cpu_hc_solver_3pt_rel_pose_w_homo_constraint(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_sols, magmaFloatComplex *h_Track,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    bool continue_from_gpu, magma_queue_t my_queue, int max_steps,
    magmaFloatComplex *h_track_numOfPred_s
  );

  // -- cpu computation of homotopy continuation for 6-point rolling shutter absolute pose estimation problem --
  real_Double_t cpu_hc_solver_r6p(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_sols, magmaFloatComplex *h_Track,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    bool continue_from_gpu, magma_queue_t my_queue, int max_steps,
    magmaFloatComplex *h_track_numOfPred_s
  );

  real_Double_t cpu_hc_solver_refractive_p5p(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_sols, magmaFloatComplex *h_Track,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    bool continue_from_gpu, magma_queue_t my_queue, int max_steps,
    magmaFloatComplex *h_track_numOfPred_s
  );

  real_Double_t cpu_hc_solver_refractive_p6p(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_sols, magmaFloatComplex *h_Track,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    bool continue_from_gpu, magma_queue_t my_queue, int max_steps,
    magmaFloatComplex *h_track_numOfPred_s
  );
}
}

#endif
