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

  // -- cpu computation of homotopy continuation for alea6 problem --
  real_Double_t cpu_hc_solver_alea6(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_Track_gpu,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    magma_queue_t my_queue, int max_steps, std::ofstream &tracks_success_file
  );

  // -- cpu computation of homotopy continuation for eco12 problem --
  real_Double_t cpu_hc_solver_echo12(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_Track_gpu,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    magma_queue_t my_queue, int max_steps, std::ofstream &tracks_success_file
  );

  // -- cpu computation of homotopy continuation for d1 problem --
  real_Double_t cpu_hc_solver_d1(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_Track_gpu,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    magma_queue_t my_queue, int max_steps, std::ofstream &tracks_success_file
  );

  // -- cpu computation of homotopy continuation for cyclic7 problem --
  real_Double_t cpu_hc_solver_cyclic7(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_Track_gpu,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    magma_queue_t my_queue, int max_steps, std::ofstream &tracks_success_file
  );

  // -- cpu computation of homotopy continuation for cyclic8 problem --
  real_Double_t cpu_hc_solver_cyclic8(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_Track_gpu,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    magma_queue_t my_queue, int max_steps, std::ofstream &tracks_success_file
  );

  // -- cpu computation of homotopy continuation for cyclic9 problem --
  real_Double_t cpu_hc_solver_cyclic9(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_Track_gpu,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    magma_queue_t my_queue, int max_steps, std::ofstream &tracks_success_file
  );

  // -- cpu computation of homotopy continuation for katsura6 problem --
  real_Double_t cpu_hc_solver_katsura6(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_Track_gpu,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    magma_queue_t my_queue, int max_steps, std::ofstream &tracks_success_file
  );

  // -- cpu computation of homotopy continuation for katsura8 problem --
  real_Double_t cpu_hc_solver_katsura7(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_Track_gpu,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    magma_queue_t my_queue, int max_steps, std::ofstream &tracks_success_file
  );

  // -- cpu computation of homotopy continuation for katsura8 problem --
  real_Double_t cpu_hc_solver_katsura8(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_Track_gpu,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    magma_queue_t my_queue, int max_steps, std::ofstream &tracks_success_file
  );

  // -- cpu computation of homotopy continuation for katsura9 problem --
  real_Double_t cpu_hc_solver_katsura9(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_Track_gpu,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    magma_queue_t my_queue, int max_steps, std::ofstream &tracks_success_file
  );

  // -- cpu computation of homotopy continuation for katsura10 problem --
  real_Double_t cpu_hc_solver_katsura10(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_Track_gpu,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    magma_queue_t my_queue, int max_steps, std::ofstream &tracks_success_file
  );

  // -- cpu computation of homotopy continuation for katsura11 problem --
  real_Double_t cpu_hc_solver_katsura11(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_Track_gpu,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    magma_queue_t my_queue, int max_steps, std::ofstream &tracks_success_file
  );

  // -- cpu computation of homotopy continuation for katsura12 problem --
  real_Double_t cpu_hc_solver_katsura12(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_Track_gpu,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    magma_queue_t my_queue, int max_steps, std::ofstream &tracks_success_file
  );

  // -- cpu computation of homotopy continuation for katsura13 problem --
  real_Double_t cpu_hc_solver_katsura13(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_Track_gpu,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    magma_queue_t my_queue, int max_steps, std::ofstream &tracks_success_file
  );

  // -- cpu computation of homotopy continuation for katsura14 problem --
  real_Double_t cpu_hc_solver_katsura14(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_Track_gpu,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    magma_queue_t my_queue, int max_steps, std::ofstream &tracks_success_file
  );

  // -- cpu computation of homotopy continuation for katsura15 problem --
  real_Double_t cpu_hc_solver_katsura15(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_Track_gpu,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    magma_queue_t my_queue, int max_steps, std::ofstream &tracks_success_file
  );

  // -- cpu computation of homotopy continuation for katsura20 problem --
  real_Double_t cpu_hc_solver_katsura20(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_Track_gpu,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    magma_queue_t my_queue, int max_steps, std::ofstream &tracks_success_file
  );

  // -- cpu computation of homotopy continuation for katsura21 problem --
  real_Double_t cpu_hc_solver_katsura21(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_Track_gpu,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    magma_queue_t my_queue, int max_steps, std::ofstream &tracks_success_file
  );

  // -- cpu computation of homotopy continuation for game6two problem --
  real_Double_t cpu_hc_solver_game6two(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_Track_gpu,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    magma_queue_t my_queue, int max_steps, std::ofstream &tracks_success_file
  );

  // -- cpu computation of homotopy continuation for game7two problem --
  real_Double_t cpu_hc_solver_game7two(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_Track_gpu,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    magma_queue_t my_queue, int max_steps, std::ofstream &tracks_success_file
  );

  // -- cpu computation of homotopy continuation for pole28sys problem --
  real_Double_t cpu_hc_solver_pole28sys(
    magmaFloatComplex *h_Track_cpu, magma_int_t *h_Track_Success,
    magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_Track_gpu,
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs,
    magmaFloatComplex *h_cgesvA, magmaFloatComplex *h_cgesvB,
    magma_int_t batchCount, magma_int_t coefsCount, magma_int_t N,
    magma_queue_t my_queue, int max_steps, std::ofstream &tracks_success_file
  );
}
}

#endif
