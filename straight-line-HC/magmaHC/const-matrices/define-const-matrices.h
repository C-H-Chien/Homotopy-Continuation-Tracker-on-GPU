#ifndef magmaHC_define_const_matrices_h
#define magmaHC_define_const_matrices_h
// ============================================================================
// Calling and defining constant matrices 
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

  // -- homotopy continuation solver - alea6 kernel --
  void const_mats::define_const_matrices(
    magmaFloatComplex *h_startCoefs, magmaFloatComplex *h_targetCoefs, int coefsCount, std::string hc_problem
  )
  {
    for(int i = 0; i < coefsCount; i++) {
      h_const_cd[i] = h_startCoefs[i] - h_targetCoefs[i];
    }
  }
}
}

#endif
