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
    magmaFloatComplex *h_startParams, magmaFloatComplex *h_targetParams, int coefsCount
  )
  {
    for(int i = 0; i < coefsCount; i++) {
      h_const_cd[i] = h_startParams[i] - h_targetParams[i];
    }
    
    // -- padding cd vector with 1 at the end of array (this is used by the brute force parametric HC method) --
    h_const_cd[coefsCount] = MAGMA_C_MAKE(1.0, 0.0);
  }
}
}

#endif
