#ifndef dev_eval_indxing_game6two_cuh_
#define dev_eval_indxing_game6two_cuh_
// ============================================================================
// Device function for parallel evaluations of Hx, Ht, and H of
// game6two problem
//
// Modifications
//    Chien  21-12-21:   Originally created
//
// ============================================================================
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>

// -- cuda included --
#include <cuda_runtime.h>

// -- magma included --
#include "flops.h"
#include "magma_v2.h"
#include "magma_lapack.h"
#include "magma_internal.h"
#undef max
#undef min
#include "magma_templates.h"
#include "sync.cuh"
#undef max
#undef min
#include "shuffle.cuh"
#undef max
#undef min
#include "batched_kernel_param.h"

namespace magmaHCWrapper {

    template<int N, int coefsCount>
    __device__ __inline__ void
    eval_cdt_game6two(
        const int tx, float t, magmaFloatComplex s_vector_cdt[coefsCount],
        magmaFloatComplex r_startCoefs[coefsCount], magmaFloatComplex r_targetCoefs[coefsCount])
    {
      #pragma unroll
      for (int i = 0; i < 31; i++) {
          s_vector_cdt[ tx + i * N ] = r_targetCoefs[tx + i * N] * t - r_startCoefs[tx + i * N] * (t-1);
      }
      if (tx < 5) {
        s_vector_cdt[ tx + 186 ] = r_targetCoefs[tx + 186] * t - r_startCoefs[tx + 186] * (t-1);
      }
    }

    // -------------------------------- L1-cached Hx idx and Ht idx ----------------------
    template<int N, int coefsCount, int size_Hx, int max_terms, int max_parts, int max_terms_parts, int N_max_terms_parts>
    __device__ __inline__ void
    eval_Jacobian_game6two_cached(
      const int tx,
      magmaFloatComplex *s_track, magmaFloatComplex r_cgesvA[N], 
      magmaFloatComplex s_cdt[coefsCount], const int *d_Hx_idx)
    {
      #pragma unroll
      for(int i = 0; i < N; i++) {
        r_cgesvA[i] = MAGMA_C_ZERO;

        #pragma unroll
        for(int j = 0; j < max_terms; j++) {
          r_cgesvA[i] += d_Hx_idx[j*max_parts + i*max_terms_parts + tx*N_max_terms_parts] 
                       * s_track[ d_Hx_idx[j*max_parts + 1 + i*max_terms_parts + tx*N_max_terms_parts] ] 
                       * s_track[ d_Hx_idx[j*max_parts + 2 + i*max_terms_parts + tx*N_max_terms_parts] ] 
                       * s_track[ d_Hx_idx[j*max_parts + 3 + i*max_terms_parts + tx*N_max_terms_parts] ] 
                       * s_track[ d_Hx_idx[j*max_parts + 4 + i*max_terms_parts + tx*N_max_terms_parts] ] 
                       * s_cdt[ d_Hx_idx[j*max_parts + 5 + i*max_terms_parts + tx*N_max_terms_parts] ];
        }
      }
    }

    template<int N, int size_Ht, int max_terms, int max_parts, int max_terms_parts>
    __device__ __inline__ void
    eval_Ht_game6two_cached(
      const int tx,
      magmaFloatComplex *s_track, magmaFloatComplex &r_cgesvB,
      const magmaFloatComplex* __restrict__ d_const_cd, const int* __restrict__ d_Ht_idx)
    {
      r_cgesvB = MAGMA_C_ZERO;
      #pragma unroll
      for (int i = 0; i < max_terms; i++) {
        r_cgesvB += d_Ht_idx[i*max_parts + tx*max_terms_parts] 
                  * s_track[ d_Ht_idx[i*max_parts + 1 + tx*max_terms_parts] ] 
                  * s_track[ d_Ht_idx[i*max_parts + 2 + tx*max_terms_parts] ] 
                  * s_track[ d_Ht_idx[i*max_parts + 3 + tx*max_terms_parts] ]
                  * s_track[ d_Ht_idx[i*max_parts + 4 + tx*max_terms_parts] ] 
                  * s_track[ d_Ht_idx[i*max_parts + 5 + tx*max_terms_parts] ]
                  * d_const_cd[ d_Ht_idx[i*max_parts + 6 + tx*max_terms_parts] ];
      }
    }

    template<int N, int coefsCount, int size_Ht, int max_terms, int max_parts, int max_terms_parts>
    __device__ __inline__ void
    eval_H_game6two_cached(
      const int tx,
      magmaFloatComplex *s_track, magmaFloatComplex &r_cgesvB,
      magmaFloatComplex s_cdt[coefsCount], const int* __restrict__ d_Ht_idx)
    {
      r_cgesvB = MAGMA_C_ZERO;
      #pragma unroll
      for (int i = 0; i < max_terms; i++) {
        r_cgesvB += d_Ht_idx[i*max_parts + tx*max_terms_parts] 
                  * s_track[ d_Ht_idx[i*max_parts + 1 + tx*max_terms_parts] ] 
                  * s_track[ d_Ht_idx[i*max_parts + 2 + tx*max_terms_parts] ] 
                  * s_track[ d_Ht_idx[i*max_parts + 3 + tx*max_terms_parts] ]
                  * s_track[ d_Ht_idx[i*max_parts + 4 + tx*max_terms_parts] ] 
                  * s_track[ d_Ht_idx[i*max_parts + 5 + tx*max_terms_parts] ]
                  * s_cdt[ d_Ht_idx[i*max_parts + 6 + tx*max_terms_parts] ];
      }
    }
}

#endif
