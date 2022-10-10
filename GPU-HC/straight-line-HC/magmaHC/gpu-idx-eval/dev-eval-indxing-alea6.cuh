#ifndef dev_eval_indxing_alea6_cuh_
#define dev_eval_indxing_alea6_cuh_
// ============================================================================
// Device function for parallel evaluations of Hx, Ht, and H of
// alea6 problem
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
    eval_cdt_alea6(
        const int tx, float t, magmaFloatComplex s_vector_cdt[coefsCount],
        magmaFloatComplex r_startCoefs[coefsCount], magmaFloatComplex r_targetCoefs[coefsCount])
    {
      #pragma unroll
      for (int i = 0; i < 4; i++) {
          s_vector_cdt[ tx + i * N ] = r_targetCoefs[tx + i * N] * t - r_startCoefs[tx + i * N] * (t-1);
      }
      if (tx < 5) {
        s_vector_cdt[ tx + 24 ] = r_targetCoefs[tx + 24] * t - r_startCoefs[tx + 24] * (t-1);
      }
    }

    // -- Improved Hx parallel indexing --
    template<int N, int coefsCount, int size_Hx, int max_terms, int max_parts>
    __device__ __inline__ void
    eval_Hx_alea6(
        magmaFloatComplex *s_track, magmaFloatComplex r_cgesvA[N], 
        magmaFloatComplex s_cdt[coefsCount], magma_int_t r_Hx_idx[size_Hx])
    {
      #pragma unroll
      for(int i = 0; i < N; i++) {
        r_cgesvA[i] = MAGMA_C_ZERO;

        #pragma unroll
        for(int j = 0; j < max_terms; j++) {
          r_cgesvA[i] += r_Hx_idx[j*max_parts + i*max_terms*max_parts] 
                       * s_track[ r_Hx_idx[j*max_parts + 1 + i*max_terms*max_parts] ]
                       * s_track[ r_Hx_idx[j*max_parts + 2 + i*max_terms*max_parts] ]
                       * s_cdt[ r_Hx_idx[j*max_parts + 3 + i*max_terms*max_parts] ];
        }
      }
    }

    // -- Improved Ht parallel indexing --
    template<int N, int size_Ht, int max_terms, int max_parts>
    __device__ __inline__ void
    eval_Ht_alea6(
        magmaFloatComplex *s_track, magmaFloatComplex &r_cgesvB,
        const magmaFloatComplex* __restrict__ d_const_cd, magma_int_t r_Ht_idx[size_Ht])
    {
      r_cgesvB = MAGMA_C_ZERO;
      #pragma unroll
      for (int i = 0; i < max_terms; i++) {
        r_cgesvB += r_Ht_idx[i*max_parts] * s_track[ r_Ht_idx[i*max_parts+1] ] * s_track[ r_Ht_idx[i*max_parts+2] ] * s_track[ r_Ht_idx[i*max_parts+3] ] * d_const_cd[ r_Ht_idx[i*max_parts+4] ];
      }
    }

    // -- Improved H parallel indexing --
    template<int N, int coefsCount, int size_Ht, int max_terms, int max_parts>
    __device__ __inline__ void
    eval_H_alea6(
        magmaFloatComplex *s_track, magmaFloatComplex &r_cgesvB,
        magmaFloatComplex s_cdt[coefsCount], magma_int_t r_Ht_idx[size_Ht])
    {
      r_cgesvB = MAGMA_C_ZERO;
      #pragma unroll
      for (int i = 0; i < max_terms; i++) {
        r_cgesvB += r_Ht_idx[i*max_parts] * s_track[ r_Ht_idx[i*max_parts+1] ] * s_track[ r_Ht_idx[i*max_parts+2] ] * s_track[ r_Ht_idx[i*max_parts+3] ] * s_cdt[ r_Ht_idx[i*max_parts+4] ];
      }
    }
}

#endif
