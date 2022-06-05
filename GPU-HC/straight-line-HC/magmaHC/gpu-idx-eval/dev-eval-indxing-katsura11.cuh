#ifndef dev_eval_indxing_katsura11_cuh_
#define dev_eval_indxing_katsura11_cuh_
// ============================================================================
// Device function for parallel evaluations of Hx, Ht, and H of
// katsura11 problem
//
// Modifications
//    Chien  22-01-21:   Originally created
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
    eval_cdt_katsura11(
        const int tx, float t, magmaFloatComplex s_vector_cdt[coefsCount],
        magmaFloatComplex r_startCoefs[coefsCount], magmaFloatComplex r_targetCoefs[coefsCount])
    {
        /*#pragma unroll
        for (int i = 0; i < 3; i++) {
            s_vector_cdt[ tx + i * N ] = r_targetCoefs[tx + i * N] * t - r_startCoefs[tx + i * N] * (t-1);
        }*/
        if (tx < 3) {
          s_vector_cdt[ tx ] = r_targetCoefs[tx] * t - r_startCoefs[tx] * (t-1);
        }
    }

    // -- Improved Hx parallel indexing --
    template<int N, int coefsCount, int size_Hx, int max_terms, int max_parts>
    __device__ __inline__ void
    eval_Jacobian_katsura11_improve(
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
                       * s_cdt[ r_Hx_idx[j*max_parts + 2 + i*max_terms*max_parts] ];
        }
      }

      /*#pragma unroll
      for(int i = 0; i < N; i++) {
        r_cgesvA[i] = MAGMA_C_ZERO;

        #pragma unroll
        for(int j = 0; j < 6; j++) {
          r_cgesvA[i] += d_const_Hx_idx[j*4+i*24+tx*192] * s_track[ d_const_Hx_idx[j*4+1+i*24+tx*192] ] * s_track[ d_const_Hx_idx[j*4+2+i*24+tx*192] ] * s_cdt[ d_const_Hx_idx[j*4+3+i*24+tx*192] ];
        }
      }*/
    }

    // -- Improved Ht parallel indexing --
    template<int N, int size_Ht, int max_terms, int max_parts>
    __device__ __inline__ void
    eval_Ht_katsura11_improve(
        magmaFloatComplex *s_track, magmaFloatComplex &r_cgesvB,
        const magmaFloatComplex* __restrict__ d_const_cd, magma_int_t r_Ht_idx[size_Ht])
    {
      r_cgesvB = MAGMA_C_ZERO;
      #pragma unroll
      for (int i = 0; i < max_terms; i++) {
        r_cgesvB += r_Ht_idx[i*max_parts] * s_track[ r_Ht_idx[i*max_parts+1] ] * s_track[ r_Ht_idx[i*max_parts+2] ] * d_const_cd[ r_Ht_idx[i*max_parts+3] ];
      }

      /*r_cgesvB = MAGMA_C_ZERO;
      #pragma unroll
      for (int i = 0; i < 9; i++) {
        r_cgesvB += d_const_Ht_idx[i*5+tx*45] * s_track[ d_const_Ht_idx[i*5+1+tx*45] ] * s_track[ d_const_Ht_idx[i*5+2+tx*45] ] * s_track[ d_const_Ht_idx[i*5+3+tx*45] ] * d_const_vec_cd[ d_const_Ht_idx[i*5+4+tx*45] ];
      }*/
    }

    // -- Improved H parallel indexing --
    template<int N, int coefsCount, int size_Ht, int max_terms, int max_parts>
    __device__ __inline__ void
    eval_H_katsura11_improve(
        magmaFloatComplex *s_track, magmaFloatComplex &r_cgesvB,
        magmaFloatComplex s_cdt[coefsCount], magma_int_t r_Ht_idx[size_Ht])
    {
      r_cgesvB = MAGMA_C_ZERO;
      #pragma unroll
      for (int i = 0; i < max_terms; i++) {
        r_cgesvB += r_Ht_idx[i*max_parts] * s_track[ r_Ht_idx[i*max_parts+1] ] * s_track[ r_Ht_idx[i*max_parts+2] ] * s_cdt[ r_Ht_idx[i*max_parts+3] ];
      }

      /*r_cgesvB = MAGMA_C_ZERO;
      #pragma unroll
      for (int i = 0; i < 9; i++) {
        r_cgesvB += d_const_Ht_idx[i*5+tx*45] * s_track[ d_const_Ht_idx[i*5+1+tx*45] ] * s_track[ d_const_Ht_idx[i*5+2+tx*45] ] * s_track[ d_const_Ht_idx[i*5+3+tx*45] ] * s_cdt[ d_const_Ht_idx[i*5+4+tx*45] ];
      }*/
    }
}

#endif
