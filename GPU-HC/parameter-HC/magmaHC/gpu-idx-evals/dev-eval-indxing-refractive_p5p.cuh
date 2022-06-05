#ifndef dev_eval_indxing_refractive_p5p_cuh_
#define dev_eval_indxing_refractive_p5p_cuh_
// ============================================================================
// Device function for evaluating the parallel indexing for Hx, Ht, and H of
// refractive_p5p problem
//
// Modifications
//    Chien  22-01-04:   Originally created
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

    // -- compute the linear interpolations of parameters of phc --
    template<int N, int N_3, int N_2>
    __device__ __inline__ void
    eval_parameter_homotopy(
        const int tx, float t, 
        magmaFloatComplex *s_phc_coeffs_Hx, magmaFloatComplex *s_phc_coeffs_Ht,
        const magmaFloatComplex __restrict__ *d_phc_coeffs_Hx,
        const magmaFloatComplex __restrict__ *d_phc_coeffs_Ht )
    {
        // =============================================================================
        // -- parameter homotopy of Hx --
        // -- floor(90/5) = 18 --
        #pragma unroll
        for (int i = 0; i < 18; i++) {
            s_phc_coeffs_Hx[ tx + i * N ] = d_phc_coeffs_Hx[ tx*3 + i*N_3] 
                                          + d_phc_coeffs_Hx[ tx*3 + 1 + i*N_3 ] * t 
                                          + d_phc_coeffs_Hx[ tx*3 + 2 + i*N_3 ] * t * t;
        }

        // =============================================================================
        // -- parameter homotopy of Ht --
        // -- floor(90/5) = 18 --
        #pragma unroll
        for (int i = 0; i < 18; i++) {
            s_phc_coeffs_Ht[ tx + i * N ] = d_phc_coeffs_Ht[ tx*2 + i*N_2] 
                                          + d_phc_coeffs_Ht[ tx*2 + 1 + i*N_2 ] * t;
        }
    }

    // -- Hx parallel indexing --
    template<int N, int max_terms, int max_parts, int max_terms_parts, int N_max_terms_parts>
    __device__ __inline__ void
    eval_Jacobian_refractive_p5p(
        const int tx, magmaFloatComplex *s_track, magmaFloatComplex r_cgesvA[N],
        const int* __restrict__ d_Hx_idx, magmaFloatComplex *s_phc_coeffs )
    {
      //#pragma unroll
      for(int i = 0; i < N; i++) {
        r_cgesvA[i] = MAGMA_C_ZERO;

        //#pragma unroll
        for(int j = 0; j < max_terms; j++) {
          r_cgesvA[i] += d_Hx_idx[j*max_parts + i*max_terms_parts + tx*N_max_terms_parts] 
                       * s_track[ d_Hx_idx[j*max_parts + 1 + i*max_terms_parts + tx*N_max_terms_parts] ] 
                       * s_track[ d_Hx_idx[j*max_parts + 2 + i*max_terms_parts + tx*N_max_terms_parts] ] 
                       * s_phc_coeffs[ d_Hx_idx[j*max_parts + 3 + i*max_terms_parts + tx*N_max_terms_parts] ];
        }
      }
    }

    // -- Ht parallel indexing --
    template<int N, int max_terms, int max_parts, int max_terms_parts>
    __device__ __inline__ void
    eval_Ht_refractive_p5p(
        const int tx, magmaFloatComplex *s_track, magmaFloatComplex &r_cgesvB,
        const int* __restrict__ d_Ht_idx, magmaFloatComplex *s_phc_coeffs)
    {
      r_cgesvB = MAGMA_C_ZERO;
      //#pragma unroll
      for (int i = 0; i < max_terms; i++) {
        r_cgesvB -= d_Ht_idx[i*max_parts + tx*max_terms_parts] 
                  * s_track[ d_Ht_idx[i*max_parts + 1 + tx*max_terms_parts] ] 
                  * s_track[ d_Ht_idx[i*max_parts + 2 + tx*max_terms_parts] ] 
                  * s_track[ d_Ht_idx[i*max_parts + 3 + tx*max_terms_parts] ] 
                  * s_phc_coeffs[ d_Ht_idx[i*max_parts + 4 + tx*max_terms_parts] ];
      }
    }

    // -- H parallel indexing --
    template<int N, int max_terms, int max_parts, int max_terms_parts>
    __device__ __inline__ void
    eval_H_refractive_p5p(
        const int tx, magmaFloatComplex *s_track, magmaFloatComplex &r_cgesvB,
        const int* __restrict__ d_Ht_idx, magmaFloatComplex *s_phc_coeffs)
    {
      r_cgesvB = MAGMA_C_ZERO;
      //#pragma unroll
      for (int i = 0; i < max_terms; i++) {
        r_cgesvB += d_Ht_idx[i*max_parts + tx*max_terms_parts] 
                  * s_track[ d_Ht_idx[i*max_parts + 1 + tx*max_terms_parts] ] 
                  * s_track[ d_Ht_idx[i*max_parts + 2 + tx*max_terms_parts] ] 
                  * s_track[ d_Ht_idx[i*max_parts + 3 + tx*max_terms_parts] ] 
                  * s_phc_coeffs[ d_Ht_idx[i*max_parts + 4 + tx*max_terms_parts] ];
      }
    }
}

#endif
