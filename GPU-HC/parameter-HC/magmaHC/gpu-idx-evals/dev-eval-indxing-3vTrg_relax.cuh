#ifndef dev_eval_indxing_3vTrg_relax_cuh_
#define dev_eval_indxing_3vTrg_relax_cuh_
// ============================================================================
// Device function for evaluating the parallel indexing for Hx, Ht, and H of
// 3vTrg_relax problem
//
// Modifications
//    Chien  22-01-02:   Originally created
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
    template<int N, int N_2>
    __device__ __inline__ void
    eval_parameter_homotopy(
        const int tx, float t, 
        magmaFloatComplex *s_phc_coeffs_Hx, magmaFloatComplex *s_phc_coeffs_Ht,
        const magmaFloatComplex __restrict__ *d_phc_coeffs_Hx,
        const magmaFloatComplex __restrict__ *d_phc_coeffs_Ht )
    {
        // =============================================================================
        // -- parameter homotopy of Hx --
        // -- floor(25/8) = 3 --
        #pragma unroll
        for (int i = 0; i < 3; i++) {
            s_phc_coeffs_Hx[ tx + i * N ] = d_phc_coeffs_Hx[ tx*2 + i*N_2] 
                                          + d_phc_coeffs_Hx[ tx*2 + 1 + i*N_2 ] * t;
        }

        // -- the remaining parts --
        // -- 25 - 8*3 = 1 --
        // -- N_2 * 3 = 48 --
        // -- last s_phc_coeffs_Hx from the above for-loop: 7 + 2*8=23
        if (tx < 1) {
          s_phc_coeffs_Hx[ tx + 24 ] = d_phc_coeffs_Hx[ tx*2 + 48] 
                                     + d_phc_coeffs_Hx[ tx*2 + 1 + 48] * t;
        }

        // =============================================================================
        // -- parameter homotopy of Ht --
        // -- floor(25/8) = 3 --
        #pragma unroll
        for (int i = 0; i < 3; i++) {
            s_phc_coeffs_Ht[ tx + i * N ] = d_phc_coeffs_Ht[ tx + i * N ];
        }

        // -- the remaining parts --
        if (tx < 7) {
          s_phc_coeffs_Ht[ tx + 24 ] = d_phc_coeffs_Ht[ tx + 24 ];
        }

    }

    // -- Hx parallel indexing --
    template<int N, int max_terms, int max_parts, int max_terms_parts, int N_max_terms_parts>
    __device__ __inline__ void
    eval_Jacobian_3vTrg_relax(
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
                       * s_phc_coeffs[ d_Hx_idx[j*max_parts + 2 + i*max_terms_parts + tx*N_max_terms_parts] ];
        }
      }
    }

    // -- Ht parallel indexing --
    template<int N, int max_terms, int max_parts, int max_terms_parts>
    __device__ __inline__ void
    eval_Ht_3vTrg_relax(
        const int tx, magmaFloatComplex *s_track, magmaFloatComplex &r_cgesvB,
        const int* __restrict__ d_Ht_idx, magmaFloatComplex *s_phc_coeffs)
    {
      r_cgesvB = MAGMA_C_ZERO;
      //#pragma unroll
      for (int i = 0; i < max_terms; i++) {
        r_cgesvB -= d_Ht_idx[i*max_parts + tx*max_terms_parts] 
                  * s_track[ d_Ht_idx[i*max_parts + 1 + tx*max_terms_parts] ] 
                  * s_track[ d_Ht_idx[i*max_parts + 2 + tx*max_terms_parts] ] 
                  * s_phc_coeffs[ d_Ht_idx[i*max_parts + 3 + tx*max_terms_parts] ];
      }
    }

    // -- H parallel indexing --
    template<int N, int max_terms, int max_parts, int max_terms_parts>
    __device__ __inline__ void
    eval_H_3vTrg_relax(
        const int tx, magmaFloatComplex *s_track, magmaFloatComplex &r_cgesvB,
        const int* __restrict__ d_Ht_idx, magmaFloatComplex *s_phc_coeffs)
    {
      r_cgesvB = MAGMA_C_ZERO;
      //#pragma unroll
      for (int i = 0; i < max_terms; i++) {
        r_cgesvB += d_Ht_idx[i*max_parts + tx*max_terms_parts] 
                  * s_track[ d_Ht_idx[i*max_parts + 1 + tx*max_terms_parts] ] 
                  * s_track[ d_Ht_idx[i*max_parts + 2 + tx*max_terms_parts] ] 
                  * s_phc_coeffs[ d_Ht_idx[i*max_parts + 3 + tx*max_terms_parts] ];
      }
    }
}

#endif
