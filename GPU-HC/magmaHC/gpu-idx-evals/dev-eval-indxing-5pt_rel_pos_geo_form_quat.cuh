#ifndef dev_eval_indxing_5pt_rel_pos_geo_form_quat_cuh_
#define dev_eval_indxing_5pt_rel_pos_geo_form_quat_cuh_
// ============================================================================
// Device function for evaluating the parallel indexing for Hx, Ht, and H of
// 5pt_rel_pos_geo_form_quat problem
//
// Modifications
//    Chien  22-11-14:   Originally created
//
// ============================================================================
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cuda_runtime.h>

//> MAGMA
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

//> Macros
#include "../definitions.hpp"

//namespace GPU_Device {

    // -- compute the linear interpolations of parameters of phc
    __device__ __inline__ void
    eval_parameter_homotopy(
        const int tx, float t, 
        magmaFloatComplex *s_phc_coeffs_Hx,
        magmaFloatComplex *s_phc_coeffs_Ht,
        const magmaFloatComplex __restrict__ *d_phc_coeffs_Hx,
        const magmaFloatComplex __restrict__ *d_phc_coeffs_Ht )
    {
        // =============================================================================
        //> parameter homotopy of Hx
        //> floor((75+1)/6) = 12
        #pragma unroll
        for (int i = 0; i < 12; i++) {
          s_phc_coeffs_Hx[ tx + i*NUM_OF_VARS ] = d_phc_coeffs_Hx[ tx*3 + i*NUM_OF_VARS*3 ] 
                                                + d_phc_coeffs_Hx[ tx*3 + 1 + i*NUM_OF_VARS*3 ] * t
                                                + (d_phc_coeffs_Hx[ tx*3 + 2 + i*NUM_OF_VARS*3 ] * t) * t;

          s_phc_coeffs_Ht[ tx + i*NUM_OF_VARS ] = d_phc_coeffs_Ht[ tx*2 + i*NUM_OF_VARS*2 ] 
                                                + d_phc_coeffs_Ht[ tx*2 + 1 + i*NUM_OF_VARS*2 ] * t;
        }

        //> the remaining part
        if (tx < 4) {
          s_phc_coeffs_Hx[ tx + 72 ] = d_phc_coeffs_Hx[ tx*3 + 216 ] 
                                     + d_phc_coeffs_Hx[ tx*3 + 216 + 1 ] * t
                                     + (d_phc_coeffs_Hx[ tx*3 + 216 + 2 ] * t) * t;

          s_phc_coeffs_Ht[ tx + 72 ] = d_phc_coeffs_Ht[ tx*2 + 144 ] 
                                     + d_phc_coeffs_Ht[ tx*2 + 144 + 1 ] * t;
        }
    }

    //> Parallel Jacobian Evaluation ∂H/∂x
    template< int max_terms_parts, int N_max_terms_parts >
    __device__ __inline__ void
    eval_Jacobian_Hx(
        const int tx, magmaFloatComplex *s_track, magmaFloatComplex r_cgesvA[NUM_OF_VARS],
        const int* __restrict__ d_Hx_idx, magmaFloatComplex *s_phc_coeffs )
    {
      #pragma unroll
      for(int i = 0; i < NUM_OF_VARS; i++) {
        r_cgesvA[i] = MAGMA_C_ZERO;

        #pragma unroll
        for(int j = 0; j < HX_MAXIMAL_TERMS; j++) {
          r_cgesvA[i] += d_Hx_idx[j*HX_MAXIMAL_PARTS + i*max_terms_parts + tx*N_max_terms_parts] 
                       * s_phc_coeffs[ d_Hx_idx[j*HX_MAXIMAL_PARTS + 1 + i*max_terms_parts + tx*N_max_terms_parts] ]
                       * s_track[ d_Hx_idx[j*HX_MAXIMAL_PARTS + 2 + i*max_terms_parts + tx*N_max_terms_parts] ]
                       * s_track[ d_Hx_idx[j*HX_MAXIMAL_PARTS + 3 + i*max_terms_parts + tx*N_max_terms_parts] ] ;
        }
      }
    }

    //> Parallel Jacobian Evaluation ∂H/∂t
    template< int max_terms_parts >
    __device__ __inline__ void
    eval_Jacobian_Ht(
        const int tx, magmaFloatComplex *s_track, magmaFloatComplex &r_cgesvB,
        const int* __restrict__ d_Ht_idx, magmaFloatComplex *s_phc_coeffs)
    {
      r_cgesvB = MAGMA_C_ZERO;
      #pragma unroll
      for (int i = 0; i < HT_MAXIMAL_TERMS; i++) {
        r_cgesvB -= d_Ht_idx[i*HT_MAXIMAL_PARTS + tx*max_terms_parts] 
                  * s_phc_coeffs[ d_Ht_idx[i*HT_MAXIMAL_PARTS + 1 + tx*max_terms_parts] ]
                  * s_track[ d_Ht_idx[i*HT_MAXIMAL_PARTS + 2 + tx*max_terms_parts] ] 
                  * s_track[ d_Ht_idx[i*HT_MAXIMAL_PARTS + 3 + tx*max_terms_parts] ] 
                  * s_track[ d_Ht_idx[i*HT_MAXIMAL_PARTS + 4 + tx*max_terms_parts] ];
      }
    }

    //> Parallel Homotopy Evaluation H
    template< int max_terms_parts >
    __device__ __inline__ void
    eval_Homotopy(
        const int tx, magmaFloatComplex *s_track, magmaFloatComplex &r_cgesvB,
        const int* __restrict__ d_Ht_idx, magmaFloatComplex *s_phc_coeffs)
    {
      r_cgesvB = MAGMA_C_ZERO;
      #pragma unroll
      for (int i = 0; i < HT_MAXIMAL_TERMS; i++) {
        r_cgesvB += d_Ht_idx[i*HT_MAXIMAL_PARTS + tx*max_terms_parts] 
                  * s_phc_coeffs[ d_Ht_idx[i*HT_MAXIMAL_PARTS + 1 + tx*max_terms_parts] ]
                  * s_track[ d_Ht_idx[i*HT_MAXIMAL_PARTS + 2 + tx*max_terms_parts] ] 
                  * s_track[ d_Ht_idx[i*HT_MAXIMAL_PARTS + 3 + tx*max_terms_parts] ] 
                  * s_track[ d_Ht_idx[i*HT_MAXIMAL_PARTS + 4 + tx*max_terms_parts] ] ;
      }
    }
//}

#endif
