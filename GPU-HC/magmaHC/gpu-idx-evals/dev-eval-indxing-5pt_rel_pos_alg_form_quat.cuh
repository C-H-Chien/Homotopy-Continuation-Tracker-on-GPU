#ifndef dev_eval_indxing_5pt_rel_pos_alg_form_quat_cuh_
#define dev_eval_indxing_5pt_rel_pos_alg_form_quat_cuh_
// =================================================================================================
// Device function for evaluating the Jacobians ∂H/∂x, ∂H/∂t, and H
//
// Modifications
//    Chiang-Heng Chien  22-11-16:   Originally created as an example problem for GPU-HC solver
//    Chiang-Heng Chien  23-12-28:   Add macros for easy variable reference
//
// =================================================================================================
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

    //> Compute the linear interpolations of parameters of phc
    template< unsigned Full_Parallel_Offset, \
              unsigned Partial_Parallel_Thread_Offset, \
              unsigned Partial_Parallel_Index_Offset, \
              unsigned Max_Order_of_t_Plus_One, \
              unsigned Partial_Parallel_Index_Offset_Hx, \
              unsigned Partial_Parallel_Index_Offset_Ht >
    __device__ __inline__ void
    eval_parameter_homotopy(
        const int tx, float t, 
        magmaFloatComplex *s_phc_coeffs_Hx,
        magmaFloatComplex *s_phc_coeffs_Ht,
        const magmaFloatComplex __restrict__ *d_phc_coeffs_Hx,
        const magmaFloatComplex __restrict__ *d_phc_coeffs_Ht )
    {
        // =============================================================================
        //> parameter homotopy for evaluating ∂H/∂x and ∂H/∂t
        #pragma unroll
        for (int i = 0; i < Full_Parallel_Offset; i++) {
          s_phc_coeffs_Hx[ tx + i*NUM_OF_VARS ] = d_phc_coeffs_Hx[ tx*Max_Order_of_t_Plus_One + i*NUM_OF_VARS*Max_Order_of_t_Plus_One ] 
                                                + d_phc_coeffs_Hx[ tx*Max_Order_of_t_Plus_One + 1 + i*NUM_OF_VARS*Max_Order_of_t_Plus_One ] * t
                                                + d_phc_coeffs_Hx[ tx*Max_Order_of_t_Plus_One + 2 + i*NUM_OF_VARS*Max_Order_of_t_Plus_One ] * t * t;

          s_phc_coeffs_Ht[ tx + i*NUM_OF_VARS ] = d_phc_coeffs_Ht[ tx*MAX_ORDER_OF_T + i*NUM_OF_VARS*MAX_ORDER_OF_T ] 
                                                + d_phc_coeffs_Ht[ tx*MAX_ORDER_OF_T + 1 + i*NUM_OF_VARS*MAX_ORDER_OF_T ] * t;
        }

        //> the remaining part
        if (tx < Partial_Parallel_Thread_Offset) {
          s_phc_coeffs_Hx[ tx + Partial_Parallel_Index_Offset ] = \
              d_phc_coeffs_Hx[ tx*Max_Order_of_t_Plus_One + Partial_Parallel_Index_Offset_Hx ] 
            + d_phc_coeffs_Hx[ tx*Max_Order_of_t_Plus_One + Partial_Parallel_Index_Offset_Hx + 1 ] * t
            + (d_phc_coeffs_Hx[ tx*Max_Order_of_t_Plus_One + Partial_Parallel_Index_Offset_Hx + 2 ] * t) * t;

          s_phc_coeffs_Ht[ tx + Partial_Parallel_Index_Offset ] = \
              d_phc_coeffs_Ht[ tx*MAX_ORDER_OF_T + Partial_Parallel_Index_Offset_Ht ] 
            + d_phc_coeffs_Ht[ tx*MAX_ORDER_OF_T + Partial_Parallel_Index_Offset_Ht + 1 ] * t;
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

    //> Parallel Homotopy Evaluation
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
