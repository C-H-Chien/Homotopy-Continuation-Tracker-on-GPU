#ifndef dev_eval_indxing_trifocal_2op1p_30_LIMUNROLL_L2CACHE_CUH_
#define dev_eval_indxing_trifocal_2op1p_30_LIMUNROLL_L2CACHE_CUH_
// ===============================================================================================
// Code Description: Device function for evaluating the parallel indexing of the Jacobians w.r.t.
//                   the unknowns x (Hx), the variable t (Ht), and the homotopy itself
// the trifocal 2op1p 30x30 problem
//
// Major Modifications
//    Chiang-Heng Chien  22-10-03:   Edited from the first version
//                                   (dev-eval-indexing-trifocal_2op1p_30.cuh)
//
// ===============================================================================================
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>

#include <cuda_runtime.h>

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

#define __HC_INLINE__ __noinline__

//> compute the parameter homotopy
template < typename T, int Num_Of_Vars >
__device__ __HC_INLINE__ void
compute_param_homotopy(
  const int tx, T t,
  magmaFloatComplex *s_param_homotopy,
  magmaFloatComplex *s_start_params,
  magmaFloatComplex *s_target_params
)
{
    s_param_homotopy[ tx ] = s_target_params[ tx ] * t + s_start_params[ tx ] * (1.0-t);

    if (tx < 3) {
      s_param_homotopy[ tx + Num_Of_Vars ] = s_target_params[ tx + Num_Of_Vars ] * t + s_start_params[ tx + Num_Of_Vars ] * (1.0-t);
    }
}

//> Jacobian \partial H / \partial x parallel evaluation
template< int Num_Of_Vars, int dHdx_Max_Terms, int dHdx_Max_Parts, int dHdx_Entry_Offset >
__device__ __HC_INLINE__ void
eval_Jacobian_Hx(
    const int tx,                                     //> thread id
    magmaFloatComplex r_cgesvA[Num_Of_Vars],          //> each row of the Jacobian matrix
    magmaFloatComplex *s_vars,                        //> variables
    magmaFloatComplex *s_start_params,                //> start parameters
    magmaFloatComplex *s_target_params,               //> target parameters
    magmaFloatComplex *s_param_homotopy,              //> parameter homotopy
    const int* __restrict__ unified_dHdx_dHdt_indices //> evaluation indices
)
{
  //> Full, explicit form of evaluation
  for(int i = 0; i < Num_Of_Vars; i++) {

    //> initialize to zero
    r_cgesvA[i] = MAGMA_C_ZERO;

    //> With transpose...
    #pragma unroll 2
    // #pragma unroll
    for(int j = 0; j < dHdx_Max_Terms; j++) {

      //> compute the element of the Jacobian matrix Hx
      r_cgesvA[i] += unified_dHdx_dHdt_indices[(i*dHdx_Entry_Offset)*Num_Of_Vars + j*dHdx_Max_Parts*Num_Of_Vars + tx]
                    * s_param_homotopy[ unified_dHdx_dHdt_indices[ (i*dHdx_Entry_Offset)*Num_Of_Vars + (j*dHdx_Max_Parts + 1)*Num_Of_Vars + tx ] ]
                    * s_param_homotopy[ unified_dHdx_dHdt_indices[ (i*dHdx_Entry_Offset)*Num_Of_Vars + (j*dHdx_Max_Parts + 2)*Num_Of_Vars + tx ] ]
                    * s_vars[           unified_dHdx_dHdt_indices[ (i*dHdx_Entry_Offset)*Num_Of_Vars + (j*dHdx_Max_Parts + 3)*Num_Of_Vars + tx ] ]
                    * s_vars[           unified_dHdx_dHdt_indices[ (i*dHdx_Entry_Offset)*Num_Of_Vars + (j*dHdx_Max_Parts + 4)*Num_Of_Vars + tx ] ];
    }
  }
}

//> Jacobian \partial H / \partial t parallel evaluation
template< int Num_Of_Vars, int dHdt_Max_Terms, int dHdt_Max_Parts, int dHdx_Index_Size >
__device__ __HC_INLINE__ void
eval_Jacobian_Ht(
    const int tx,                                       //> thread id
    magmaFloatComplex &r_cgesvB,                        //> each row of the Jacobian matrix
    magmaFloatComplex *s_vars,                          //> variables
    magmaFloatComplex *s_start_params,                  //> start parameters
    magmaFloatComplex *s_target_params,                 //> target parameters
    magmaFloatComplex *s_param_homotopy,                //> parameter homotopy
    const int* __restrict__ unified_dHdx_dHdt_indices,  //> evaluation indices
    magmaFloatComplex *s_diffParams
)
{
  //> initialize each element to 0
  r_cgesvB = MAGMA_C_ZERO;

  #pragma unroll 2
  // #pragma unroll
  for (int i = 0; i < dHdt_Max_Terms; i++) {

    //> With transpose...
    r_cgesvB -= unified_dHdx_dHdt_indices[dHdx_Index_Size + i*dHdt_Max_Parts*Num_Of_Vars + tx]
              * (s_diffParams[unified_dHdx_dHdt_indices[ dHdx_Index_Size + (i*dHdt_Max_Parts + 1)*Num_Of_Vars + tx ]] * s_param_homotopy[ unified_dHdx_dHdt_indices[ dHdx_Index_Size + (i*dHdt_Max_Parts + 2)*Num_Of_Vars + tx ] ]
               + s_diffParams[unified_dHdx_dHdt_indices[ dHdx_Index_Size + (i*dHdt_Max_Parts + 2)*Num_Of_Vars + tx ]] * s_param_homotopy[ unified_dHdx_dHdt_indices[ dHdx_Index_Size + (i*dHdt_Max_Parts + 1)*Num_Of_Vars + tx ] ] )
              * s_vars[       unified_dHdx_dHdt_indices[ dHdx_Index_Size + (i*dHdt_Max_Parts + 3)*Num_Of_Vars + tx ] ]
              * s_vars[       unified_dHdx_dHdt_indices[ dHdx_Index_Size + (i*dHdt_Max_Parts + 4)*Num_Of_Vars + tx ] ]
              * s_vars[       unified_dHdx_dHdt_indices[ dHdx_Index_Size + (i*dHdt_Max_Parts + 5)*Num_Of_Vars + tx ] ];
  }
}

//> Homotopy evaluation
template< int Num_Of_Vars, int dHdt_Max_Terms, int dHdt_Max_Parts, int dHdx_Index_Size >
__device__ __HC_INLINE__ void
eval_Homotopy(
    const int tx,                                       //> thread id
    magmaFloatComplex &r_cgesvB,                        //> each row of the parameter homotopy
    magmaFloatComplex *s_vars,                          //> variables
    magmaFloatComplex *s_start_params,                  //> start parameters
    magmaFloatComplex *s_target_params,                 //> target parameters
    magmaFloatComplex *s_param_homotopy,                //> parameter homotopy
    const int* __restrict__ unified_dHdx_dHdt_indices   //> evaluation indices
)
{
  //> initialize each element to 0
  r_cgesvB = MAGMA_C_ZERO;

  #pragma unroll 2
  // #pragma unroll
  for (int i = 0; i < dHdt_Max_Terms; i++) {

    r_cgesvB += unified_dHdx_dHdt_indices[ dHdx_Index_Size + i*dHdt_Max_Parts*Num_Of_Vars + tx ]
              * s_param_homotopy[ unified_dHdx_dHdt_indices[ dHdx_Index_Size + (i*dHdt_Max_Parts + 1)*Num_Of_Vars + tx ] ]
              * s_param_homotopy[ unified_dHdx_dHdt_indices[ dHdx_Index_Size + (i*dHdt_Max_Parts + 2)*Num_Of_Vars + tx ] ]
              * s_vars[           unified_dHdx_dHdt_indices[ dHdx_Index_Size + (i*dHdt_Max_Parts + 3)*Num_Of_Vars + tx ] ]
              * s_vars[           unified_dHdx_dHdt_indices[ dHdx_Index_Size + (i*dHdt_Max_Parts + 4)*Num_Of_Vars + tx ] ]
              * s_vars[           unified_dHdx_dHdt_indices[ dHdx_Index_Size + (i*dHdt_Max_Parts + 5)*Num_Of_Vars + tx ] ];
  }
}

#endif
