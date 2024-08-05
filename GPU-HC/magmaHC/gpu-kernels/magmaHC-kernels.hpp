#ifndef magmaHC_kernels_HPP
#define magmaHC_kernels_HPP
// ============================================================================
// Header file declaring all kernels
//
// Modifications
//    Chiang-Heng Chien  22-10-31:   Initially Created
//    Chiang-Heng Chien  24-06-12:   Add kernels of different settings
//
//> (c) LEMS, Brown University
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// ============================================================================
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <chrono>

#include "magma_v2.h"

//> (#1) PH-(x), RKL-(x), inline(x), LimUnroll(x), L2Cache-(x), TrunPaths-(x) (Naive GPU-HC approach)
real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_P2C(
  magma_queue_t         my_queue,
  int                   sub_RANSAC_iters,
  int                   HC_max_steps, 
  int                   HC_max_correction_steps, 
  int                   HC_delta_t_incremental_steps,
  magmaFloatComplex**   d_startSols_array, 
  magmaFloatComplex**   d_Track_array,
  magma_int_t*          d_Hx_idx_array,
  magma_int_t*          d_Ht_idx_array, 
  magmaFloatComplex_ptr d_phc_coeffs_Hx, 
  magmaFloatComplex_ptr d_phc_coeffs_Ht,
  bool*                 d_is_GPU_HC_Sol_Converge,
  bool*                 d_is_GPU_HC_Sol_Infinity,
  magmaFloatComplex*    d_Debug_Purpose
);

//> (#2) PH-(v), RKL-(x), inline(x), LimUnroll(x), L2Cache-(x), TrunPaths-(x)
real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_PH(
  magma_queue_t         my_queue,
  int                   sub_RANSAC_iters,
  int                   HC_max_steps, 
  int                   HC_max_correction_steps, 
  int                   HC_delta_t_incremental_steps,
  magmaFloatComplex**   d_startSols_array,  
  magmaFloatComplex**   d_Track_array,
  magmaFloatComplex*    d_startParams,      
  magmaFloatComplex*    d_targetParams,
  magmaFloatComplex*    d_diffParams,
  int*                  d_dHdx_indx, 
  int*                  d_dHdt_indx,
  bool*                 d_is_GPU_HC_Sol_Converge,        
  bool*                 d_is_GPU_HC_Sol_Infinity,
  magmaFloatComplex*    d_Debug_Purpose
);

//> (#3) PH-(v), RKL-(v), inline(x), LimUnroll(x), L2Cache-(x), TrunPaths-(x)
real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_PH_RKL(
  magma_queue_t         my_queue,
  int                   sub_RANSAC_iters,
  int                   HC_max_steps, 
  int                   HC_max_correction_steps, 
  int                   HC_delta_t_incremental_steps,
  magmaFloatComplex**   d_startSols_array,  
  magmaFloatComplex**   d_Track_array,
  magmaFloatComplex*    d_startParams,      
  magmaFloatComplex*    d_targetParams,
  magmaFloatComplex*    d_diffParams,
  int*                  d_dHdx_indx, 
  int*                  d_dHdt_indx,
  bool*                 d_is_GPU_HC_Sol_Converge,        
  bool*                 d_is_GPU_HC_Sol_Infinity,
  magmaFloatComplex*    d_Debug_Purpose
);

//> (#4) PH-(v), RKL-(v), inline(v), LimUnroll(x), L2Cache-(x), TrunPaths-(x)
real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_PH_RKL_inline(
  magma_queue_t         my_queue,
  int                   sub_RANSAC_iters,
  int                   HC_max_steps, 
  int                   HC_max_correction_steps, 
  int                   HC_delta_t_incremental_steps,
  magmaFloatComplex**   d_startSols_array,  
  magmaFloatComplex**   d_Track_array,
  magmaFloatComplex*    d_startParams,      
  magmaFloatComplex*    d_targetParams,
  magmaFloatComplex*    d_diffParams,
  int*                  d_dHdx_indx, 
  int*                  d_dHdt_indx,
  bool*                 d_is_GPU_HC_Sol_Converge,        
  bool*                 d_is_GPU_HC_Sol_Infinity,
  magmaFloatComplex*    d_Debug_Purpose
);

//> (#5) PH-(v), RKL-(v), inline(v), LimUnroll(v), L2Cache-(x), TrunPaths-(x)
real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_PH_RKL_inline_LimUnroll(
  magma_queue_t         my_queue,
  int                   sub_RANSAC_iters,
  int                   HC_max_steps, 
  int                   HC_max_correction_steps, 
  int                   HC_delta_t_incremental_steps,
  magmaFloatComplex**   d_startSols_array,  
  magmaFloatComplex**   d_Track_array,
  magmaFloatComplex*    d_startParams,      
  magmaFloatComplex*    d_targetParams,
  magmaFloatComplex*    d_diffParams,
  int*                  d_dHdx_indx, 
  int*                  d_dHdt_indx,
  bool*                 d_is_GPU_HC_Sol_Converge,        
  bool*                 d_is_GPU_HC_Sol_Infinity,
  magmaFloatComplex*    d_Debug_Purpose
);

//> (#6) PH-(v), RKL-(v), inline(x), LimUnroll(v), L2Cache-(v), TrunPaths-(x)
real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_PH_RKL_inline_LimUnroll_L2Cache(
  magma_queue_t         my_queue,
  int                   sub_RANSAC_iters,
  int                   HC_max_steps, 
  int                   HC_max_correction_steps, 
  int                   HC_delta_t_incremental_steps,
  magmaFloatComplex**   d_startSols_array,  
  magmaFloatComplex**   d_Track_array,
  magmaFloatComplex*    d_startParams,      
  magmaFloatComplex*    d_targetParams,
  magmaFloatComplex*    d_diffParams,
  int*                  d_unified_dHdx_dHdt_Index,
  bool*                 d_is_GPU_HC_Sol_Converge,        
  bool*                 d_is_GPU_HC_Sol_Infinity,
  magmaFloatComplex*    d_Debug_Purpose
);

//> (#7) PH-(v), RKL-(v), inline(x), LimUnroll(v), L2Cache-(v), TrunPaths-(v)
real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_PH_RKL_inline_LimUnroll_L2Cache_TrunPaths(
  magma_queue_t       my_queue,
  int                 sub_RANSAC_iters,
  int                 HC_max_steps, 
  int                 HC_max_correction_steps, 
  int                 HC_delta_t_incremental_steps,
  magmaFloatComplex** d_startSols_array, 
  magmaFloatComplex** d_Track_array,
  magmaFloatComplex*  d_startParams,
  magmaFloatComplex*  d_targetParams,
  magmaFloatComplex*  d_diffParams,
  int*                d_unified_dHdx_dHdt_Index,
  bool*               d_is_GPU_HC_Sol_Converge,
  bool*               d_is_GPU_HC_Sol_Infinity,
  magmaFloatComplex*  d_Debug_Purpose
);

//> (#8) PH-(v), RKL-(v), inline(x), LimUnroll(v), L2Cache-(x), TrunPaths-(x)
real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_PH_RKL_LimUnroll(
  magma_queue_t         my_queue,
  int                   sub_RANSAC_iters,
  int                   HC_max_steps, 
  int                   HC_max_correction_steps, 
  int                   HC_delta_t_incremental_steps,
  magmaFloatComplex**   d_startSols_array,  
  magmaFloatComplex**   d_Track_array,
  magmaFloatComplex*    d_startParams,      
  magmaFloatComplex*    d_targetParams,
  magmaFloatComplex*    d_diffParams,
  int*                  d_dHdx_indx, 
  int*                  d_dHdt_indx,
  bool*                 d_is_GPU_HC_Sol_Converge,        
  bool*                 d_is_GPU_HC_Sol_Infinity,
  magmaFloatComplex*    d_Debug_Purpose
);

//> (#9) PH-(v), RKL-(v), inline(x), LimUnroll(v), L2Cache-(v), TrunPaths-(x)
real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_PH_RKL_LimUnroll_L2Cache(
  magma_queue_t       my_queue,
  int                 sub_RANSAC_iters,
  int                 HC_max_steps, 
  int                 HC_max_correction_steps, 
  int                 HC_delta_t_incremental_steps,
  magmaFloatComplex** d_startSols_array, 
  magmaFloatComplex** d_Track_array,
  magmaFloatComplex*  d_startParams,
  magmaFloatComplex*  d_targetParams,
  magmaFloatComplex*  d_diffParams,
  int*                d_unified_dHdx_dHdt_Index,
  bool*               d_is_GPU_HC_Sol_Converge,
  bool*               d_is_GPU_HC_Sol_Infinity,
  magmaFloatComplex*  d_Debug_Purpose
);

//> (#10) PH-(v), RKL-(v), inline(x), LimUnroll(v), L2Cache-(v), TrunPaths-(v)
real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_PH_RKL_LimUnroll_L2Cache_TrunPaths(
  magma_queue_t       my_queue,
  int                 sub_RANSAC_iters,
  int                 HC_max_steps, 
  int                 HC_max_correction_steps, 
  int                 HC_delta_t_incremental_steps,
  magmaFloatComplex** d_startSols_array, 
  magmaFloatComplex** d_Track_array,
  magmaFloatComplex*  d_startParams,
  magmaFloatComplex*  d_targetParams,
  magmaFloatComplex*  d_diffParams,
  int*                d_unified_dHdx_dHdt_Index,
  bool*               d_is_GPU_HC_Sol_Converge,
  bool*               d_is_GPU_HC_Sol_Infinity,
  magmaFloatComplex*  d_Debug_Purpose
);

#endif
