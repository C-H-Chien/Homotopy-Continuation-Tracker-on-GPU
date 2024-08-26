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

//> (#1) PH-(x), CodeOpt-(x), TrunPaths-(x), TrunRANSAC-(x) (Naive GPU-HC approach)
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

//> (#2) PH-(v), CodeOpt-(x), TrunPaths-(x), TrunRANSAC-(x)
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

//> (#3-Ampere and above) PH-(v), CodeOpt-(v), TrunPaths-(x), TrunRANSAC-(x)
real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_PH_CodeOpt(
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

//> (#3-Volta) PH-(v), CodeOpt-(v), TrunPaths-(x), TrunRANSAC-(x)
real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_PH_CodeOpt_Volta(
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
  int*                d_dHdx_indx, 
  int*                d_dHdt_indx,
  bool*               d_is_GPU_HC_Sol_Converge,
  bool*               d_is_GPU_HC_Sol_Infinity,
  magmaFloatComplex*  d_Debug_Purpose
);

//> (#4-Ampere and above) PH-(v), CodeOpt-(v), TrunPaths-(v), TrunRANSAC-(x)
real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_PH_CodeOpt_TrunPaths(
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

//> (#4-Volta) PH-(v), CodeOpt-(v), TrunPaths-(v), TrunRANSAC-(x)
real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_PH_CodeOpt_TrunPaths_Volta(
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
  int*                d_dHdx_indx, 
  int*                d_dHdt_indx,
  bool*               d_is_GPU_HC_Sol_Converge,
  bool*               d_is_GPU_HC_Sol_Infinity,
  magmaFloatComplex*  d_Debug_Purpose
);

//> (#5-Ampere and above) PH-(v), CodeOpt-(v), TrunPaths-(v), TrunRANSAC-(v)
real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_PH_CodeOpt_TrunPaths_TrunRANSAC(
  magma_queue_t       my_queue,
  int                 sub_RANSAC_iters,
  int                 Num_Of_Triplet_Edgels,
  int                 HC_max_steps, 
  int                 HC_max_correction_steps, 
  int                 HC_delta_t_incremental_steps,
  magmaFloatComplex** d_startSols_array, 
  magmaFloatComplex** d_Track_array,
  magmaFloatComplex*  d_startParams,
  magmaFloatComplex*  d_targetParams,
  magmaFloatComplex*  d_diffParams,
  int*                d_unified_dHdx_dHdt_Index,
  float*              d_Triplet_Edge_Locations,
  float*              d_Intrinsic_Matrix,
  bool*               d_is_GPU_HC_Sol_Converge,
  bool*               d_is_GPU_HC_Sol_Infinity,
  magmaFloatComplex*  d_Debug_Purpose,
  bool*               d_Found_Trifocal_Sols, 
  int*                d_Trifocal_Sols_Batch_Index
);

//> (#5-Volta) PH-(v), CodeOpt-(v), TrunPaths-(v), TrunRANSAC-(v)
real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_PH_CodeOpt_TrunPaths_TrunRANSAC_Volta(
  magma_queue_t       my_queue,
  int                 sub_RANSAC_iters,
  int                 Num_Of_Triplet_Edgels,
  int                 HC_max_steps, 
  int                 HC_max_correction_steps, 
  int                 HC_delta_t_incremental_steps,
  magmaFloatComplex** d_startSols_array, 
  magmaFloatComplex** d_Track_array,
  magmaFloatComplex*  d_startParams,
  magmaFloatComplex*  d_targetParams,
  magmaFloatComplex*  d_diffParams,
  int*                d_dHdx_indx, 
  int*                d_dHdt_indx,
  float*              d_Triplet_Edge_Locations,
  float*              d_Intrinsic_Matrix,
  bool*               d_is_GPU_HC_Sol_Converge,
  bool*               d_is_GPU_HC_Sol_Infinity,
  magmaFloatComplex*  d_Debug_Purpose,
  bool*               d_Found_Trifocal_Sols, 
  int*                d_Trifocal_Sols_Batch_Index
);

#endif
