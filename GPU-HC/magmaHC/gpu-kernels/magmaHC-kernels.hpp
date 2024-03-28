#ifndef magmaHC_kernels_HPP
#define magmaHC_kernels_HPP
// ============================================================================
// Header file declaring all kernels
//
// Modifications
//    Chiang-Heng Chien  22-10-31:   Initially Created (Copied from other repos)
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


//> Geometric Form
real_Double_t kernel_HC_Solver_5pt_rel_pos_geo_form_quat(                      
  magma_queue_t         my_queue, 
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

//> Algebraic Form
real_Double_t kernel_HC_Solver_5pt_rel_pos_alg_form_quat(                      
  magma_queue_t         my_queue, 
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

//> Trifocal Relative Pose from Lines at Points (30x30 Form)
real_Double_t kernel_HC_Solver_trifocal_2op1p_30x30(
  magma_queue_t         my_queue,
  int                   HC_max_steps, 
  int                   HC_max_correction_steps, 
  int                   HC_delta_t_incremental_steps,
  magmaFloatComplex**   d_startSols_array,  
  magmaFloatComplex**   d_Track_array,
  magmaFloatComplex*    d_startParams,      
  magmaFloatComplex*    d_targetParams,
  magmaFloatComplex*    d_diffParams,
  int*                  d_Hx_indx, 
  int*                  d_Ht_indx,
  bool*                 d_is_GPU_HC_Sol_Converge,        
  bool*                 d_is_GPU_HC_Sol_Infinity,
  magmaFloatComplex*    d_Debug_Purpose
);

real_Double_t kernel_HC_Solver_trifocal_2op1p_30x30(
  magma_queue_t         my_queue,
  int                   HC_max_steps, 
  int                   HC_max_correction_steps, 
  int                   HC_delta_t_incremental_steps,
  magmaFloatComplex**   d_startSols_array,  
  magmaFloatComplex**   d_Track_array,
  magmaFloatComplex*    d_startParams,      
  magmaFloatComplex*    d_targetParams,
  magmaFloatComplex*    d_diffParams,
  char*        d_Hx_indx, 
  char*        d_Ht_indx,
  bool*                 d_is_GPU_HC_Sol_Converge,        
  bool*                 d_is_GPU_HC_Sol_Infinity,
  magmaFloatComplex*    d_Debug_Purpose
);

//> Trifocal Relative Pose from Lines at Points (30x30 Form) P2C Version
real_Double_t kernel_HC_Solver_trifocal_2op1p_30x30_P2C (
  magma_queue_t         my_queue,
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

#endif
