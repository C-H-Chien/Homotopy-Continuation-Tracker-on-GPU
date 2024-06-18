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

//> P2C, Naive GPU-HC approach
real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_P2C(
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

//> 1) 32-bit indices in global memory + no Runge-Kutta in a loop + no inlined device functions + no limited loop unroll + no truncated HC paths
real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_GM32b(
  magma_queue_t         my_queue,
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

real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_GM32b_RKL(
  magma_queue_t         my_queue,
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

//> 2) 32-bit indices in global memory + Runge-Kutta in a loop + inlined device functions + no limited loop unroll + no truncated HC paths
real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_GM32b_RKL_inline(
  magma_queue_t         my_queue,
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

//> 3) 8-bit indices in global memory + Runge-Kutta in a loop + inlined device functions + no limited loop unroll + no truncated HC paths
real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_GM8b_RKL_inline(
  magma_queue_t         my_queue,
  int                   HC_max_steps, 
  int                   HC_max_correction_steps, 
  int                   HC_delta_t_incremental_steps,
  magmaFloatComplex**   d_startSols_array,  
  magmaFloatComplex**   d_Track_array,
  magmaFloatComplex*    d_startParams,      
  magmaFloatComplex*    d_targetParams,
  magmaFloatComplex*    d_diffParams,
  char*                 d_dHdx_indx, 
  char*                 d_dHdt_indx,
  bool*                 d_is_GPU_HC_Sol_Converge,        
  bool*                 d_is_GPU_HC_Sol_Infinity,
  magmaFloatComplex*    d_Debug_Purpose
);

//> 4) 8-bit indices in shared memory + Runge-Kutta in a loop + inlined device functions + no limited loop unroll + no truncated HC paths
real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_SM8b_RKL_inline(
  magma_queue_t         my_queue,
  int                   HC_max_steps, 
  int                   HC_max_correction_steps, 
  int                   HC_delta_t_incremental_steps,
  magmaFloatComplex**   d_startSols_array,  
  magmaFloatComplex**   d_Track_array,
  magmaFloatComplex*    d_startParams,      
  magmaFloatComplex*    d_targetParams,
  magmaFloatComplex*    d_diffParams,
  char*                 d_dHdx_indx, 
  char*                 d_dHdt_indx,
  bool*                 d_is_GPU_HC_Sol_Converge,        
  bool*                 d_is_GPU_HC_Sol_Infinity,
  magmaFloatComplex*    d_Debug_Purpose
);

//> 5) 8-bit indices in shared memory + Runge-Kutta in a loop + no inlined device functions + no limited loop unroll + no truncated HC paths
real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_SM8b_RKL(
  magma_queue_t         my_queue,
  int                   HC_max_steps, 
  int                   HC_max_correction_steps, 
  int                   HC_delta_t_incremental_steps,
  magmaFloatComplex**   d_startSols_array,  
  magmaFloatComplex**   d_Track_array,
  magmaFloatComplex*    d_startParams,      
  magmaFloatComplex*    d_targetParams,
  magmaFloatComplex*    d_diffParams,
  char*                 d_dHdx_indx, 
  char*                 d_dHdt_indx,
  bool*                 d_is_GPU_HC_Sol_Converge,        
  bool*                 d_is_GPU_HC_Sol_Infinity,
  magmaFloatComplex*    d_Debug_Purpose
);

//> 6) 8-bit indices in shared memory + Runge-Kutta in a loop + no inlined device functions + limited loop unroll + no truncated HC paths
real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_SM8b_RKL_LimUnroll(
  magma_queue_t         my_queue,
  int                   HC_max_steps, 
  int                   HC_max_correction_steps, 
  int                   HC_delta_t_incremental_steps,
  magmaFloatComplex**   d_startSols_array,  
  magmaFloatComplex**   d_Track_array,
  magmaFloatComplex*    d_startParams,      
  magmaFloatComplex*    d_targetParams,
  magmaFloatComplex*    d_diffParams,
  char*                  d_dHdx_indx, 
  char*                  d_dHdt_indx,
  bool*                 d_is_GPU_HC_Sol_Converge,        
  bool*                 d_is_GPU_HC_Sol_Infinity,
  magmaFloatComplex*    d_Debug_Purpose
);

//> 7) 32-bit indices in global memory + Runge-Kutta in a loop + no inlined device functions + limited loop unroll + no truncated HC paths
real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_GM32b_RKL_LimUnroll(
  magma_queue_t         my_queue,
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

//> 8) 32-bit indices in global memory + Runge-Kutta in a loop + inlined device functions + limited loop unroll + no truncated HC paths
real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_GM32b_RKL_inline_LimUnroll(
  magma_queue_t         my_queue,
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

//> 9) 32-bit indices in global memory + Runge-Kutta in a loop + inlined device functions + limited loop unroll + truncated HC paths
real_Double_t kernel_GPUHC_trifocal_2op1p_30x30_GM32b_RKL_inline_LimUnroll_TrunPaths(
  magma_queue_t         my_queue,
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

#endif
