#ifndef CPU_EVAL_INDX_TRICOAL_2OP1P_30X30_HPP
#define CPU_EVAL_INDX_TRICOAL_2OP1P_30X30_HPP
// ============================================================================
// Evaluations of Jacobians dH/dX and dH/dt for the trifocal 2op1p 30x30 problem
// This is provided by Macauley2 first and reformulated from MATLAB scripts
// reformateEvalFromM2.m and reformat_y.m
//
// Modifications
//    Chien  24-09-19:   Initially created.
//
// ============================================================================
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>

#include "magma_v2.h"
#include "magma_lapack.h"
#include "magma_internal.h"
#undef max
#undef min

inline void cpu_eval_indx_dHdX_trifocal_2op1p_30(
  const int Num_Of_Vars,       const int dHdx_Max_Terms, 
  const int dHdx_Entry_Offset, const int dHdx_Max_Parts, 
  const int* dHdx_Index,
  magmaFloatComplex* variable, magmaFloatComplex* param_homotopy,
  magmaFloatComplex* cgesvA )
{
  //> dH/dX
  for(int i = 0; i < Num_Of_Vars*Num_Of_Vars; i++) {

    int row_offset = floor(i / Num_Of_Vars);
    int column_offset = i % Num_Of_Vars;
    int index = row_offset + column_offset*Num_Of_Vars;

    cgesvA[index] = MAGMA_C_ZERO;
    for(int j = 0; j < dHdx_Max_Terms; j++) {

      //> compute the element of the Jacobian matrix Hx
      cgesvA[index] += dHdx_Index[(column_offset*dHdx_Entry_Offset)*Num_Of_Vars + j*dHdx_Max_Parts*Num_Of_Vars + row_offset]
                 * param_homotopy[ dHdx_Index[ (column_offset*dHdx_Entry_Offset)*Num_Of_Vars + (j*dHdx_Max_Parts + 1)*Num_Of_Vars + row_offset ] ]
                 * param_homotopy[ dHdx_Index[ (column_offset*dHdx_Entry_Offset)*Num_Of_Vars + (j*dHdx_Max_Parts + 2)*Num_Of_Vars + row_offset ] ]
                 * variable[       dHdx_Index[ (column_offset*dHdx_Entry_Offset)*Num_Of_Vars + (j*dHdx_Max_Parts + 3)*Num_Of_Vars + row_offset ] ]
                 * variable[       dHdx_Index[ (column_offset*dHdx_Entry_Offset)*Num_Of_Vars + (j*dHdx_Max_Parts + 4)*Num_Of_Vars + row_offset ] ];
    }
  }
}

inline void cpu_eval_indx_dHdt_trifocal_2op1p_30(
  const int Num_Of_Vars,       const int dHdt_Max_Terms, 
  const int dHdt_Max_Parts,    const int* h_dHdt_Index,
  magmaFloatComplex* variable, magmaFloatComplex* param_homotopy,
  magmaFloatComplex* cgesvB,   magmaFloatComplex* diff_params )
{
  //> dH/dt
  for (int i = 0; i < Num_Of_Vars; i++) {

    cgesvB[i] = MAGMA_C_ZERO;
    for (int j = 0; j < dHdt_Max_Terms; j++) {

      cgesvB[i] -= h_dHdt_Index[j*dHdt_Max_Parts*Num_Of_Vars + i]
                    * (diff_params[h_dHdt_Index[ (j*dHdt_Max_Parts + 1)*Num_Of_Vars + i ]] * param_homotopy[ h_dHdt_Index[ (j*dHdt_Max_Parts + 2)*Num_Of_Vars + i ] ]
                     + diff_params[h_dHdt_Index[ (j*dHdt_Max_Parts + 2)*Num_Of_Vars + i ]] * param_homotopy[ h_dHdt_Index[ (j*dHdt_Max_Parts + 1)*Num_Of_Vars + i ] ] )
                    * variable[    h_dHdt_Index[ (j*dHdt_Max_Parts + 3)*Num_Of_Vars + i ] ]
                    * variable[    h_dHdt_Index[ (j*dHdt_Max_Parts + 4)*Num_Of_Vars + i ] ]
                    * variable[    h_dHdt_Index[ (j*dHdt_Max_Parts + 5)*Num_Of_Vars + i ] ];
    }
  }
}

inline void cpu_eval_indx_H_trifocal_2op1p_30(
  const int Num_Of_Vars,       const int dHdt_Max_Terms, 
  const int dHdt_Max_Parts,    const int* h_dHdt_Index,
  magmaFloatComplex* variable, magmaFloatComplex* param_homotopy, 
  magmaFloatComplex* cgesvB )
{
  for (int i = 0; i < Num_Of_Vars; i++) {

    cgesvB[i] = MAGMA_C_ZERO;
    for (int j = 0; j < dHdt_Max_Terms; j++) {
      cgesvB[i] += h_dHdt_Index[j*dHdt_Max_Parts*Num_Of_Vars + (i) ]
              * param_homotopy[ h_dHdt_Index[ (j*dHdt_Max_Parts + 1)*Num_Of_Vars + (i) ] ]
              * param_homotopy[ h_dHdt_Index[ (j*dHdt_Max_Parts + 2)*Num_Of_Vars + (i) ] ]
              * variable[       h_dHdt_Index[ (j*dHdt_Max_Parts + 3)*Num_Of_Vars + (i) ] ]
              * variable[       h_dHdt_Index[ (j*dHdt_Max_Parts + 4)*Num_Of_Vars + (i) ] ]
              * variable[       h_dHdt_Index[ (j*dHdt_Max_Parts + 5)*Num_Of_Vars + (i) ] ];
    }
  }
}

#endif
