#ifndef P2C_DISTORED_2VIEW_TRIANGULATION_H
#define P2C_DISTORED_2VIEW_TRIANGULATION_H
// ==========================================================================================================================
//
// Modifications
//    Chiang-Heng Chien  24-05-29:   Built from generalized three-view relative pose problem.
//
//> (c) LEMS, Brown University
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// ===========================================================================================================================
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <chrono>

#include "magma_v2.h"
#include "magma_lapack.h"
#include "magma_internal.h"
#undef max
#undef min

namespace magmaHCWrapper {

  void p2c_distorted_2view_triangulation(
    magmaFloatComplex *h_targetParams, magmaFloatComplex *h_startParams, 
    magmaFloatComplex *h_phc_coeffs_Hx, magmaFloatComplex *h_phc_coeffs_Ht
  )
  {
    magmaFloatComplex p1 = h_targetParams[0];
    magmaFloatComplex p2 = h_targetParams[1];
    magmaFloatComplex p3 = h_targetParams[2];
    magmaFloatComplex p4 = h_targetParams[3];
    magmaFloatComplex p5 = h_targetParams[4];
    magmaFloatComplex p6 = h_targetParams[5];
    magmaFloatComplex p7 = h_targetParams[6];
    magmaFloatComplex p8 = h_targetParams[7];
    magmaFloatComplex p9 = h_targetParams[8];
    magmaFloatComplex p10 = h_targetParams[9];
    magmaFloatComplex p11 = h_targetParams[10];
    magmaFloatComplex p12 = h_targetParams[11];
    magmaFloatComplex p13 = h_targetParams[12];
    magmaFloatComplex p14 = h_targetParams[13];
    magmaFloatComplex p15 = h_targetParams[14];
    
    magmaFloatComplex q1 = h_startParams[0];
    magmaFloatComplex q2 = h_startParams[1];
    magmaFloatComplex q3 = h_startParams[2];
    magmaFloatComplex q4 = h_startParams[3];
    magmaFloatComplex q5 = h_startParams[4];
    magmaFloatComplex q6 = h_startParams[5];
    magmaFloatComplex q7 = h_startParams[6];
    magmaFloatComplex q8 = h_startParams[7];
    magmaFloatComplex q9 = h_startParams[8];
    magmaFloatComplex q10 = h_startParams[9];
    magmaFloatComplex q11 = h_startParams[10];
    magmaFloatComplex q12 = h_startParams[11];
    magmaFloatComplex q13 = h_startParams[12];
    magmaFloatComplex q14 = h_startParams[13];
    magmaFloatComplex q15 = h_startParams[14];
    
    h_phc_coeffs_Hx[0]=q1;
    h_phc_coeffs_Hx[1]=p1 - q1;
    h_phc_coeffs_Hx[2] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[3] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[4]=q10;
    h_phc_coeffs_Hx[5]=p10 - q10;
    h_phc_coeffs_Hx[6] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[7] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[8]=q10*q14;
    h_phc_coeffs_Hx[9]=p10*q14 + p14*q10 - 2*q10*q14;
    h_phc_coeffs_Hx[10]=p10*p14 - p10*q14 - p14*q10 + q10*q14;
    h_phc_coeffs_Hx[11] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[12]=q10*q15;
    h_phc_coeffs_Hx[13]=p10*q15 + p15*q10 - 2*q10*q15;
    h_phc_coeffs_Hx[14]=p10*p15 - p10*q15 - p15*q10 + q10*q15;
    h_phc_coeffs_Hx[15] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[16]=q11;
    h_phc_coeffs_Hx[17]=p11 - q11;
    h_phc_coeffs_Hx[18] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[19] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[20]=q11*q14;
    h_phc_coeffs_Hx[21]=p11*q14 + p14*q11 - 2*q11*q14;
    h_phc_coeffs_Hx[22]=p11*p14 - p11*q14 - p14*q11 + q11*q14;
    h_phc_coeffs_Hx[23] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[24]=q11*q15;
    h_phc_coeffs_Hx[25]=p11*q15 + p15*q11 - 2*q11*q15;
    h_phc_coeffs_Hx[26]=p11*p15 - p11*q15 - p15*q11 + q11*q15;
    h_phc_coeffs_Hx[27] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[28]=q12;
    h_phc_coeffs_Hx[29]=p12 - q12;
    h_phc_coeffs_Hx[30] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[31] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[32]=q12*q14;
    h_phc_coeffs_Hx[33]=p12*q14 + p14*q12 - 2*q12*q14;
    h_phc_coeffs_Hx[34]=p12*p14 - p12*q14 - p14*q12 + q12*q14;
    h_phc_coeffs_Hx[35] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[36]=q12*q15;
    h_phc_coeffs_Hx[37]=p12*q15 + p15*q12 - 2*q12*q15;
    h_phc_coeffs_Hx[38]=p12*p15 - p12*q15 - p15*q12 + q12*q15;
    h_phc_coeffs_Hx[39] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[40]=q13;
    h_phc_coeffs_Hx[41]=p13 - q13;
    h_phc_coeffs_Hx[42] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[43] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[44]=q13*q14;
    h_phc_coeffs_Hx[45]=p13*q14 + p14*q13 - 2*q13*q14;
    h_phc_coeffs_Hx[46]=p13*p14 - p13*q14 - p14*q13 + q13*q14;
    h_phc_coeffs_Hx[47] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[48]=q13*q14*q15;
    h_phc_coeffs_Hx[49]=p13*q14*q15 + p14*q13*q15 + p15*q13*q14 - 3*q13*q14*q15;
    h_phc_coeffs_Hx[50]=p13*p14*q15 + p13*p15*q14 + p14*p15*q13 - 2*p13*q14*q15 - 2*p14*q13*q15 - 2*p15*q13*q14 + 3*q13*q14*q15;
    h_phc_coeffs_Hx[51]=p13*p14*p15 - p13*p14*q15 - p13*p15*q14 - p14*p15*q13 + p13*q14*q15 + p14*q13*q15 + p15*q13*q14 - q13*q14*q15;
    h_phc_coeffs_Hx[52]=q13*q15;
    h_phc_coeffs_Hx[53]=p13*q15 + p15*q13 - 2*q13*q15;
    h_phc_coeffs_Hx[54]=p13*p15 - p13*q15 - p15*q13 + q13*q15;
    h_phc_coeffs_Hx[55] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[56]=q2;
    h_phc_coeffs_Hx[57]=p2 - q2;
    h_phc_coeffs_Hx[58] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[59] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[60]=q3;
    h_phc_coeffs_Hx[61]=p3 - q3;
    h_phc_coeffs_Hx[62] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[63] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[64]=q4;
    h_phc_coeffs_Hx[65]=p4 - q4;
    h_phc_coeffs_Hx[66] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[67] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[68]=q5;
    h_phc_coeffs_Hx[69]=p5 - q5;
    h_phc_coeffs_Hx[70] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[71] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[72]=q6;
    h_phc_coeffs_Hx[73]=p6 - q6;
    h_phc_coeffs_Hx[74] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[75] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[76]=q7;
    h_phc_coeffs_Hx[77]=p7 - q7;
    h_phc_coeffs_Hx[78] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[79] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[80]=q7*q14;
    h_phc_coeffs_Hx[81]=p7*q14 + p14*q7 - 2*q7*q14;
    h_phc_coeffs_Hx[82]=p7*p14 - p7*q14 - p14*q7 + q7*q14;
    h_phc_coeffs_Hx[83] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[84]=q7*q15;
    h_phc_coeffs_Hx[85]=p7*q15 + p15*q7 - 2*q7*q15;
    h_phc_coeffs_Hx[86]=p7*p15 - p7*q15 - p15*q7 + q7*q15;
    h_phc_coeffs_Hx[87] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[88]=q8;
    h_phc_coeffs_Hx[89]=p8 - q8;
    h_phc_coeffs_Hx[90] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[91] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[92]=q9;
    h_phc_coeffs_Hx[93]=p9 - q9;
    h_phc_coeffs_Hx[94] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[95] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[96] = MAGMA_C_ONE;
    h_phc_coeffs_Hx[97] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[98] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[99] = MAGMA_C_ZERO;

    h_phc_coeffs_Ht[0]=p1 - q1;
    h_phc_coeffs_Ht[1]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[2]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[3]=p10 - q10;
    h_phc_coeffs_Ht[4]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[5]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[6]=p10*q14 + p14*q10 - 2*q10*q14;
    h_phc_coeffs_Ht[7]=2*p10*p14 - 2*p10*q14 - 2*p14*q10 + 2*q10*q14;
    h_phc_coeffs_Ht[8]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[9]=p10*q15 + p15*q10 - 2*q10*q15;
    h_phc_coeffs_Ht[10]=2*p10*p15 - 2*p10*q15 - 2*p15*q10 + 2*q10*q15;
    h_phc_coeffs_Ht[11]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[12]=p11 - q11;
    h_phc_coeffs_Ht[13]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[14]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[15]=p11*q14 + p14*q11 - 2*q11*q14;
    h_phc_coeffs_Ht[16]=2*p11*p14 - 2*p11*q14 - 2*p14*q11 + 2*q11*q14;
    h_phc_coeffs_Ht[17]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[18]=p11*q15 + p15*q11 - 2*q11*q15;
    h_phc_coeffs_Ht[19]=2*p11*p15 - 2*p11*q15 - 2*p15*q11 + 2*q11*q15;
    h_phc_coeffs_Ht[20]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[21]=p12 - q12;
    h_phc_coeffs_Ht[22]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[23]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[24]=p12*q14 + p14*q12 - 2*q12*q14;
    h_phc_coeffs_Ht[25]=2*p12*p14 - 2*p12*q14 - 2*p14*q12 + 2*q12*q14;
    h_phc_coeffs_Ht[26]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[27]=p12*q15 + p15*q12 - 2*q12*q15;
    h_phc_coeffs_Ht[28]=2*p12*p15 - 2*p12*q15 - 2*p15*q12 + 2*q12*q15;
    h_phc_coeffs_Ht[29]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[30]=p13 - q13;
    h_phc_coeffs_Ht[31]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[32]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[33]=p13*q14 + p14*q13 - 2*q13*q14;
    h_phc_coeffs_Ht[34]=2*p13*p14 - 2*p13*q14 - 2*p14*q13 + 2*q13*q14;
    h_phc_coeffs_Ht[35]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[36]=p13*q14*q15 + p14*q13*q15 + p15*q13*q14 - 3*q13*q14*q15;
    h_phc_coeffs_Ht[37]=2*p13*p14*q15 + 2*p13*p15*q14 + 2*p14*p15*q13 - 4*p13*q14*q15 - 4*p14*q13*q15 - 4*p15*q13*q14 + 6*q13*q14*q15;
    h_phc_coeffs_Ht[38]=3*p13*p14*p15 - 3*p13*p14*q15 - 3*p13*p15*q14 - 3*p14*p15*q13 + 3*p13*q14*q15 + 3*p14*q13*q15 + 3*p15*q13*q14 - 3*q13*q14*q15;
    h_phc_coeffs_Ht[39]=p13*q15 + p15*q13 - 2*q13*q15;
    h_phc_coeffs_Ht[40]=2*p13*p15 - 2*p13*q15 - 2*p15*q13 + 2*q13*q15;
    h_phc_coeffs_Ht[41]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[42]=p2 - q2;
    h_phc_coeffs_Ht[43]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[44]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[45]=p3 - q3;
    h_phc_coeffs_Ht[46]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[47]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[48]=p4 - q4;
    h_phc_coeffs_Ht[49]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[50]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[51]=p5 - q5;
    h_phc_coeffs_Ht[52]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[53]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[54]=p6 - q6;
    h_phc_coeffs_Ht[55]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[56]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[57]=p7 - q7;
    h_phc_coeffs_Ht[58]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[59]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[60]=p7*q14 + p14*q7 - 2*q7*q14;
    h_phc_coeffs_Ht[61]=2*p7*p14 - 2*p7*q14 - 2*p14*q7 + 2*q7*q14;
    h_phc_coeffs_Ht[62]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[63]=p7*q15 + p15*q7 - 2*q7*q15;
    h_phc_coeffs_Ht[64]=2*p7*p15 - 2*p7*q15 - 2*p15*q7 + 2*q7*q15;
    h_phc_coeffs_Ht[65]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[66]=p8 - q8;
    h_phc_coeffs_Ht[67]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[68]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[69]=p9 - q9;
    h_phc_coeffs_Ht[70]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[71]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[72]= MAGMA_C_ONE;
    h_phc_coeffs_Ht[73]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[74]= MAGMA_C_ZERO;

  }
} // end of namespace

#endif
