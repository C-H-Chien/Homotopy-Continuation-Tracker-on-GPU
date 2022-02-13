#ifndef p2c_3vTrg_relax_h
#define p2c_3vTrg_relax_h
// =======================================================================
//
// Modifications
//    Chien  22-01-02:   Initially Created
//
// =======================================================================
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <chrono>

// -- magma --
#include "flops.h"
#include "magma_v2.h"
#include "magma_lapack.h"
#include "magma_internal.h"
#undef max
#undef min

namespace magmaHCWrapper {

  void p2c_3vTrg_relax(
    magmaFloatComplex *h_targetParams, magmaFloatComplex *h_startParams, 
    magmaFloatComplex *h_phc_coeffs_H, magmaFloatComplex *h_phc_coeffs_Ht
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
    magmaFloatComplex p16 = h_targetParams[15];
    magmaFloatComplex p17 = h_targetParams[16];
    magmaFloatComplex p18 = h_targetParams[17];
    magmaFloatComplex p19 = h_targetParams[18];
    magmaFloatComplex p20 = h_targetParams[19];
    magmaFloatComplex p21 = h_targetParams[20];
    magmaFloatComplex p22 = h_targetParams[21];
    magmaFloatComplex p23 = h_targetParams[22];
    magmaFloatComplex p24 = h_targetParams[23];

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
    magmaFloatComplex q16 = h_startParams[15];
    magmaFloatComplex q17 = h_startParams[16];
    magmaFloatComplex q18 = h_startParams[17];
    magmaFloatComplex q19 = h_startParams[18];
    magmaFloatComplex q20 = h_startParams[19];
    magmaFloatComplex q21 = h_startParams[20];
    magmaFloatComplex q22 = h_startParams[21];
    magmaFloatComplex q23 = h_startParams[22];
    magmaFloatComplex q24 = h_startParams[23];

    h_phc_coeffs_H[0]=MAGMA_C_ONE + MAGMA_C_ONE;
    h_phc_coeffs_H[1]= MAGMA_C_ZERO;
    h_phc_coeffs_H[2]=q1;
    h_phc_coeffs_H[3]=p1 - q1;
    h_phc_coeffs_H[4]=q4;
    h_phc_coeffs_H[5]=p4 - q4;
    h_phc_coeffs_H[6]=q7;
    h_phc_coeffs_H[7]=p7 - q7;
    h_phc_coeffs_H[8]=-2*q19;
    h_phc_coeffs_H[9]=2*q19 - 2*p19;
    h_phc_coeffs_H[10]=q2;
    h_phc_coeffs_H[11]=p2 - q2;
    h_phc_coeffs_H[12]=q5;
    h_phc_coeffs_H[13]=p5 - q5;
    h_phc_coeffs_H[14]=q8;
    h_phc_coeffs_H[15]=p8 - q8;
    h_phc_coeffs_H[16]=-2*q20;
    h_phc_coeffs_H[17]=2*q20 - 2*p20;
    h_phc_coeffs_H[18]=q10;
    h_phc_coeffs_H[19]=p10 - q10;
    h_phc_coeffs_H[20]=q13;
    h_phc_coeffs_H[21]=p13 - q13;
    h_phc_coeffs_H[22]=q3;
    h_phc_coeffs_H[23]=p3 - q3;
    h_phc_coeffs_H[24]=q16;
    h_phc_coeffs_H[25]=p16 - q16;
    h_phc_coeffs_H[26]=-2*q21;
    h_phc_coeffs_H[27]=2*q21 - 2*p21;
    h_phc_coeffs_H[28]=q11;
    h_phc_coeffs_H[29]=p11 - q11;
    h_phc_coeffs_H[30]=q14;
    h_phc_coeffs_H[31]=p14 - q14;
    h_phc_coeffs_H[32]=q6;
    h_phc_coeffs_H[33]=p6 - q6;
    h_phc_coeffs_H[34]=q17;
    h_phc_coeffs_H[35]=p17 - q17;
    h_phc_coeffs_H[36]=-2*q22;
    h_phc_coeffs_H[37]=2*q22 - 2*p22;
    h_phc_coeffs_H[38]=q12;
    h_phc_coeffs_H[39]=p12 - q12;
    h_phc_coeffs_H[40]=-2*q23;
    h_phc_coeffs_H[41]=2*q23 - 2*p23;
    h_phc_coeffs_H[42]=q15;
    h_phc_coeffs_H[43]=p15 - q15;
    h_phc_coeffs_H[44]=-2*q24;
    h_phc_coeffs_H[45]=2*q24 - 2*p24;
    h_phc_coeffs_H[46]=q9;
    h_phc_coeffs_H[47]=p9 - q9;
    h_phc_coeffs_H[48]=q18;
    h_phc_coeffs_H[49]=p18 - q18;

    h_phc_coeffs_Ht[0]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[1]=p1 - q1;
    h_phc_coeffs_Ht[2]=p4 - q4;
    h_phc_coeffs_Ht[3]=p7 - q7;
    h_phc_coeffs_Ht[4]=2*q19 - 2*p19;
    h_phc_coeffs_Ht[5]=p2 - q2;
    h_phc_coeffs_Ht[6]=p5 - q5;
    h_phc_coeffs_Ht[7]=p8 - q8;
    h_phc_coeffs_Ht[8]=2*q20 - 2*p20;
    h_phc_coeffs_Ht[9]=p10 - q10;
    h_phc_coeffs_Ht[10]=p13 - q13;
    h_phc_coeffs_Ht[11]=p3 - q3;
    h_phc_coeffs_Ht[12]=p16 - q16;
    h_phc_coeffs_Ht[13]=2*q21 - 2*p21;
    h_phc_coeffs_Ht[14]=p11 - q11;
    h_phc_coeffs_Ht[15]=p14 - q14;
    h_phc_coeffs_Ht[16]=p6 - q6;
    h_phc_coeffs_Ht[17]=p17 - q17;
    h_phc_coeffs_Ht[18]=2*q22 - 2*p22;
    h_phc_coeffs_Ht[19]=p12 - q12;
    h_phc_coeffs_Ht[20]=2*q23 - 2*p23;
    h_phc_coeffs_Ht[21]=p15 - q15;
    h_phc_coeffs_Ht[22]=2*q24 - 2*p24;
    h_phc_coeffs_Ht[23]=p9 - q9;
    h_phc_coeffs_Ht[24]=p18 - q18;
  }
} // end of namespace

#endif
