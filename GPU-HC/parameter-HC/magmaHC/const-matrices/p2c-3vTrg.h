#ifndef p2c_3vTrg_h
#define p2c_3vTrg_h
// =======================================================================
//
// Modifications
//    Chien  21-12-30:   Initially Created
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

  void p2c_3vTrg(
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
    magmaFloatComplex p25 = h_targetParams[24];
    magmaFloatComplex p26 = h_targetParams[25];
    magmaFloatComplex p27 = h_targetParams[26];
    magmaFloatComplex p28 = h_targetParams[27];
    magmaFloatComplex p29 = h_targetParams[28];
    magmaFloatComplex p30 = h_targetParams[29];
    magmaFloatComplex p31 = h_targetParams[30];
    magmaFloatComplex p32 = h_targetParams[31];
    magmaFloatComplex p33 = h_targetParams[32];

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
    magmaFloatComplex q25 = h_startParams[24];
    magmaFloatComplex q26 = h_startParams[25];
    magmaFloatComplex q27 = h_startParams[26];
    magmaFloatComplex q28 = h_startParams[27];
    magmaFloatComplex q29 = h_startParams[28];
    magmaFloatComplex q30 = h_startParams[29];
    magmaFloatComplex q31 = h_startParams[30];
    magmaFloatComplex q32 = h_startParams[31];
    magmaFloatComplex q33 = h_startParams[32];

    h_phc_coeffs_H[0]=MAGMA_C_ONE+MAGMA_C_ONE;
    h_phc_coeffs_H[1]= MAGMA_C_ZERO;
    h_phc_coeffs_H[2]=q1;
    h_phc_coeffs_H[3]=p1 - q1;
    h_phc_coeffs_H[4]=q2;
    h_phc_coeffs_H[5]=p2 - q2;
    h_phc_coeffs_H[6]=q19;
    h_phc_coeffs_H[7]=p19 - q19;
    h_phc_coeffs_H[8]=q20;
    h_phc_coeffs_H[9]=p20 - q20;
    h_phc_coeffs_H[10]=q3;
    h_phc_coeffs_H[11]=p3 - q3;
    h_phc_coeffs_H[12]=q21;
    h_phc_coeffs_H[13]=p21 - q21;
    h_phc_coeffs_H[14]=-2*q28;
    h_phc_coeffs_H[15]=2*q28 - 2*p28;
    h_phc_coeffs_H[16]=q4;
    h_phc_coeffs_H[17]=p4 - q4;
    h_phc_coeffs_H[18]=q5;
    h_phc_coeffs_H[19]=p5 - q5;
    h_phc_coeffs_H[20]=q22;
    h_phc_coeffs_H[21]=p22 - q22;
    h_phc_coeffs_H[22]=q23;
    h_phc_coeffs_H[23]=p23 - q23;
    h_phc_coeffs_H[24]=q6;
    h_phc_coeffs_H[25]=p6 - q6;
    h_phc_coeffs_H[26]=q24;
    h_phc_coeffs_H[27]=p24 - q24;
    h_phc_coeffs_H[28]=-2*q29;
    h_phc_coeffs_H[29]=2*q29 - 2*p29;
    h_phc_coeffs_H[30]=q10;
    h_phc_coeffs_H[31]=p10 - q10;
    h_phc_coeffs_H[32]=q11;
    h_phc_coeffs_H[33]=p11 - q11;
    h_phc_coeffs_H[34]=q7;
    h_phc_coeffs_H[35]=p7 - q7;
    h_phc_coeffs_H[36]=q12;
    h_phc_coeffs_H[37]=p12 - q12;
    h_phc_coeffs_H[38]=-2*q30;
    h_phc_coeffs_H[39]=2*q30 - 2*p30;
    h_phc_coeffs_H[40]=q13;
    h_phc_coeffs_H[41]=p13 - q13;
    h_phc_coeffs_H[42]=q14;
    h_phc_coeffs_H[43]=p14 - q14;
    h_phc_coeffs_H[44]=q8;
    h_phc_coeffs_H[45]=p8 - q8;
    h_phc_coeffs_H[46]=q15;
    h_phc_coeffs_H[47]=p15 - q15;
    h_phc_coeffs_H[48]=-2*q31;
    h_phc_coeffs_H[49]=2*q31 - 2*p31;
    h_phc_coeffs_H[50]=q16;
    h_phc_coeffs_H[51]=p16 - q16;
    h_phc_coeffs_H[52]=q25;
    h_phc_coeffs_H[53]=p25 - q25;
    h_phc_coeffs_H[54]=-2*q32;
    h_phc_coeffs_H[55]=2*q32 - 2*p32;
    h_phc_coeffs_H[56]=q17;
    h_phc_coeffs_H[57]=p17 - q17;
    h_phc_coeffs_H[58]=q26;
    h_phc_coeffs_H[59]=p26 - q26;
    h_phc_coeffs_H[60]=-2*q33;
    h_phc_coeffs_H[61]=2*q33 - 2*p33;
    h_phc_coeffs_H[62]=q9;
    h_phc_coeffs_H[63]=p9 - q9;
    h_phc_coeffs_H[64]=q18;
    h_phc_coeffs_H[65]=p18 - q18;
    h_phc_coeffs_H[66]=q27;
    h_phc_coeffs_H[67]=p27 - q27;

    h_phc_coeffs_Ht[0]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[1]=p1 - q1;
    h_phc_coeffs_Ht[2]=p2 - q2;
    h_phc_coeffs_Ht[3]=p19 - q19;
    h_phc_coeffs_Ht[4]=p20 - q20;
    h_phc_coeffs_Ht[5]=p3 - q3;
    h_phc_coeffs_Ht[6]=p21 - q21;
    h_phc_coeffs_Ht[7]=2*q28 - 2*p28;
    h_phc_coeffs_Ht[8]=p4 - q4;
    h_phc_coeffs_Ht[9]=p5 - q5;
    h_phc_coeffs_Ht[10]=p22 - q22;
    h_phc_coeffs_Ht[11]=p23 - q23;
    h_phc_coeffs_Ht[12]=p6 - q6;
    h_phc_coeffs_Ht[13]=p24 - q24;
    h_phc_coeffs_Ht[14]=2*q29 - 2*p29;
    h_phc_coeffs_Ht[15]=p10 - q10;
    h_phc_coeffs_Ht[16]=p11 - q11;
    h_phc_coeffs_Ht[17]=p7 - q7;
    h_phc_coeffs_Ht[18]=p12 - q12;
    h_phc_coeffs_Ht[19]=2*q30 - 2*p30;
    h_phc_coeffs_Ht[20]=p13 - q13;
    h_phc_coeffs_Ht[21]=p14 - q14;
    h_phc_coeffs_Ht[22]=p8 - q8;
    h_phc_coeffs_Ht[23]=p15 - q15;
    h_phc_coeffs_Ht[24]=2*q31 - 2*p31;
    h_phc_coeffs_Ht[25]=p16 - q16;
    h_phc_coeffs_Ht[26]=p25 - q25;
    h_phc_coeffs_Ht[27]=2*q32 - 2*p32;
    h_phc_coeffs_Ht[28]=p17 - q17;
    h_phc_coeffs_Ht[29]=p26 - q26;
    h_phc_coeffs_Ht[30]=2*q33 - 2*p33;
    h_phc_coeffs_Ht[31]=p9 - q9;
    h_phc_coeffs_Ht[32]=p18 - q18;
    h_phc_coeffs_Ht[33]=p27 - q27;
  }
} // end of namespace

#endif
