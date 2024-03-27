#ifndef p2c_trifocal_2op1p_30x30_h
#define p2c_trifocal_2op1p_30x30_h
// =======================================================================
//
// Modifications
//    Chien  23-06-29:   Shifted from other private repository
//
// =======================================================================
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <chrono>

#include "flops.h"
#include "magma_v2.h"
#include "magma_lapack.h"
#include "magma_internal.h"
#undef max
#undef min

namespace magmaHCWrapper {

  void p2c_trifocal_2op1p_30x30(
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

    h_phc_coeffs_Hx[0]=q1*q31;
    h_phc_coeffs_Hx[1]=p1*q31 + p31*q1 - 2*q1*q31;
    h_phc_coeffs_Hx[2]=p1*p31 - p1*q31 - p31*q1 + q1*q31;
    h_phc_coeffs_Hx[3]=q1*q32;
    h_phc_coeffs_Hx[4]=p1*q32 + p32*q1 - 2*q1*q32;
    h_phc_coeffs_Hx[5]=p1*p32 - p1*q32 - p32*q1 + q1*q32;
    h_phc_coeffs_Hx[6]=q10;
    h_phc_coeffs_Hx[7]=p10 - q10;
    h_phc_coeffs_Hx[8] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[9]=q11;
    h_phc_coeffs_Hx[10]=p11 - q11;
    h_phc_coeffs_Hx[11] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[12]=q12;
    h_phc_coeffs_Hx[13]=p12 - q12;
    h_phc_coeffs_Hx[14] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[15]=q13;
    h_phc_coeffs_Hx[16]=p13 - q13;
    h_phc_coeffs_Hx[17] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[18]=q14;
    h_phc_coeffs_Hx[19]=p14 - q14;
    h_phc_coeffs_Hx[20] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[21]=q15;
    h_phc_coeffs_Hx[22]=p15 - q15;
    h_phc_coeffs_Hx[23] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[24]=q16;
    h_phc_coeffs_Hx[25]=p16 - q16;
    h_phc_coeffs_Hx[26] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[27]=q17;
    h_phc_coeffs_Hx[28]=p17 - q17;
    h_phc_coeffs_Hx[29] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[30]=q18;
    h_phc_coeffs_Hx[31]=p18 - q18;
    h_phc_coeffs_Hx[32] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[33]=q19;
    h_phc_coeffs_Hx[34]=p19 - q19;
    h_phc_coeffs_Hx[35] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[36]=q2*q31;
    h_phc_coeffs_Hx[37]=p2*q31 + p31*q2 - 2*q2*q31;
    h_phc_coeffs_Hx[38]=p2*p31 - p2*q31 - p31*q2 + q2*q31;
    h_phc_coeffs_Hx[39]=q2*q32;
    h_phc_coeffs_Hx[40]=p2*q32 + p32*q2 - 2*q2*q32;
    h_phc_coeffs_Hx[41]=p2*p32 - p2*q32 - p32*q2 + q2*q32;
    h_phc_coeffs_Hx[42]=q20;
    h_phc_coeffs_Hx[43]=p20 - q20;
    h_phc_coeffs_Hx[44] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[45]=q21;
    h_phc_coeffs_Hx[46]=p21 - q21;
    h_phc_coeffs_Hx[47] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[48]=q22;
    h_phc_coeffs_Hx[49]=p22 - q22;
    h_phc_coeffs_Hx[50] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[51]=q23;
    h_phc_coeffs_Hx[52]=p23 - q23;
    h_phc_coeffs_Hx[53] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[54]=q24;
    h_phc_coeffs_Hx[55]=p24 - q24;
    h_phc_coeffs_Hx[56] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[57]=q25;
    h_phc_coeffs_Hx[58]=p25 - q25;
    h_phc_coeffs_Hx[59] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[60]=q26;
    h_phc_coeffs_Hx[61]=p26 - q26;
    h_phc_coeffs_Hx[62] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[63]=q27;
    h_phc_coeffs_Hx[64]=p27 - q27;
    h_phc_coeffs_Hx[65] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[66]=q28;
    h_phc_coeffs_Hx[67]=p28 - q28;
    h_phc_coeffs_Hx[68] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[69]=q29;
    h_phc_coeffs_Hx[70]=p29 - q29;
    h_phc_coeffs_Hx[71] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[72]=q3;
    h_phc_coeffs_Hx[73]=p3 - q3;
    h_phc_coeffs_Hx[74] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[75]=q30;
    h_phc_coeffs_Hx[76]=p30 - q30;
    h_phc_coeffs_Hx[77] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[78]=q31;
    h_phc_coeffs_Hx[79]=p31 - q31;
    h_phc_coeffs_Hx[80] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[81]=q32;
    h_phc_coeffs_Hx[82]=p32 - q32;
    h_phc_coeffs_Hx[83] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[84]=q33;
    h_phc_coeffs_Hx[85]=p33 - q33;
    h_phc_coeffs_Hx[86] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[87]=q4;
    h_phc_coeffs_Hx[88]=p4 - q4;
    h_phc_coeffs_Hx[89] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[90]=q5;
    h_phc_coeffs_Hx[91]=p5 - q5;
    h_phc_coeffs_Hx[92] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[93]=q6;
    h_phc_coeffs_Hx[94]=p6 - q6;
    h_phc_coeffs_Hx[95] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[96]=q7;
    h_phc_coeffs_Hx[97]=p7 - q7;
    h_phc_coeffs_Hx[98] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[99]=q7*q33;
    h_phc_coeffs_Hx[100]=p7*q33 + p33*q7 - 2*q7*q33;
    h_phc_coeffs_Hx[101]=p7*p33 - p7*q33 - p33*q7 + q7*q33;
    h_phc_coeffs_Hx[102]=q8;
    h_phc_coeffs_Hx[103]=p8 - q8;
    h_phc_coeffs_Hx[104] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[105]=q8*q33;
    h_phc_coeffs_Hx[106]=p8*q33 + p33*q8 - 2*q8*q33;
    h_phc_coeffs_Hx[107]=p8*p33 - p8*q33 - p33*q8 + q8*q33;
    h_phc_coeffs_Hx[108]=q9;
    h_phc_coeffs_Hx[109]=p9 - q9;
    h_phc_coeffs_Hx[110] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[111] = MAGMA_C_ONE;
    h_phc_coeffs_Hx[112] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[113] = MAGMA_C_ZERO;

    h_phc_coeffs_Ht[0]=p1*q31 + p31*q1 - 2*q1*q31;
    h_phc_coeffs_Ht[1]=2*p1*p31 - 2*p1*q31 - 2*p31*q1 + 2*q1*q31;
    h_phc_coeffs_Ht[2]=p1*q32 + p32*q1 - 2*q1*q32;
    h_phc_coeffs_Ht[3]=2*p1*p32 - 2*p1*q32 - 2*p32*q1 + 2*q1*q32;
    h_phc_coeffs_Ht[4]=p10 - q10;
    h_phc_coeffs_Ht[5]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[6]=p11 - q11;
    h_phc_coeffs_Ht[7]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[8]=p12 - q12;
    h_phc_coeffs_Ht[9]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[10]=p13 - q13;
    h_phc_coeffs_Ht[11]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[12]=p14 - q14;
    h_phc_coeffs_Ht[13]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[14]=p15 - q15;
    h_phc_coeffs_Ht[15]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[16]=p16 - q16;
    h_phc_coeffs_Ht[17]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[18]=p17 - q17;
    h_phc_coeffs_Ht[19]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[20]=p18 - q18;
    h_phc_coeffs_Ht[21]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[22]=p19 - q19;
    h_phc_coeffs_Ht[23]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[24]=p2*q31 + p31*q2 - 2*q2*q31;
    h_phc_coeffs_Ht[25]=2*p2*p31 - 2*p2*q31 - 2*p31*q2 + 2*q2*q31;
    h_phc_coeffs_Ht[26]=p2*q32 + p32*q2 - 2*q2*q32;
    h_phc_coeffs_Ht[27]=2*p2*p32 - 2*p2*q32 - 2*p32*q2 + 2*q2*q32;
    h_phc_coeffs_Ht[28]=p20 - q20;
    h_phc_coeffs_Ht[29]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[30]=p21 - q21;
    h_phc_coeffs_Ht[31]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[32]=p22 - q22;
    h_phc_coeffs_Ht[33]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[34]=p23 - q23;
    h_phc_coeffs_Ht[35]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[36]=p24 - q24;
    h_phc_coeffs_Ht[37]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[38]=p25 - q25;
    h_phc_coeffs_Ht[39]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[40]=p26 - q26;
    h_phc_coeffs_Ht[41]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[42]=p27 - q27;
    h_phc_coeffs_Ht[43]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[44]=p28 - q28;
    h_phc_coeffs_Ht[45]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[46]=p29 - q29;
    h_phc_coeffs_Ht[47]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[48]=p3 - q3;
    h_phc_coeffs_Ht[49]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[50]=p30 - q30;
    h_phc_coeffs_Ht[51]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[52]=p31 - q31;
    h_phc_coeffs_Ht[53]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[54]=p32 - q32;
    h_phc_coeffs_Ht[55]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[56]=p33 - q33;
    h_phc_coeffs_Ht[57]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[58]=p4 - q4;
    h_phc_coeffs_Ht[59]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[60]=p5 - q5;
    h_phc_coeffs_Ht[61]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[62]=p6 - q6;
    h_phc_coeffs_Ht[63]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[64]=p7 - q7;
    h_phc_coeffs_Ht[65]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[66]=p7*q33 + p33*q7 - 2*q7*q33;
    h_phc_coeffs_Ht[67]=2*p7*p33 - 2*p7*q33 - 2*p33*q7 + 2*q7*q33;
    h_phc_coeffs_Ht[68]=p8 - q8;
    h_phc_coeffs_Ht[69]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[70]=p8*q33 + p33*q8 - 2*q8*q33;
    h_phc_coeffs_Ht[71]=2*p8*p33 - 2*p8*q33 - 2*p33*q8 + 2*q8*q33;
    h_phc_coeffs_Ht[72]=p9 - q9;
    h_phc_coeffs_Ht[73]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[74]= MAGMA_C_ONE;
    h_phc_coeffs_Ht[75]= MAGMA_C_ZERO;
  }
} // end of namespace

#endif
