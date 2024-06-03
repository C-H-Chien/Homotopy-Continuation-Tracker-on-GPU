#ifndef P2C_3PT_REL_POS_HOMOG_H
#define P2C_3PT_REL_POS_HOMOG_H
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

  void p2c_3pt_rel_pos_homog(
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

    h_phc_coeffs_Hx[0]=q1*q10;
    h_phc_coeffs_Hx[1]=p1*q10 + p10*q1 - 2*q1*q10;
    h_phc_coeffs_Hx[2]=p1*p10 - p1*q10 - p10*q1 + q1*q10;
    h_phc_coeffs_Hx[3]=q1*q13;
    h_phc_coeffs_Hx[4]=p1*q13 + p13*q1 - 2*q1*q13;
    h_phc_coeffs_Hx[5]=p1*p13 - p1*q13 - p13*q1 + q1*q13;
    h_phc_coeffs_Hx[6]=q1*q16;
    h_phc_coeffs_Hx[7]=p1*q16 + p16*q1 - 2*q1*q16;
    h_phc_coeffs_Hx[8]=p1*p16 - p1*q16 - p16*q1 + q1*q16;
    h_phc_coeffs_Hx[9]=q2*q11;
    h_phc_coeffs_Hx[10]=p2*q11 + p11*q2 - 2*q2*q11;
    h_phc_coeffs_Hx[11]=p2*p11 - p2*q11 - p11*q2 + q2*q11;
    h_phc_coeffs_Hx[12]=q2*q14;
    h_phc_coeffs_Hx[13]=p2*q14 + p14*q2 - 2*q2*q14;
    h_phc_coeffs_Hx[14]=p2*p14 - p2*q14 - p14*q2 + q2*q14;
    h_phc_coeffs_Hx[15]=q2*q17;
    h_phc_coeffs_Hx[16]=p2*q17 + p17*q2 - 2*q2*q17;
    h_phc_coeffs_Hx[17]=p2*p17 - p2*q17 - p17*q2 + q2*q17;
    h_phc_coeffs_Hx[18]=q3*q12;
    h_phc_coeffs_Hx[19]=p3*q12 + p12*q3 - 2*q3*q12;
    h_phc_coeffs_Hx[20]=p3*p12 - p3*q12 - p12*q3 + q3*q12;
    h_phc_coeffs_Hx[21]=q3*q15;
    h_phc_coeffs_Hx[22]=p3*q15 + p15*q3 - 2*q3*q15;
    h_phc_coeffs_Hx[23]=p3*p15 - p3*q15 - p15*q3 + q3*q15;
    h_phc_coeffs_Hx[24]=q3*q18;
    h_phc_coeffs_Hx[25]=p3*q18 + p18*q3 - 2*q3*q18;
    h_phc_coeffs_Hx[26]=p3*p18 - p3*q18 - p18*q3 + q3*q18;
    h_phc_coeffs_Hx[27]=q4*q10;
    h_phc_coeffs_Hx[28]=p4*q10 + p10*q4 - 2*q4*q10;
    h_phc_coeffs_Hx[29]=p4*p10 - p4*q10 - p10*q4 + q4*q10;
    h_phc_coeffs_Hx[30]=q4*q13;
    h_phc_coeffs_Hx[31]=p4*q13 + p13*q4 - 2*q4*q13;
    h_phc_coeffs_Hx[32]=p4*p13 - p4*q13 - p13*q4 + q4*q13;
    h_phc_coeffs_Hx[33]=q4*q16;
    h_phc_coeffs_Hx[34]=p4*q16 + p16*q4 - 2*q4*q16;
    h_phc_coeffs_Hx[35]=p4*p16 - p4*q16 - p16*q4 + q4*q16;
    h_phc_coeffs_Hx[36]=q5*q11;
    h_phc_coeffs_Hx[37]=p5*q11 + p11*q5 - 2*q5*q11;
    h_phc_coeffs_Hx[38]=p5*p11 - p5*q11 - p11*q5 + q5*q11;
    h_phc_coeffs_Hx[39]=q5*q14;
    h_phc_coeffs_Hx[40]=p5*q14 + p14*q5 - 2*q5*q14;
    h_phc_coeffs_Hx[41]=p5*p14 - p5*q14 - p14*q5 + q5*q14;
    h_phc_coeffs_Hx[42]=q5*q17;
    h_phc_coeffs_Hx[43]=p5*q17 + p17*q5 - 2*q5*q17;
    h_phc_coeffs_Hx[44]=p5*p17 - p5*q17 - p17*q5 + q5*q17;
    h_phc_coeffs_Hx[45]=q6*q12;
    h_phc_coeffs_Hx[46]=p6*q12 + p12*q6 - 2*q6*q12;
    h_phc_coeffs_Hx[47]=p6*p12 - p6*q12 - p12*q6 + q6*q12;
    h_phc_coeffs_Hx[48]=q6*q15;
    h_phc_coeffs_Hx[49]=p6*q15 + p15*q6 - 2*q6*q15;
    h_phc_coeffs_Hx[50]=p6*p15 - p6*q15 - p15*q6 + q6*q15;
    h_phc_coeffs_Hx[51]=q6*q18;
    h_phc_coeffs_Hx[52]=p6*q18 + p18*q6 - 2*q6*q18;
    h_phc_coeffs_Hx[53]=p6*p18 - p6*q18 - p18*q6 + q6*q18;
    h_phc_coeffs_Hx[54]=q7*q10;
    h_phc_coeffs_Hx[55]=p7*q10 + p10*q7 - 2*q7*q10;
    h_phc_coeffs_Hx[56]=p7*p10 - p7*q10 - p10*q7 + q7*q10;
    h_phc_coeffs_Hx[57]=q7*q13;
    h_phc_coeffs_Hx[58]=p7*q13 + p13*q7 - 2*q7*q13;
    h_phc_coeffs_Hx[59]=p7*p13 - p7*q13 - p13*q7 + q7*q13;
    h_phc_coeffs_Hx[60]=q7*q16;
    h_phc_coeffs_Hx[61]=p7*q16 + p16*q7 - 2*q7*q16;
    h_phc_coeffs_Hx[62]=p7*p16 - p7*q16 - p16*q7 + q7*q16;
    h_phc_coeffs_Hx[63]=q8*q11;
    h_phc_coeffs_Hx[64]=p8*q11 + p11*q8 - 2*q8*q11;
    h_phc_coeffs_Hx[65]=p8*p11 - p8*q11 - p11*q8 + q8*q11;
    h_phc_coeffs_Hx[66]=q8*q14;
    h_phc_coeffs_Hx[67]=p8*q14 + p14*q8 - 2*q8*q14;
    h_phc_coeffs_Hx[68]=p8*p14 - p8*q14 - p14*q8 + q8*q14;
    h_phc_coeffs_Hx[69]=q8*q17;
    h_phc_coeffs_Hx[70]=p8*q17 + p17*q8 - 2*q8*q17;
    h_phc_coeffs_Hx[71]=p8*p17 - p8*q17 - p17*q8 + q8*q17;
    h_phc_coeffs_Hx[72]=q9*q12;
    h_phc_coeffs_Hx[73]=p9*q12 + p12*q9 - 2*q9*q12;
    h_phc_coeffs_Hx[74]=p9*p12 - p9*q12 - p12*q9 + q9*q12;
    h_phc_coeffs_Hx[75]=q9*q15;
    h_phc_coeffs_Hx[76]=p9*q15 + p15*q9 - 2*q9*q15;
    h_phc_coeffs_Hx[77]=p9*p15 - p9*q15 - p15*q9 + q9*q15;
    h_phc_coeffs_Hx[78]=q9*q18;
    h_phc_coeffs_Hx[79]=p9*q18 + p18*q9 - 2*q9*q18;
    h_phc_coeffs_Hx[80]=p9*p18 - p9*q18 - p18*q9 + q9*q18;
    h_phc_coeffs_Hx[81]=MAGMA_C_ONE;
    h_phc_coeffs_Hx[82]=MAGMA_C_ZERO;
    h_phc_coeffs_Hx[83]=MAGMA_C_ZERO;

    h_phc_coeffs_Ht[0]=p1*q10 + p10*q1 - 2*q1*q10;
    h_phc_coeffs_Ht[1]=2*p1*p10 - 2*p1*q10 - 2*p10*q1 + 2*q1*q10;
    h_phc_coeffs_Ht[2]=p1*q13 + p13*q1 - 2*q1*q13;
    h_phc_coeffs_Ht[3]=2*p1*p13 - 2*p1*q13 - 2*p13*q1 + 2*q1*q13;
    h_phc_coeffs_Ht[4]=p1*q16 + p16*q1 - 2*q1*q16;
    h_phc_coeffs_Ht[5]=2*p1*p16 - 2*p1*q16 - 2*p16*q1 + 2*q1*q16;
    h_phc_coeffs_Ht[6]=p2*q11 + p11*q2 - 2*q2*q11;
    h_phc_coeffs_Ht[7]=2*p2*p11 - 2*p2*q11 - 2*p11*q2 + 2*q2*q11;
    h_phc_coeffs_Ht[8]=p2*q14 + p14*q2 - 2*q2*q14;
    h_phc_coeffs_Ht[9]=2*p2*p14 - 2*p2*q14 - 2*p14*q2 + 2*q2*q14;
    h_phc_coeffs_Ht[10]=p2*q17 + p17*q2 - 2*q2*q17;
    h_phc_coeffs_Ht[11]=2*p2*p17 - 2*p2*q17 - 2*p17*q2 + 2*q2*q17;
    h_phc_coeffs_Ht[12]=p3*q12 + p12*q3 - 2*q3*q12;
    h_phc_coeffs_Ht[13]=2*p3*p12 - 2*p3*q12 - 2*p12*q3 + 2*q3*q12;
    h_phc_coeffs_Ht[14]=p3*q15 + p15*q3 - 2*q3*q15;
    h_phc_coeffs_Ht[15]=2*p3*p15 - 2*p3*q15 - 2*p15*q3 + 2*q3*q15;
    h_phc_coeffs_Ht[16]=p3*q18 + p18*q3 - 2*q3*q18;
    h_phc_coeffs_Ht[17]=2*p3*p18 - 2*p3*q18 - 2*p18*q3 + 2*q3*q18;
    h_phc_coeffs_Ht[18]=p4*q10 + p10*q4 - 2*q4*q10;
    h_phc_coeffs_Ht[19]=2*p4*p10 - 2*p4*q10 - 2*p10*q4 + 2*q4*q10;
    h_phc_coeffs_Ht[20]=p4*q13 + p13*q4 - 2*q4*q13;
    h_phc_coeffs_Ht[21]=2*p4*p13 - 2*p4*q13 - 2*p13*q4 + 2*q4*q13;
    h_phc_coeffs_Ht[22]=p4*q16 + p16*q4 - 2*q4*q16;
    h_phc_coeffs_Ht[23]=2*p4*p16 - 2*p4*q16 - 2*p16*q4 + 2*q4*q16;
    h_phc_coeffs_Ht[24]=p5*q11 + p11*q5 - 2*q5*q11;
    h_phc_coeffs_Ht[25]=2*p5*p11 - 2*p5*q11 - 2*p11*q5 + 2*q5*q11;
    h_phc_coeffs_Ht[26]=p5*q14 + p14*q5 - 2*q5*q14;
    h_phc_coeffs_Ht[27]=2*p5*p14 - 2*p5*q14 - 2*p14*q5 + 2*q5*q14;
    h_phc_coeffs_Ht[28]=p5*q17 + p17*q5 - 2*q5*q17;
    h_phc_coeffs_Ht[29]=2*p5*p17 - 2*p5*q17 - 2*p17*q5 + 2*q5*q17;
    h_phc_coeffs_Ht[30]=p6*q12 + p12*q6 - 2*q6*q12;
    h_phc_coeffs_Ht[31]=2*p6*p12 - 2*p6*q12 - 2*p12*q6 + 2*q6*q12;
    h_phc_coeffs_Ht[32]=p6*q15 + p15*q6 - 2*q6*q15;
    h_phc_coeffs_Ht[33]=2*p6*p15 - 2*p6*q15 - 2*p15*q6 + 2*q6*q15;
    h_phc_coeffs_Ht[34]=p6*q18 + p18*q6 - 2*q6*q18;
    h_phc_coeffs_Ht[35]=2*p6*p18 - 2*p6*q18 - 2*p18*q6 + 2*q6*q18;
    h_phc_coeffs_Ht[36]=p7*q10 + p10*q7 - 2*q7*q10;
    h_phc_coeffs_Ht[37]=2*p7*p10 - 2*p7*q10 - 2*p10*q7 + 2*q7*q10;
    h_phc_coeffs_Ht[38]=p7*q13 + p13*q7 - 2*q7*q13;
    h_phc_coeffs_Ht[39]=2*p7*p13 - 2*p7*q13 - 2*p13*q7 + 2*q7*q13;
    h_phc_coeffs_Ht[40]=p7*q16 + p16*q7 - 2*q7*q16;
    h_phc_coeffs_Ht[41]=2*p7*p16 - 2*p7*q16 - 2*p16*q7 + 2*q7*q16;
    h_phc_coeffs_Ht[42]=p8*q11 + p11*q8 - 2*q8*q11;
    h_phc_coeffs_Ht[43]=2*p8*p11 - 2*p8*q11 - 2*p11*q8 + 2*q8*q11;
    h_phc_coeffs_Ht[44]=p8*q14 + p14*q8 - 2*q8*q14;
    h_phc_coeffs_Ht[45]=2*p8*p14 - 2*p8*q14 - 2*p14*q8 + 2*q8*q14;
    h_phc_coeffs_Ht[46]=p8*q17 + p17*q8 - 2*q8*q17;
    h_phc_coeffs_Ht[47]=2*p8*p17 - 2*p8*q17 - 2*p17*q8 + 2*q8*q17;
    h_phc_coeffs_Ht[48]=p9*q12 + p12*q9 - 2*q9*q12;
    h_phc_coeffs_Ht[49]=2*p9*p12 - 2*p9*q12 - 2*p12*q9 + 2*q9*q12;
    h_phc_coeffs_Ht[50]=p9*q15 + p15*q9 - 2*q9*q15;
    h_phc_coeffs_Ht[51]=2*p9*p15 - 2*p9*q15 - 2*p15*q9 + 2*q9*q15;
    h_phc_coeffs_Ht[52]=p9*q18 + p18*q9 - 2*q9*q18;
    h_phc_coeffs_Ht[53]=2*p9*p18 - 2*p9*q18 - 2*p18*q9 + 2*q9*q18;
    h_phc_coeffs_Ht[54]=MAGMA_C_ONE;
    h_phc_coeffs_Ht[55]=MAGMA_C_ZERO;

  }
} // end of namespace

#endif
