#ifndef p2c_3pt_rel_pose_w_homo_constraint_h
#define p2c_3pt_rel_pose_w_homo_constraint_h
// =======================================================================
// Modifications
//    Chien  22-01-02:   Initially Created
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

  void p2c_3pt_rel_pose_w_homo_constraint(
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

    h_phc_coeffs_H[0]=q1*q16;
    h_phc_coeffs_H[1]=p1*q16 + p16*q1 - 2*q1*q16;
    h_phc_coeffs_H[2]=p1*p16 - p1*q16 - p16*q1 + q1*q16;
    h_phc_coeffs_H[3]=-q1*q13;
    h_phc_coeffs_H[4]=2*q1*q13 - p13*q1 - p1*q13;
    h_phc_coeffs_H[5]=p1*q13 - p1*p13 + p13*q1 - q1*q13;
    h_phc_coeffs_H[6]=q4*q16;
    h_phc_coeffs_H[7]=p4*q16 + p16*q4 - 2*q4*q16;
    h_phc_coeffs_H[8]=p4*p16 - p4*q16 - p16*q4 + q4*q16;
    h_phc_coeffs_H[9]=-q4*q13;
    h_phc_coeffs_H[10]=2*q4*q13 - p13*q4 - p4*q13;
    h_phc_coeffs_H[11]=p4*q13 - p4*p13 + p13*q4 - q4*q13;
    h_phc_coeffs_H[12]=q7*q16;
    h_phc_coeffs_H[13]=p7*q16 + p16*q7 - 2*q7*q16;
    h_phc_coeffs_H[14]=p7*p16 - p7*q16 - p16*q7 + q7*q16;
    h_phc_coeffs_H[15]=-q7*q13;
    h_phc_coeffs_H[16]=2*q7*q13 - p13*q7 - p7*q13;
    h_phc_coeffs_H[17]=p7*q13 - p7*p13 + p13*q7 - q7*q13;
    h_phc_coeffs_H[18]=-q4*q16;
    h_phc_coeffs_H[19]=2*q4*q16 - p16*q4 - p4*q16;
    h_phc_coeffs_H[20]=p4*q16 - p4*p16 + p16*q4 - q4*q16;
    h_phc_coeffs_H[21]=-q1*q16;
    h_phc_coeffs_H[22]=2*q1*q16 - p16*q1 - p1*q16;
    h_phc_coeffs_H[23]=p1*q16 - p1*p16 + p16*q1 - q1*q16;
    h_phc_coeffs_H[24]=q7*q13;
    h_phc_coeffs_H[25]=p7*q13 + p13*q7 - 2*q7*q13;
    h_phc_coeffs_H[26]=p7*p13 - p7*q13 - p13*q7 + q7*q13;
    h_phc_coeffs_H[27]=q1*q10;
    h_phc_coeffs_H[28]=p1*q10 + p10*q1 - 2*q1*q10;
    h_phc_coeffs_H[29]=p1*p10 - p1*q10 - p10*q1 + q1*q10;
    h_phc_coeffs_H[30]=q4*q10;
    h_phc_coeffs_H[31]=p4*q10 + p10*q4 - 2*q4*q10;
    h_phc_coeffs_H[32]=p4*p10 - p4*q10 - p10*q4 + q4*q10;
    h_phc_coeffs_H[33]=-q7*q16;
    h_phc_coeffs_H[34]=2*q7*q16 - p16*q7 - p7*q16;
    h_phc_coeffs_H[35]=p7*q16 - p7*p16 + p16*q7 - q7*q16;
    h_phc_coeffs_H[36]=q7*q10;
    h_phc_coeffs_H[37]=p7*q10 + p10*q7 - 2*q7*q10;
    h_phc_coeffs_H[38]=p7*p10 - p7*q10 - p10*q7 + q7*q10;
    h_phc_coeffs_H[39]=-q7*q10;
    h_phc_coeffs_H[40]=2*q7*q10 - p10*q7 - p7*q10;
    h_phc_coeffs_H[41]=p7*q10 - p7*p10 + p10*q7 - q7*q10;
    h_phc_coeffs_H[42]=q2*q17;
    h_phc_coeffs_H[43]=p2*q17 + p17*q2 - 2*q2*q17;
    h_phc_coeffs_H[44]=p2*p17 - p2*q17 - p17*q2 + q2*q17;
    h_phc_coeffs_H[45]=-q2*q14;
    h_phc_coeffs_H[46]=2*q2*q14 - p14*q2 - p2*q14;
    h_phc_coeffs_H[47]=p2*q14 - p2*p14 + p14*q2 - q2*q14;
    h_phc_coeffs_H[48]=q5*q17;
    h_phc_coeffs_H[49]=p5*q17 + p17*q5 - 2*q5*q17;
    h_phc_coeffs_H[50]=p5*p17 - p5*q17 - p17*q5 + q5*q17;
    h_phc_coeffs_H[51]=-q5*q14;
    h_phc_coeffs_H[52]=2*q5*q14 - p14*q5 - p5*q14;
    h_phc_coeffs_H[53]=p5*q14 - p5*p14 + p14*q5 - q5*q14;
    h_phc_coeffs_H[54]=q8*q17;
    h_phc_coeffs_H[55]=p8*q17 + p17*q8 - 2*q8*q17;
    h_phc_coeffs_H[56]=p8*p17 - p8*q17 - p17*q8 + q8*q17;
    h_phc_coeffs_H[57]=-q8*q14;
    h_phc_coeffs_H[58]=2*q8*q14 - p14*q8 - p8*q14;
    h_phc_coeffs_H[59]=p8*q14 - p8*p14 + p14*q8 - q8*q14;
    h_phc_coeffs_H[60]=-q5*q17;
    h_phc_coeffs_H[61]=2*q5*q17 - p17*q5 - p5*q17;
    h_phc_coeffs_H[62]=p5*q17 - p5*p17 + p17*q5 - q5*q17;
    h_phc_coeffs_H[63]=-q2*q17;
    h_phc_coeffs_H[64]=2*q2*q17 - p17*q2 - p2*q17;
    h_phc_coeffs_H[65]=p2*q17 - p2*p17 + p17*q2 - q2*q17;
    h_phc_coeffs_H[66]=q8*q14;
    h_phc_coeffs_H[67]=p8*q14 + p14*q8 - 2*q8*q14;
    h_phc_coeffs_H[68]=p8*p14 - p8*q14 - p14*q8 + q8*q14;
    h_phc_coeffs_H[69]=q2*q11;
    h_phc_coeffs_H[70]=p2*q11 + p11*q2 - 2*q2*q11;
    h_phc_coeffs_H[71]=p2*p11 - p2*q11 - p11*q2 + q2*q11;
    h_phc_coeffs_H[72]=q5*q11;
    h_phc_coeffs_H[73]=p5*q11 + p11*q5 - 2*q5*q11;
    h_phc_coeffs_H[74]=p5*p11 - p5*q11 - p11*q5 + q5*q11;
    h_phc_coeffs_H[75]=-q8*q17;
    h_phc_coeffs_H[76]=2*q8*q17 - p17*q8 - p8*q17;
    h_phc_coeffs_H[77]=p8*q17 - p8*p17 + p17*q8 - q8*q17;
    h_phc_coeffs_H[78]=q8*q11;
    h_phc_coeffs_H[79]=p8*q11 + p11*q8 - 2*q8*q11;
    h_phc_coeffs_H[80]=p8*p11 - p8*q11 - p11*q8 + q8*q11;
    h_phc_coeffs_H[81]=-q8*q11;
    h_phc_coeffs_H[82]=2*q8*q11 - p11*q8 - p8*q11;
    h_phc_coeffs_H[83]=p8*q11 - p8*p11 + p11*q8 - q8*q11;
    h_phc_coeffs_H[84]=q3*q18;
    h_phc_coeffs_H[85]=p3*q18 + p18*q3 - 2*q3*q18;
    h_phc_coeffs_H[86]=p3*p18 - p3*q18 - p18*q3 + q3*q18;
    h_phc_coeffs_H[87]=-q3*q15;
    h_phc_coeffs_H[88]=2*q3*q15 - p15*q3 - p3*q15;
    h_phc_coeffs_H[89]=p3*q15 - p3*p15 + p15*q3 - q3*q15;
    h_phc_coeffs_H[90]=q6*q18;
    h_phc_coeffs_H[91]=p6*q18 + p18*q6 - 2*q6*q18;
    h_phc_coeffs_H[92]=p6*p18 - p6*q18 - p18*q6 + q6*q18;
    h_phc_coeffs_H[93]=-q6*q15;
    h_phc_coeffs_H[94]=2*q6*q15 - p15*q6 - p6*q15;
    h_phc_coeffs_H[95]=p6*q15 - p6*p15 + p15*q6 - q6*q15;
    h_phc_coeffs_H[96]=q9*q18;
    h_phc_coeffs_H[97]=p9*q18 + p18*q9 - 2*q9*q18;
    h_phc_coeffs_H[98]=p9*p18 - p9*q18 - p18*q9 + q9*q18;
    h_phc_coeffs_H[99]=-q9*q15;
    h_phc_coeffs_H[100]=2*q9*q15 - p15*q9 - p9*q15;
    h_phc_coeffs_H[101]=p9*q15 - p9*p15 + p15*q9 - q9*q15;
    h_phc_coeffs_H[102]=-q6*q18;
    h_phc_coeffs_H[103]=2*q6*q18 - p18*q6 - p6*q18;
    h_phc_coeffs_H[104]=p6*q18 - p6*p18 + p18*q6 - q6*q18;
    h_phc_coeffs_H[105]=-q3*q18;
    h_phc_coeffs_H[106]=2*q3*q18 - p18*q3 - p3*q18;
    h_phc_coeffs_H[107]=p3*q18 - p3*p18 + p18*q3 - q3*q18;
    h_phc_coeffs_H[108]=q9*q15;
    h_phc_coeffs_H[109]=p9*q15 + p15*q9 - 2*q9*q15;
    h_phc_coeffs_H[110]=p9*p15 - p9*q15 - p15*q9 + q9*q15;
    h_phc_coeffs_H[111]=q3*q12;
    h_phc_coeffs_H[112]=p3*q12 + p12*q3 - 2*q3*q12;
    h_phc_coeffs_H[113]=p3*p12 - p3*q12 - p12*q3 + q3*q12;
    h_phc_coeffs_H[114]=q6*q12;
    h_phc_coeffs_H[115]=p6*q12 + p12*q6 - 2*q6*q12;
    h_phc_coeffs_H[116]=p6*p12 - p6*q12 - p12*q6 + q6*q12;
    h_phc_coeffs_H[117]=-q9*q18;
    h_phc_coeffs_H[118]=2*q9*q18 - p18*q9 - p9*q18;
    h_phc_coeffs_H[119]=p9*q18 - p9*p18 + p18*q9 - q9*q18;
    h_phc_coeffs_H[120]=q9*q12;
    h_phc_coeffs_H[121]=p9*q12 + p12*q9 - 2*q9*q12;
    h_phc_coeffs_H[122]=p9*p12 - p9*q12 - p12*q9 + q9*q12;
    h_phc_coeffs_H[123]=-q9*q12;
    h_phc_coeffs_H[124]=2*q9*q12 - p12*q9 - p9*q12;
    h_phc_coeffs_H[125]=p9*q12 - p9*p12 + p12*q9 - q9*q12;
    h_phc_coeffs_H[126]=MAGMA_C_ONE;
    h_phc_coeffs_H[127]= MAGMA_C_ZERO;
    h_phc_coeffs_H[128]= MAGMA_C_ZERO;
    h_phc_coeffs_H[129]=MAGMA_C_NEG_ONE;
    h_phc_coeffs_H[130]= MAGMA_C_ZERO;
    h_phc_coeffs_H[131]= MAGMA_C_ZERO;

    h_phc_coeffs_Ht[0]=p1*q16 + p16*q1 - 2*q1*q16;
    h_phc_coeffs_Ht[1]=2*p1*p16 - 2*p1*q16 - 2*p16*q1 + 2*q1*q16;
    h_phc_coeffs_Ht[2]=2*q1*q13 - p13*q1 - p1*q13;
    h_phc_coeffs_Ht[3]=2*p1*q13 - 2*p1*p13 + 2*p13*q1 - 2*q1*q13;
    h_phc_coeffs_Ht[4]=p4*q16 + p16*q4 - 2*q4*q16;
    h_phc_coeffs_Ht[5]=2*p4*p16 - 2*p4*q16 - 2*p16*q4 + 2*q4*q16;
    h_phc_coeffs_Ht[6]=2*q4*q13 - p13*q4 - p4*q13;
    h_phc_coeffs_Ht[7]=2*p4*q13 - 2*p4*p13 + 2*p13*q4 - 2*q4*q13;
    h_phc_coeffs_Ht[8]=p7*q16 + p16*q7 - 2*q7*q16;
    h_phc_coeffs_Ht[9]=2*p7*p16 - 2*p7*q16 - 2*p16*q7 + 2*q7*q16;
    h_phc_coeffs_Ht[10]=2*q7*q13 - p13*q7 - p7*q13;
    h_phc_coeffs_Ht[11]=2*p7*q13 - 2*p7*p13 + 2*p13*q7 - 2*q7*q13;
    h_phc_coeffs_Ht[12]=2*q4*q16 - p16*q4 - p4*q16;
    h_phc_coeffs_Ht[13]=2*p4*q16 - 2*p4*p16 + 2*p16*q4 - 2*q4*q16;
    h_phc_coeffs_Ht[14]=2*q1*q16 - p16*q1 - p1*q16;
    h_phc_coeffs_Ht[15]=2*p1*q16 - 2*p1*p16 + 2*p16*q1 - 2*q1*q16;
    h_phc_coeffs_Ht[16]=p7*q13 + p13*q7 - 2*q7*q13;
    h_phc_coeffs_Ht[17]=2*p7*p13 - 2*p7*q13 - 2*p13*q7 + 2*q7*q13;
    h_phc_coeffs_Ht[18]=p1*q10 + p10*q1 - 2*q1*q10;
    h_phc_coeffs_Ht[19]=2*p1*p10 - 2*p1*q10 - 2*p10*q1 + 2*q1*q10;
    h_phc_coeffs_Ht[20]=p4*q10 + p10*q4 - 2*q4*q10;
    h_phc_coeffs_Ht[21]=2*p4*p10 - 2*p4*q10 - 2*p10*q4 + 2*q4*q10;
    h_phc_coeffs_Ht[22]=2*q7*q16 - p16*q7 - p7*q16;
    h_phc_coeffs_Ht[23]=2*p7*q16 - 2*p7*p16 + 2*p16*q7 - 2*q7*q16;
    h_phc_coeffs_Ht[24]=p7*q10 + p10*q7 - 2*q7*q10;
    h_phc_coeffs_Ht[25]=2*p7*p10 - 2*p7*q10 - 2*p10*q7 + 2*q7*q10;
    h_phc_coeffs_Ht[26]=2*q7*q10 - p10*q7 - p7*q10;
    h_phc_coeffs_Ht[27]=2*p7*q10 - 2*p7*p10 + 2*p10*q7 - 2*q7*q10;
    h_phc_coeffs_Ht[28]=p2*q17 + p17*q2 - 2*q2*q17;
    h_phc_coeffs_Ht[29]=2*p2*p17 - 2*p2*q17 - 2*p17*q2 + 2*q2*q17;
    h_phc_coeffs_Ht[30]=2*q2*q14 - p14*q2 - p2*q14;
    h_phc_coeffs_Ht[31]=2*p2*q14 - 2*p2*p14 + 2*p14*q2 - 2*q2*q14;
    h_phc_coeffs_Ht[32]=p5*q17 + p17*q5 - 2*q5*q17;
    h_phc_coeffs_Ht[33]=2*p5*p17 - 2*p5*q17 - 2*p17*q5 + 2*q5*q17;
    h_phc_coeffs_Ht[34]=2*q5*q14 - p14*q5 - p5*q14;
    h_phc_coeffs_Ht[35]=2*p5*q14 - 2*p5*p14 + 2*p14*q5 - 2*q5*q14;
    h_phc_coeffs_Ht[36]=p8*q17 + p17*q8 - 2*q8*q17;
    h_phc_coeffs_Ht[37]=2*p8*p17 - 2*p8*q17 - 2*p17*q8 + 2*q8*q17;
    h_phc_coeffs_Ht[38]=2*q8*q14 - p14*q8 - p8*q14;
    h_phc_coeffs_Ht[39]=2*p8*q14 - 2*p8*p14 + 2*p14*q8 - 2*q8*q14;
    h_phc_coeffs_Ht[40]=2*q5*q17 - p17*q5 - p5*q17;
    h_phc_coeffs_Ht[41]=2*p5*q17 - 2*p5*p17 + 2*p17*q5 - 2*q5*q17;
    h_phc_coeffs_Ht[42]=2*q2*q17 - p17*q2 - p2*q17;
    h_phc_coeffs_Ht[43]=2*p2*q17 - 2*p2*p17 + 2*p17*q2 - 2*q2*q17;
    h_phc_coeffs_Ht[44]=p8*q14 + p14*q8 - 2*q8*q14;
    h_phc_coeffs_Ht[45]=2*p8*p14 - 2*p8*q14 - 2*p14*q8 + 2*q8*q14;
    h_phc_coeffs_Ht[46]=p2*q11 + p11*q2 - 2*q2*q11;
    h_phc_coeffs_Ht[47]=2*p2*p11 - 2*p2*q11 - 2*p11*q2 + 2*q2*q11;
    h_phc_coeffs_Ht[48]=p5*q11 + p11*q5 - 2*q5*q11;
    h_phc_coeffs_Ht[49]=2*p5*p11 - 2*p5*q11 - 2*p11*q5 + 2*q5*q11;
    h_phc_coeffs_Ht[50]=2*q8*q17 - p17*q8 - p8*q17;
    h_phc_coeffs_Ht[51]=2*p8*q17 - 2*p8*p17 + 2*p17*q8 - 2*q8*q17;
    h_phc_coeffs_Ht[52]=p8*q11 + p11*q8 - 2*q8*q11;
    h_phc_coeffs_Ht[53]=2*p8*p11 - 2*p8*q11 - 2*p11*q8 + 2*q8*q11;
    h_phc_coeffs_Ht[54]=2*q8*q11 - p11*q8 - p8*q11;
    h_phc_coeffs_Ht[55]=2*p8*q11 - 2*p8*p11 + 2*p11*q8 - 2*q8*q11;
    h_phc_coeffs_Ht[56]=p3*q18 + p18*q3 - 2*q3*q18;
    h_phc_coeffs_Ht[57]=2*p3*p18 - 2*p3*q18 - 2*p18*q3 + 2*q3*q18;
    h_phc_coeffs_Ht[58]=2*q3*q15 - p15*q3 - p3*q15;
    h_phc_coeffs_Ht[59]=2*p3*q15 - 2*p3*p15 + 2*p15*q3 - 2*q3*q15;
    h_phc_coeffs_Ht[60]=p6*q18 + p18*q6 - 2*q6*q18;
    h_phc_coeffs_Ht[61]=2*p6*p18 - 2*p6*q18 - 2*p18*q6 + 2*q6*q18;
    h_phc_coeffs_Ht[62]=2*q6*q15 - p15*q6 - p6*q15;
    h_phc_coeffs_Ht[63]=2*p6*q15 - 2*p6*p15 + 2*p15*q6 - 2*q6*q15;
    h_phc_coeffs_Ht[64]=p9*q18 + p18*q9 - 2*q9*q18;
    h_phc_coeffs_Ht[65]=2*p9*p18 - 2*p9*q18 - 2*p18*q9 + 2*q9*q18;
    h_phc_coeffs_Ht[66]=2*q9*q15 - p15*q9 - p9*q15;
    h_phc_coeffs_Ht[67]=2*p9*q15 - 2*p9*p15 + 2*p15*q9 - 2*q9*q15;
    h_phc_coeffs_Ht[68]=2*q6*q18 - p18*q6 - p6*q18;
    h_phc_coeffs_Ht[69]=2*p6*q18 - 2*p6*p18 + 2*p18*q6 - 2*q6*q18;
    h_phc_coeffs_Ht[70]=2*q3*q18 - p18*q3 - p3*q18;
    h_phc_coeffs_Ht[71]=2*p3*q18 - 2*p3*p18 + 2*p18*q3 - 2*q3*q18;
    h_phc_coeffs_Ht[72]=p9*q15 + p15*q9 - 2*q9*q15;
    h_phc_coeffs_Ht[73]=2*p9*p15 - 2*p9*q15 - 2*p15*q9 + 2*q9*q15;
    h_phc_coeffs_Ht[74]=p3*q12 + p12*q3 - 2*q3*q12;
    h_phc_coeffs_Ht[75]=2*p3*p12 - 2*p3*q12 - 2*p12*q3 + 2*q3*q12;
    h_phc_coeffs_Ht[76]=p6*q12 + p12*q6 - 2*q6*q12;
    h_phc_coeffs_Ht[77]=2*p6*p12 - 2*p6*q12 - 2*p12*q6 + 2*q6*q12;
    h_phc_coeffs_Ht[78]=2*q9*q18 - p18*q9 - p9*q18;
    h_phc_coeffs_Ht[79]=2*p9*q18 - 2*p9*p18 + 2*p18*q9 - 2*q9*q18;
    h_phc_coeffs_Ht[80]=p9*q12 + p12*q9 - 2*q9*q12;
    h_phc_coeffs_Ht[81]=2*p9*p12 - 2*p9*q12 - 2*p12*q9 + 2*q9*q12;
    h_phc_coeffs_Ht[82]=2*q9*q12 - p12*q9 - p9*q12;
    h_phc_coeffs_Ht[83]=2*p9*q12 - 2*p9*p12 + 2*p12*q9 - 2*q9*q12;
    h_phc_coeffs_Ht[84]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[85]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[86]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[87]= MAGMA_C_ZERO;
  }
} // end of namespace

#endif
