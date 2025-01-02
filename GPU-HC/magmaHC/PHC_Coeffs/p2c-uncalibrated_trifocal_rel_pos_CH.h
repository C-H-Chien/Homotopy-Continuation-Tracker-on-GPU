#ifndef P2C_UNCALIBRATED_TRIFOCAL_REL_POS_CH_H
#define P2C_UNCALIBRATED_TRIFOCAL_REL_POS_CH_H
// ==========================================================================================================================
//
// Modifications
//    Chiang-Heng Chien  24-05-29:   Initially created. Built on top of generalized three-view cameras minimal problems.
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

  void p2c_uncalibrated_trifocal_rel_pos_CH(
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

    h_phc_coeffs_Hx[0]=q1;
    h_phc_coeffs_Hx[1]=p1 - q1;
    h_phc_coeffs_Hx[2] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[3]=q1 + q17;
    h_phc_coeffs_Hx[4]=p1 + p17 - q1 - q17;
    h_phc_coeffs_Hx[5] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[6]=q1 + q9;
    h_phc_coeffs_Hx[7]=p1 + p9 - q1 - q9;
    h_phc_coeffs_Hx[8] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[9]=q1 - q17;
    h_phc_coeffs_Hx[10]=p1 - p17 - q1 + q17;
    h_phc_coeffs_Hx[11] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[12]=q1 - q9;
    h_phc_coeffs_Hx[13]=p1 - p9 - q1 + q9;
    h_phc_coeffs_Hx[14] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[15]=q1*q13;
    h_phc_coeffs_Hx[16]=p1*q13 + p13*q1 - 2*q1*q13;
    h_phc_coeffs_Hx[17]=p1*p13 - p1*q13 - p13*q1 + q1*q13;
    h_phc_coeffs_Hx[18]=q1*q17;
    h_phc_coeffs_Hx[19]=p1*q17 + p17*q1 - 2*q1*q17;
    h_phc_coeffs_Hx[20]=p1*p17 - p1*q17 - p17*q1 + q1*q17;
    h_phc_coeffs_Hx[21]=q1*q21;
    h_phc_coeffs_Hx[22]=p1*q21 + p21*q1 - 2*q1*q21;
    h_phc_coeffs_Hx[23]=p1*p21 - p1*q21 - p21*q1 + q1*q21;
    h_phc_coeffs_Hx[24]=q1*q9;
    h_phc_coeffs_Hx[25]=p1*q9 + p9*q1 - 2*q1*q9;
    h_phc_coeffs_Hx[26]=p1*p9 - p1*q9 - p9*q1 + q1*q9;
    h_phc_coeffs_Hx[27]=q10;
    h_phc_coeffs_Hx[28]=p10 - q10;
    h_phc_coeffs_Hx[29] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[30]=q11;
    h_phc_coeffs_Hx[31]=p11 - q11;
    h_phc_coeffs_Hx[32] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[33]=q12;
    h_phc_coeffs_Hx[34]=p12 - q12;
    h_phc_coeffs_Hx[35] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[36]=q13;
    h_phc_coeffs_Hx[37]=p13 - q13;
    h_phc_coeffs_Hx[38] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[39]=q14;
    h_phc_coeffs_Hx[40]=p14 - q14;
    h_phc_coeffs_Hx[41] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[42]=q15;
    h_phc_coeffs_Hx[43]=p15 - q15;
    h_phc_coeffs_Hx[44] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[45]=q16;
    h_phc_coeffs_Hx[46]=p16 - q16;
    h_phc_coeffs_Hx[47] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[48]=q17;
    h_phc_coeffs_Hx[49]=p17 - q17;
    h_phc_coeffs_Hx[50] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[51]=q18;
    h_phc_coeffs_Hx[52]=p18 - q18;
    h_phc_coeffs_Hx[53] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[54]=q19;
    h_phc_coeffs_Hx[55]=p19 - q19;
    h_phc_coeffs_Hx[56] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[57]=q2;
    h_phc_coeffs_Hx[58]=p2 - q2;
    h_phc_coeffs_Hx[59] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[60]=q2 + q10;
    h_phc_coeffs_Hx[61]=p2 + p10 - q2 - q10;
    h_phc_coeffs_Hx[62] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[63]=q2 + q18;
    h_phc_coeffs_Hx[64]=p2 + p18 - q2 - q18;
    h_phc_coeffs_Hx[65] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[66]=q2 - q10;
    h_phc_coeffs_Hx[67]=p2 - p10 - q2 + q10;
    h_phc_coeffs_Hx[68] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[69]=q2 - q18;
    h_phc_coeffs_Hx[70]=p2 - p18 - q2 + q18;
    h_phc_coeffs_Hx[71] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[72]=q2*q10;
    h_phc_coeffs_Hx[73]=p2*q10 + p10*q2 - 2*q2*q10;
    h_phc_coeffs_Hx[74]=p2*p10 - p2*q10 - p10*q2 + q2*q10;
    h_phc_coeffs_Hx[75]=q2*q14;
    h_phc_coeffs_Hx[76]=p2*q14 + p14*q2 - 2*q2*q14;
    h_phc_coeffs_Hx[77]=p2*p14 - p2*q14 - p14*q2 + q2*q14;
    h_phc_coeffs_Hx[78]=q2*q18;
    h_phc_coeffs_Hx[79]=p2*q18 + p18*q2 - 2*q2*q18;
    h_phc_coeffs_Hx[80]=p2*p18 - p2*q18 - p18*q2 + q2*q18;
    h_phc_coeffs_Hx[81]=q2*q22;
    h_phc_coeffs_Hx[82]=p2*q22 + p22*q2 - 2*q2*q22;
    h_phc_coeffs_Hx[83]=p2*p22 - p2*q22 - p22*q2 + q2*q22;
    h_phc_coeffs_Hx[84]=q20;
    h_phc_coeffs_Hx[85]=p20 - q20;
    h_phc_coeffs_Hx[86] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[87]=q21;
    h_phc_coeffs_Hx[88]=p21 - q21;
    h_phc_coeffs_Hx[89] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[90]=q22;
    h_phc_coeffs_Hx[91]=p22 - q22;
    h_phc_coeffs_Hx[92] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[93]=q23;
    h_phc_coeffs_Hx[94]=p23 - q23;
    h_phc_coeffs_Hx[95] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[96]=q24;
    h_phc_coeffs_Hx[97]=p24 - q24;
    h_phc_coeffs_Hx[98] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[99]=q3;
    h_phc_coeffs_Hx[100]=p3 - q3;
    h_phc_coeffs_Hx[101] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[102]=q3 + q11;
    h_phc_coeffs_Hx[103]=p3 + p11 - q3 - q11;
    h_phc_coeffs_Hx[104] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[105]=q3 + q19;
    h_phc_coeffs_Hx[106]=p3 + p19 - q3 - q19;
    h_phc_coeffs_Hx[107] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[108]=q3 - q11;
    h_phc_coeffs_Hx[109]=p3 - p11 - q3 + q11;
    h_phc_coeffs_Hx[110] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[111]=q3 - q19;
    h_phc_coeffs_Hx[112]=p3 - p19 - q3 + q19;
    h_phc_coeffs_Hx[113] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[114]=q3*q11;
    h_phc_coeffs_Hx[115]=p3*q11 + p11*q3 - 2*q3*q11;
    h_phc_coeffs_Hx[116]=p3*p11 - p3*q11 - p11*q3 + q3*q11;
    h_phc_coeffs_Hx[117]=q3*q15;
    h_phc_coeffs_Hx[118]=p3*q15 + p15*q3 - 2*q3*q15;
    h_phc_coeffs_Hx[119]=p3*p15 - p3*q15 - p15*q3 + q3*q15;
    h_phc_coeffs_Hx[120]=q3*q19;
    h_phc_coeffs_Hx[121]=p3*q19 + p19*q3 - 2*q3*q19;
    h_phc_coeffs_Hx[122]=p3*p19 - p3*q19 - p19*q3 + q3*q19;
    h_phc_coeffs_Hx[123]=q3*q23;
    h_phc_coeffs_Hx[124]=p3*q23 + p23*q3 - 2*q3*q23;
    h_phc_coeffs_Hx[125]=p3*p23 - p3*q23 - p23*q3 + q3*q23;
    h_phc_coeffs_Hx[126]=q4;
    h_phc_coeffs_Hx[127]=p4 - q4;
    h_phc_coeffs_Hx[128] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[129]=q4 + q12;
    h_phc_coeffs_Hx[130]=p4 + p12 - q4 - q12;
    h_phc_coeffs_Hx[131] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[132]=q4 + q20;
    h_phc_coeffs_Hx[133]=p4 + p20 - q4 - q20;
    h_phc_coeffs_Hx[134] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[135]=q4 - q12;
    h_phc_coeffs_Hx[136]=p4 - p12 - q4 + q12;
    h_phc_coeffs_Hx[137] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[138]=q4 - q20;
    h_phc_coeffs_Hx[139]=p4 - p20 - q4 + q20;
    h_phc_coeffs_Hx[140] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[141]=q4*q12;
    h_phc_coeffs_Hx[142]=p4*q12 + p12*q4 - 2*q4*q12;
    h_phc_coeffs_Hx[143]=p4*p12 - p4*q12 - p12*q4 + q4*q12;
    h_phc_coeffs_Hx[144]=q4*q16;
    h_phc_coeffs_Hx[145]=p4*q16 + p16*q4 - 2*q4*q16;
    h_phc_coeffs_Hx[146]=p4*p16 - p4*q16 - p16*q4 + q4*q16;
    h_phc_coeffs_Hx[147]=q4*q20;
    h_phc_coeffs_Hx[148]=p4*q20 + p20*q4 - 2*q4*q20;
    h_phc_coeffs_Hx[149]=p4*p20 - p4*q20 - p20*q4 + q4*q20;
    h_phc_coeffs_Hx[150]=q4*q24;
    h_phc_coeffs_Hx[151]=p4*q24 + p24*q4 - 2*q4*q24;
    h_phc_coeffs_Hx[152]=p4*p24 - p4*q24 - p24*q4 + q4*q24;
    h_phc_coeffs_Hx[153]=q5;
    h_phc_coeffs_Hx[154]=p5 - q5;
    h_phc_coeffs_Hx[155] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[156]=q5 + q13;
    h_phc_coeffs_Hx[157]=p5 + p13 - q5 - q13;
    h_phc_coeffs_Hx[158] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[159]=q5 + q21;
    h_phc_coeffs_Hx[160]=p5 + p21 - q5 - q21;
    h_phc_coeffs_Hx[161] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[162]=q5 - q13;
    h_phc_coeffs_Hx[163]=p5 - p13 - q5 + q13;
    h_phc_coeffs_Hx[164] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[165]=q5 - q21;
    h_phc_coeffs_Hx[166]=p5 - p21 - q5 + q21;
    h_phc_coeffs_Hx[167] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[168]=q5*q13;
    h_phc_coeffs_Hx[169]=p5*q13 + p13*q5 - 2*q5*q13;
    h_phc_coeffs_Hx[170]=p5*p13 - p5*q13 - p13*q5 + q5*q13;
    h_phc_coeffs_Hx[171]=q5*q17;
    h_phc_coeffs_Hx[172]=p5*q17 + p17*q5 - 2*q5*q17;
    h_phc_coeffs_Hx[173]=p5*p17 - p5*q17 - p17*q5 + q5*q17;
    h_phc_coeffs_Hx[174]=q5*q21;
    h_phc_coeffs_Hx[175]=p5*q21 + p21*q5 - 2*q5*q21;
    h_phc_coeffs_Hx[176]=p5*p21 - p5*q21 - p21*q5 + q5*q21;
    h_phc_coeffs_Hx[177]=q5*q9;
    h_phc_coeffs_Hx[178]=p5*q9 + p9*q5 - 2*q5*q9;
    h_phc_coeffs_Hx[179]=p5*p9 - p5*q9 - p9*q5 + q5*q9;
    h_phc_coeffs_Hx[180]=q6;
    h_phc_coeffs_Hx[181]=p6 - q6;
    h_phc_coeffs_Hx[182] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[183]=q6 + q14;
    h_phc_coeffs_Hx[184]=p6 + p14 - q6 - q14;
    h_phc_coeffs_Hx[185] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[186]=q6 + q22;
    h_phc_coeffs_Hx[187]=p6 + p22 - q6 - q22;
    h_phc_coeffs_Hx[188] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[189]=q6 - q14;
    h_phc_coeffs_Hx[190]=p6 - p14 - q6 + q14;
    h_phc_coeffs_Hx[191] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[192]=q6 - q22;
    h_phc_coeffs_Hx[193]=p6 - p22 - q6 + q22;
    h_phc_coeffs_Hx[194] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[195]=q6*q10;
    h_phc_coeffs_Hx[196]=p6*q10 + p10*q6 - 2*q6*q10;
    h_phc_coeffs_Hx[197]=p6*p10 - p6*q10 - p10*q6 + q6*q10;
    h_phc_coeffs_Hx[198]=q6*q14;
    h_phc_coeffs_Hx[199]=p6*q14 + p14*q6 - 2*q6*q14;
    h_phc_coeffs_Hx[200]=p6*p14 - p6*q14 - p14*q6 + q6*q14;
    h_phc_coeffs_Hx[201]=q6*q18;
    h_phc_coeffs_Hx[202]=p6*q18 + p18*q6 - 2*q6*q18;
    h_phc_coeffs_Hx[203]=p6*p18 - p6*q18 - p18*q6 + q6*q18;
    h_phc_coeffs_Hx[204]=q6*q22;
    h_phc_coeffs_Hx[205]=p6*q22 + p22*q6 - 2*q6*q22;
    h_phc_coeffs_Hx[206]=p6*p22 - p6*q22 - p22*q6 + q6*q22;
    h_phc_coeffs_Hx[207]=q7;
    h_phc_coeffs_Hx[208]=p7 - q7;
    h_phc_coeffs_Hx[209] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[210]=q7 + q15;
    h_phc_coeffs_Hx[211]=p7 + p15 - q7 - q15;
    h_phc_coeffs_Hx[212] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[213]=q7 + q23;
    h_phc_coeffs_Hx[214]=p7 + p23 - q7 - q23;
    h_phc_coeffs_Hx[215] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[216]=q7 - q15;
    h_phc_coeffs_Hx[217]=p7 - p15 - q7 + q15;
    h_phc_coeffs_Hx[218] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[219]=q7 - q23;
    h_phc_coeffs_Hx[220]=p7 - p23 - q7 + q23;
    h_phc_coeffs_Hx[221] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[222]=q7*q11;
    h_phc_coeffs_Hx[223]=p7*q11 + p11*q7 - 2*q7*q11;
    h_phc_coeffs_Hx[224]=p7*p11 - p7*q11 - p11*q7 + q7*q11;
    h_phc_coeffs_Hx[225]=q7*q15;
    h_phc_coeffs_Hx[226]=p7*q15 + p15*q7 - 2*q7*q15;
    h_phc_coeffs_Hx[227]=p7*p15 - p7*q15 - p15*q7 + q7*q15;
    h_phc_coeffs_Hx[228]=q7*q19;
    h_phc_coeffs_Hx[229]=p7*q19 + p19*q7 - 2*q7*q19;
    h_phc_coeffs_Hx[230]=p7*p19 - p7*q19 - p19*q7 + q7*q19;
    h_phc_coeffs_Hx[231]=q7*q23;
    h_phc_coeffs_Hx[232]=p7*q23 + p23*q7 - 2*q7*q23;
    h_phc_coeffs_Hx[233]=p7*p23 - p7*q23 - p23*q7 + q7*q23;
    h_phc_coeffs_Hx[234]=q8;
    h_phc_coeffs_Hx[235]=p8 - q8;
    h_phc_coeffs_Hx[236] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[237]=q8 + q16;
    h_phc_coeffs_Hx[238]=p8 + p16 - q8 - q16;
    h_phc_coeffs_Hx[239] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[240]=q8 + q24;
    h_phc_coeffs_Hx[241]=p8 + p24 - q8 - q24;
    h_phc_coeffs_Hx[242] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[243]=q8 - q16;
    h_phc_coeffs_Hx[244]=p8 - p16 - q8 + q16;
    h_phc_coeffs_Hx[245] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[246]=q8 - q24;
    h_phc_coeffs_Hx[247]=p8 - p24 - q8 + q24;
    h_phc_coeffs_Hx[248] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[249]=q8*q12;
    h_phc_coeffs_Hx[250]=p8*q12 + p12*q8 - 2*q8*q12;
    h_phc_coeffs_Hx[251]=p8*p12 - p8*q12 - p12*q8 + q8*q12;
    h_phc_coeffs_Hx[252]=q8*q16;
    h_phc_coeffs_Hx[253]=p8*q16 + p16*q8 - 2*q8*q16;
    h_phc_coeffs_Hx[254]=p8*p16 - p8*q16 - p16*q8 + q8*q16;
    h_phc_coeffs_Hx[255]=q8*q20;
    h_phc_coeffs_Hx[256]=p8*q20 + p20*q8 - 2*q8*q20;
    h_phc_coeffs_Hx[257]=p8*p20 - p8*q20 - p20*q8 + q8*q20;
    h_phc_coeffs_Hx[258]=q8*q24;
    h_phc_coeffs_Hx[259]=p8*q24 + p24*q8 - 2*q8*q24;
    h_phc_coeffs_Hx[260]=p8*p24 - p8*q24 - p24*q8 + q8*q24;
    h_phc_coeffs_Hx[261]=q9;
    h_phc_coeffs_Hx[262]=p9 - q9;
    h_phc_coeffs_Hx[263] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[264] = MAGMA_C_ONE;
    h_phc_coeffs_Hx[265] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[266] = MAGMA_C_ZERO;

    h_phc_coeffs_Ht[0]=p1 - q1;
    h_phc_coeffs_Ht[1]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[2]=p1 + p17 - q1 - q17;
    h_phc_coeffs_Ht[3]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[4]=p1 + p9 - q1 - q9;
    h_phc_coeffs_Ht[5]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[6]=p1 - p17 - q1 + q17;
    h_phc_coeffs_Ht[7]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[8]=p1 - p9 - q1 + q9;
    h_phc_coeffs_Ht[9]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[10]=p1*q13 + p13*q1 - 2*q1*q13;
    h_phc_coeffs_Ht[11]=2*p1*p13 - 2*p1*q13 - 2*p13*q1 + 2*q1*q13;
    h_phc_coeffs_Ht[12]=p1*q17 + p17*q1 - 2*q1*q17;
    h_phc_coeffs_Ht[13]=2*p1*p17 - 2*p1*q17 - 2*p17*q1 + 2*q1*q17;
    h_phc_coeffs_Ht[14]=p1*q21 + p21*q1 - 2*q1*q21;
    h_phc_coeffs_Ht[15]=2*p1*p21 - 2*p1*q21 - 2*p21*q1 + 2*q1*q21;
    h_phc_coeffs_Ht[16]=p1*q9 + p9*q1 - 2*q1*q9;
    h_phc_coeffs_Ht[17]=2*p1*p9 - 2*p1*q9 - 2*p9*q1 + 2*q1*q9;
    h_phc_coeffs_Ht[18]=p10 - q10;
    h_phc_coeffs_Ht[19]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[20]=p11 - q11;
    h_phc_coeffs_Ht[21]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[22]=p12 - q12;
    h_phc_coeffs_Ht[23]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[24]=p13 - q13;
    h_phc_coeffs_Ht[25]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[26]=p14 - q14;
    h_phc_coeffs_Ht[27]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[28]=p15 - q15;
    h_phc_coeffs_Ht[29]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[30]=p16 - q16;
    h_phc_coeffs_Ht[31]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[32]=p17 - q17;
    h_phc_coeffs_Ht[33]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[34]=p18 - q18;
    h_phc_coeffs_Ht[35]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[36]=p19 - q19;
    h_phc_coeffs_Ht[37]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[38]=p2 - q2;
    h_phc_coeffs_Ht[39]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[40]=p2 + p10 - q2 - q10;
    h_phc_coeffs_Ht[41]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[42]=p2 + p18 - q2 - q18;
    h_phc_coeffs_Ht[43]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[44]=p2 - p10 - q2 + q10;
    h_phc_coeffs_Ht[45]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[46]=p2 - p18 - q2 + q18;
    h_phc_coeffs_Ht[47]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[48]=p2*q10 + p10*q2 - 2*q2*q10;
    h_phc_coeffs_Ht[49]=2*p2*p10 - 2*p2*q10 - 2*p10*q2 + 2*q2*q10;
    h_phc_coeffs_Ht[50]=p2*q14 + p14*q2 - 2*q2*q14;
    h_phc_coeffs_Ht[51]=2*p2*p14 - 2*p2*q14 - 2*p14*q2 + 2*q2*q14;
    h_phc_coeffs_Ht[52]=p2*q18 + p18*q2 - 2*q2*q18;
    h_phc_coeffs_Ht[53]=2*p2*p18 - 2*p2*q18 - 2*p18*q2 + 2*q2*q18;
    h_phc_coeffs_Ht[54]=p2*q22 + p22*q2 - 2*q2*q22;
    h_phc_coeffs_Ht[55]=2*p2*p22 - 2*p2*q22 - 2*p22*q2 + 2*q2*q22;
    h_phc_coeffs_Ht[56]=p20 - q20;
    h_phc_coeffs_Ht[57]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[58]=p21 - q21;
    h_phc_coeffs_Ht[59]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[60]=p22 - q22;
    h_phc_coeffs_Ht[61]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[62]=p23 - q23;
    h_phc_coeffs_Ht[63]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[64]=p24 - q24;
    h_phc_coeffs_Ht[65]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[66]=p3 - q3;
    h_phc_coeffs_Ht[67]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[68]=p3 + p11 - q3 - q11;
    h_phc_coeffs_Ht[69]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[70]=p3 + p19 - q3 - q19;
    h_phc_coeffs_Ht[71]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[72]=p3 - p11 - q3 + q11;
    h_phc_coeffs_Ht[73]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[74]=p3 - p19 - q3 + q19;
    h_phc_coeffs_Ht[75]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[76]=p3*q11 + p11*q3 - 2*q3*q11;
    h_phc_coeffs_Ht[77]=2*p3*p11 - 2*p3*q11 - 2*p11*q3 + 2*q3*q11;
    h_phc_coeffs_Ht[78]=p3*q15 + p15*q3 - 2*q3*q15;
    h_phc_coeffs_Ht[79]=2*p3*p15 - 2*p3*q15 - 2*p15*q3 + 2*q3*q15;
    h_phc_coeffs_Ht[80]=p3*q19 + p19*q3 - 2*q3*q19;
    h_phc_coeffs_Ht[81]=2*p3*p19 - 2*p3*q19 - 2*p19*q3 + 2*q3*q19;
    h_phc_coeffs_Ht[82]=p3*q23 + p23*q3 - 2*q3*q23;
    h_phc_coeffs_Ht[83]=2*p3*p23 - 2*p3*q23 - 2*p23*q3 + 2*q3*q23;
    h_phc_coeffs_Ht[84]=p4 - q4;
    h_phc_coeffs_Ht[85]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[86]=p4 + p12 - q4 - q12;
    h_phc_coeffs_Ht[87]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[88]=p4 + p20 - q4 - q20;
    h_phc_coeffs_Ht[89]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[90]=p4 - p12 - q4 + q12;
    h_phc_coeffs_Ht[91]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[92]=p4 - p20 - q4 + q20;
    h_phc_coeffs_Ht[93]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[94]=p4*q12 + p12*q4 - 2*q4*q12;
    h_phc_coeffs_Ht[95]=2*p4*p12 - 2*p4*q12 - 2*p12*q4 + 2*q4*q12;
    h_phc_coeffs_Ht[96]=p4*q16 + p16*q4 - 2*q4*q16;
    h_phc_coeffs_Ht[97]=2*p4*p16 - 2*p4*q16 - 2*p16*q4 + 2*q4*q16;
    h_phc_coeffs_Ht[98]=p4*q20 + p20*q4 - 2*q4*q20;
    h_phc_coeffs_Ht[99]=2*p4*p20 - 2*p4*q20 - 2*p20*q4 + 2*q4*q20;
    h_phc_coeffs_Ht[100]=p4*q24 + p24*q4 - 2*q4*q24;
    h_phc_coeffs_Ht[101]=2*p4*p24 - 2*p4*q24 - 2*p24*q4 + 2*q4*q24;
    h_phc_coeffs_Ht[102]=p5 - q5;
    h_phc_coeffs_Ht[103]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[104]=p5 + p13 - q5 - q13;
    h_phc_coeffs_Ht[105]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[106]=p5 + p21 - q5 - q21;
    h_phc_coeffs_Ht[107]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[108]=p5 - p13 - q5 + q13;
    h_phc_coeffs_Ht[109]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[110]=p5 - p21 - q5 + q21;
    h_phc_coeffs_Ht[111]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[112]=p5*q13 + p13*q5 - 2*q5*q13;
    h_phc_coeffs_Ht[113]=2*p5*p13 - 2*p5*q13 - 2*p13*q5 + 2*q5*q13;
    h_phc_coeffs_Ht[114]=p5*q17 + p17*q5 - 2*q5*q17;
    h_phc_coeffs_Ht[115]=2*p5*p17 - 2*p5*q17 - 2*p17*q5 + 2*q5*q17;
    h_phc_coeffs_Ht[116]=p5*q21 + p21*q5 - 2*q5*q21;
    h_phc_coeffs_Ht[117]=2*p5*p21 - 2*p5*q21 - 2*p21*q5 + 2*q5*q21;
    h_phc_coeffs_Ht[118]=p5*q9 + p9*q5 - 2*q5*q9;
    h_phc_coeffs_Ht[119]=2*p5*p9 - 2*p5*q9 - 2*p9*q5 + 2*q5*q9;
    h_phc_coeffs_Ht[120]=p6 - q6;
    h_phc_coeffs_Ht[121]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[122]=p6 + p14 - q6 - q14;
    h_phc_coeffs_Ht[123]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[124]=p6 + p22 - q6 - q22;
    h_phc_coeffs_Ht[125]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[126]=p6 - p14 - q6 + q14;
    h_phc_coeffs_Ht[127]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[128]=p6 - p22 - q6 + q22;
    h_phc_coeffs_Ht[129]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[130]=p6*q10 + p10*q6 - 2*q6*q10;
    h_phc_coeffs_Ht[131]=2*p6*p10 - 2*p6*q10 - 2*p10*q6 + 2*q6*q10;
    h_phc_coeffs_Ht[132]=p6*q14 + p14*q6 - 2*q6*q14;
    h_phc_coeffs_Ht[133]=2*p6*p14 - 2*p6*q14 - 2*p14*q6 + 2*q6*q14;
    h_phc_coeffs_Ht[134]=p6*q18 + p18*q6 - 2*q6*q18;
    h_phc_coeffs_Ht[135]=2*p6*p18 - 2*p6*q18 - 2*p18*q6 + 2*q6*q18;
    h_phc_coeffs_Ht[136]=p6*q22 + p22*q6 - 2*q6*q22;
    h_phc_coeffs_Ht[137]=2*p6*p22 - 2*p6*q22 - 2*p22*q6 + 2*q6*q22;
    h_phc_coeffs_Ht[138]=p7 - q7;
    h_phc_coeffs_Ht[139]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[140]=p7 + p15 - q7 - q15;
    h_phc_coeffs_Ht[141]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[142]=p7 + p23 - q7 - q23;
    h_phc_coeffs_Ht[143]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[144]=p7 - p15 - q7 + q15;
    h_phc_coeffs_Ht[145]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[146]=p7 - p23 - q7 + q23;
    h_phc_coeffs_Ht[147]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[148]=p7*q11 + p11*q7 - 2*q7*q11;
    h_phc_coeffs_Ht[149]=2*p7*p11 - 2*p7*q11 - 2*p11*q7 + 2*q7*q11;
    h_phc_coeffs_Ht[150]=p7*q15 + p15*q7 - 2*q7*q15;
    h_phc_coeffs_Ht[151]=2*p7*p15 - 2*p7*q15 - 2*p15*q7 + 2*q7*q15;
    h_phc_coeffs_Ht[152]=p7*q19 + p19*q7 - 2*q7*q19;
    h_phc_coeffs_Ht[153]=2*p7*p19 - 2*p7*q19 - 2*p19*q7 + 2*q7*q19;
    h_phc_coeffs_Ht[154]=p7*q23 + p23*q7 - 2*q7*q23;
    h_phc_coeffs_Ht[155]=2*p7*p23 - 2*p7*q23 - 2*p23*q7 + 2*q7*q23;
    h_phc_coeffs_Ht[156]=p8 - q8;
    h_phc_coeffs_Ht[157]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[158]=p8 + p16 - q8 - q16;
    h_phc_coeffs_Ht[159]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[160]=p8 + p24 - q8 - q24;
    h_phc_coeffs_Ht[161]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[162]=p8 - p16 - q8 + q16;
    h_phc_coeffs_Ht[163]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[164]=p8 - p24 - q8 + q24;
    h_phc_coeffs_Ht[165]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[166]=p8*q12 + p12*q8 - 2*q8*q12;
    h_phc_coeffs_Ht[167]=2*p8*p12 - 2*p8*q12 - 2*p12*q8 + 2*q8*q12;
    h_phc_coeffs_Ht[168]=p8*q16 + p16*q8 - 2*q8*q16;
    h_phc_coeffs_Ht[169]=2*p8*p16 - 2*p8*q16 - 2*p16*q8 + 2*q8*q16;
    h_phc_coeffs_Ht[170]=p8*q20 + p20*q8 - 2*q8*q20;
    h_phc_coeffs_Ht[171]=2*p8*p20 - 2*p8*q20 - 2*p20*q8 + 2*q8*q20;
    h_phc_coeffs_Ht[172]=p8*q24 + p24*q8 - 2*q8*q24;
    h_phc_coeffs_Ht[173]=2*p8*p24 - 2*p8*q24 - 2*p24*q8 + 2*q8*q24;
    h_phc_coeffs_Ht[174]=p9 - q9;
    h_phc_coeffs_Ht[175]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[176]= MAGMA_C_ONE;
    h_phc_coeffs_Ht[177]= MAGMA_C_ZERO;


  }
} // end of namespace

#endif