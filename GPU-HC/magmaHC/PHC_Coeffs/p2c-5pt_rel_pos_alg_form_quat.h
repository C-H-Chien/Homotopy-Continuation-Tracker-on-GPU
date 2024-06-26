#ifndef p2c_5pt_rel_pos_alg_form_quat_h
#define p2c_5pt_rel_pos_alg_form_quat_h

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

namespace magmaHCWrapper{

void p2c_5pt_rel_pos_alg_form_quat(
  magmaFloatComplex *h_targetParams,  magmaFloatComplex *h_startParams, 
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

    h_phc_coeffs_Hx[0]=q1;
    h_phc_coeffs_Hx[1]=p1 - q1;
    h_phc_coeffs_Hx[2] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[3]=q1 - q3;
    h_phc_coeffs_Hx[4]=p1 - p3 - q1 + q3;
    h_phc_coeffs_Hx[5] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[6]=q1*q3 + 1;
    h_phc_coeffs_Hx[7]=p1*q3 + p3*q1 - 2*q1*q3;
    h_phc_coeffs_Hx[8]=(p1 - q1)*(p3 - q3);
    h_phc_coeffs_Hx[9]=q1*q3 + q2*q4;
    h_phc_coeffs_Hx[10]=p1*q3 + p3*q1 + p2*q4 + p4*q2 - 2*q1*q3 - 2*q2*q4;
    h_phc_coeffs_Hx[11]=p1*p3 + p2*p4 - p1*q3 - p3*q1 - p2*q4 - p4*q2 + q1*q3 + q2*q4;
    h_phc_coeffs_Hx[12]=q1*q3 - 1;
    h_phc_coeffs_Hx[13]=p1*q3 + p3*q1 - 2*q1*q3;
    h_phc_coeffs_Hx[14]=(p1 - q1)*(p3 - q3);
    h_phc_coeffs_Hx[15]=q1*q3 - q2*q4;
    h_phc_coeffs_Hx[16]=p1*q3 + p3*q1 - p2*q4 - p4*q2 - 2*q1*q3 + 2*q2*q4;
    h_phc_coeffs_Hx[17]=p1*p3 - p2*p4 - p1*q3 - p3*q1 + p2*q4 + p4*q2 + q1*q3 - q2*q4;
    h_phc_coeffs_Hx[18]=q1*q4;
    h_phc_coeffs_Hx[19]=p1*q4 + p4*q1 - 2*q1*q4;
    h_phc_coeffs_Hx[20]=(p1 - q1)*(p4 - q4);
    h_phc_coeffs_Hx[21]=q1*q4 - q2*q3;
    h_phc_coeffs_Hx[22]=p1*q4 - p2*q3 - p3*q2 + p4*q1 - 2*q1*q4 + 2*q2*q3;
    h_phc_coeffs_Hx[23]=p1*p4 - p2*p3 - p1*q4 + p2*q3 + p3*q2 - p4*q1 + q1*q4 - q2*q3;
    h_phc_coeffs_Hx[24]=q10;
    h_phc_coeffs_Hx[25]=p10 - q10;
    h_phc_coeffs_Hx[26] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[27]=q10 - q12;
    h_phc_coeffs_Hx[28]=p10 - p12 - q10 + q12;
    h_phc_coeffs_Hx[29] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[30]=q10*q11;
    h_phc_coeffs_Hx[31]=p10*q11 + p11*q10 - 2*q10*q11;
    h_phc_coeffs_Hx[32]=(p10 - q10)*(p11 - q11);
    h_phc_coeffs_Hx[33]=q10*q12 + 1;
    h_phc_coeffs_Hx[34]=p10*q12 + p12*q10 - 2*q10*q12;
    h_phc_coeffs_Hx[35]=(p10 - q10)*(p12 - q12);
    h_phc_coeffs_Hx[36]=q10*q12 - 1;
    h_phc_coeffs_Hx[37]=p10*q12 + p12*q10 - 2*q10*q12;
    h_phc_coeffs_Hx[38]=(p10 - q10)*(p12 - q12);
    h_phc_coeffs_Hx[39]=q11;
    h_phc_coeffs_Hx[40]=p11 - q11;
    h_phc_coeffs_Hx[41] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[42]=q12;
    h_phc_coeffs_Hx[43]=p12 - q12;
    h_phc_coeffs_Hx[44] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[45]=q13;
    h_phc_coeffs_Hx[46]=p13 - q13;
    h_phc_coeffs_Hx[47] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[48]=q13 - q15;
    h_phc_coeffs_Hx[49]=p13 - p15 - q13 + q15;
    h_phc_coeffs_Hx[50] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[51]=q13*q15 + 1;
    h_phc_coeffs_Hx[52]=p13*q15 + p15*q13 - 2*q13*q15;
    h_phc_coeffs_Hx[53]=(p13 - q13)*(p15 - q15);
    h_phc_coeffs_Hx[54]=q13*q15 + q14*q16;
    h_phc_coeffs_Hx[55]=p13*q15 + p15*q13 + p14*q16 + p16*q14 - 2*q13*q15 - 2*q14*q16;
    h_phc_coeffs_Hx[56]=p13*p15 + p14*p16 - p13*q15 - p15*q13 - p14*q16 - p16*q14 + q13*q15 + q14*q16;
    h_phc_coeffs_Hx[57]=q13*q15 - 1;
    h_phc_coeffs_Hx[58]=p13*q15 + p15*q13 - 2*q13*q15;
    h_phc_coeffs_Hx[59]=(p13 - q13)*(p15 - q15);
    h_phc_coeffs_Hx[60]=q13*q15 - q14*q16;
    h_phc_coeffs_Hx[61]=p13*q15 + p15*q13 - p14*q16 - p16*q14 - 2*q13*q15 + 2*q14*q16;
    h_phc_coeffs_Hx[62]=p13*p15 - p14*p16 - p13*q15 - p15*q13 + p14*q16 + p16*q14 + q13*q15 - q14*q16;
    h_phc_coeffs_Hx[63]=q13*q16;
    h_phc_coeffs_Hx[64]=p13*q16 + p16*q13 - 2*q13*q16;
    h_phc_coeffs_Hx[65]=(p13 - q13)*(p16 - q16);
    h_phc_coeffs_Hx[66]=q13*q16 - q14*q15;
    h_phc_coeffs_Hx[67]=p13*q16 - p14*q15 - p15*q14 + p16*q13 - 2*q13*q16 + 2*q14*q15;
    h_phc_coeffs_Hx[68]=p13*p16 - p14*p15 - p13*q16 + p14*q15 + p15*q14 - p16*q13 + q13*q16 - q14*q15;
    h_phc_coeffs_Hx[69]=q14;
    h_phc_coeffs_Hx[70]=p14 - q14;
    h_phc_coeffs_Hx[71] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[72]=q14 - q16;
    h_phc_coeffs_Hx[73]=p14 - p16 - q14 + q16;
    h_phc_coeffs_Hx[74] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[75]=q14*q15;
    h_phc_coeffs_Hx[76]=p14*q15 + p15*q14 - 2*q14*q15;
    h_phc_coeffs_Hx[77]=(p14 - q14)*(p15 - q15);
    h_phc_coeffs_Hx[78]=q14*q16 + 1;
    h_phc_coeffs_Hx[79]=p14*q16 + p16*q14 - 2*q14*q16;
    h_phc_coeffs_Hx[80]=(p14 - q14)*(p16 - q16);
    h_phc_coeffs_Hx[81]=q14*q16 - 1;
    h_phc_coeffs_Hx[82]=p14*q16 + p16*q14 - 2*q14*q16;
    h_phc_coeffs_Hx[83]=(p14 - q14)*(p16 - q16);
    h_phc_coeffs_Hx[84]=q15;
    h_phc_coeffs_Hx[85]=p15 - q15;
    h_phc_coeffs_Hx[86] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[87]=q16;
    h_phc_coeffs_Hx[88]=p16 - q16;
    h_phc_coeffs_Hx[89] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[90]=q17;
    h_phc_coeffs_Hx[91]=p17 - q17;
    h_phc_coeffs_Hx[92] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[93]=q17 - q19;
    h_phc_coeffs_Hx[94]=p17 - p19 - q17 + q19;
    h_phc_coeffs_Hx[95] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[96]=q17*q19 + 1;
    h_phc_coeffs_Hx[97]=p17*q19 + p19*q17 - 2*q17*q19;
    h_phc_coeffs_Hx[98]=(p17 - q17)*(p19 - q19);
    h_phc_coeffs_Hx[99]=q17*q19 + q18*q20;
    h_phc_coeffs_Hx[100]=p17*q19 + p19*q17 + p18*q20 + p20*q18 - 2*q17*q19 - 2*q18*q20;
    h_phc_coeffs_Hx[101]=p17*p19 + p18*p20 - p17*q19 - p19*q17 - p18*q20 - p20*q18 + q17*q19 + q18*q20;
    h_phc_coeffs_Hx[102]=q17*q19 - 1;
    h_phc_coeffs_Hx[103]=p17*q19 + p19*q17 - 2*q17*q19;
    h_phc_coeffs_Hx[104]=(p17 - q17)*(p19 - q19);
    h_phc_coeffs_Hx[105]=q17*q19 - q18*q20;
    h_phc_coeffs_Hx[106]=p17*q19 + p19*q17 - p18*q20 - p20*q18 - 2*q17*q19 + 2*q18*q20;
    h_phc_coeffs_Hx[107]=p17*p19 - p18*p20 - p17*q19 - p19*q17 + p18*q20 + p20*q18 + q17*q19 - q18*q20;
    h_phc_coeffs_Hx[108]=q17*q20;
    h_phc_coeffs_Hx[109]=p17*q20 + p20*q17 - 2*q17*q20;
    h_phc_coeffs_Hx[110]=(p17 - q17)*(p20 - q20);
    h_phc_coeffs_Hx[111]=q17*q20 - q18*q19;
    h_phc_coeffs_Hx[112]=p17*q20 - p18*q19 - p19*q18 + p20*q17 - 2*q17*q20 + 2*q18*q19;
    h_phc_coeffs_Hx[113]=p17*p20 - p18*p19 - p17*q20 + p18*q19 + p19*q18 - p20*q17 + q17*q20 - q18*q19;
    h_phc_coeffs_Hx[114]=q18;
    h_phc_coeffs_Hx[115]=p18 - q18;
    h_phc_coeffs_Hx[116] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[117]=q18 - q20;
    h_phc_coeffs_Hx[118]=p18 - p20 - q18 + q20;
    h_phc_coeffs_Hx[119] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[120]=q18*q19;
    h_phc_coeffs_Hx[121]=p18*q19 + p19*q18 - 2*q18*q19;
    h_phc_coeffs_Hx[122]=(p18 - q18)*(p19 - q19);
    h_phc_coeffs_Hx[123]=q18*q20 + 1;
    h_phc_coeffs_Hx[124]=p18*q20 + p20*q18 - 2*q18*q20;
    h_phc_coeffs_Hx[125]=(p18 - q18)*(p20 - q20);
    h_phc_coeffs_Hx[126]=q18*q20 - 1;
    h_phc_coeffs_Hx[127]=p18*q20 + p20*q18 - 2*q18*q20;
    h_phc_coeffs_Hx[128]=(p18 - q18)*(p20 - q20);
    h_phc_coeffs_Hx[129]=q19;
    h_phc_coeffs_Hx[130]=p19 - q19;
    h_phc_coeffs_Hx[131] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[132]=q2;
    h_phc_coeffs_Hx[133]=p2 - q2;
    h_phc_coeffs_Hx[134] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[135]=q2 - q4;
    h_phc_coeffs_Hx[136]=p2 - p4 - q2 + q4;
    h_phc_coeffs_Hx[137] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[138]=q2*q3;
    h_phc_coeffs_Hx[139]=p2*q3 + p3*q2 - 2*q2*q3;
    h_phc_coeffs_Hx[140]=(p2 - q2)*(p3 - q3);
    h_phc_coeffs_Hx[141]=q2*q4 + 1;
    h_phc_coeffs_Hx[142]=p2*q4 + p4*q2 - 2*q2*q4;
    h_phc_coeffs_Hx[143]=(p2 - q2)*(p4 - q4);
    h_phc_coeffs_Hx[144]=q2*q4 - 1;
    h_phc_coeffs_Hx[145]=p2*q4 + p4*q2 - 2*q2*q4;
    h_phc_coeffs_Hx[146]=(p2 - q2)*(p4 - q4);
    h_phc_coeffs_Hx[147]=q20;
    h_phc_coeffs_Hx[148]=p20 - q20;
    h_phc_coeffs_Hx[149] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[150]=q3;
    h_phc_coeffs_Hx[151]=p3 - q3;
    h_phc_coeffs_Hx[152] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[153]=q4;
    h_phc_coeffs_Hx[154]=p4 - q4;
    h_phc_coeffs_Hx[155] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[156]=q5;
    h_phc_coeffs_Hx[157]=p5 - q5;
    h_phc_coeffs_Hx[158] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[159]=q5 - q7;
    h_phc_coeffs_Hx[160]=p5 - p7 - q5 + q7;
    h_phc_coeffs_Hx[161] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[162]=q5*q7 + 1;
    h_phc_coeffs_Hx[163]=p5*q7 + p7*q5 - 2*q5*q7;
    h_phc_coeffs_Hx[164]=(p5 - q5)*(p7 - q7);
    h_phc_coeffs_Hx[165]=q5*q7 + q6*q8;
    h_phc_coeffs_Hx[166]=p5*q7 + p7*q5 + p6*q8 + p8*q6 - 2*q5*q7 - 2*q6*q8;
    h_phc_coeffs_Hx[167]=p5*p7 + p6*p8 - p5*q7 - p7*q5 - p6*q8 - p8*q6 + q5*q7 + q6*q8;
    h_phc_coeffs_Hx[168]=q5*q7 - 1;
    h_phc_coeffs_Hx[169]=p5*q7 + p7*q5 - 2*q5*q7;
    h_phc_coeffs_Hx[170]=(p5 - q5)*(p7 - q7);
    h_phc_coeffs_Hx[171]=q5*q7 - q6*q8;
    h_phc_coeffs_Hx[172]=p5*q7 + p7*q5 - p6*q8 - p8*q6 - 2*q5*q7 + 2*q6*q8;
    h_phc_coeffs_Hx[173]=p5*p7 - p6*p8 - p5*q7 - p7*q5 + p6*q8 + p8*q6 + q5*q7 - q6*q8;
    h_phc_coeffs_Hx[174]=q5*q8;
    h_phc_coeffs_Hx[175]=p5*q8 + p8*q5 - 2*q5*q8;
    h_phc_coeffs_Hx[176]=(p5 - q5)*(p8 - q8);
    h_phc_coeffs_Hx[177]=q5*q8 - q6*q7;
    h_phc_coeffs_Hx[178]=p5*q8 - p6*q7 - p7*q6 + p8*q5 - 2*q5*q8 + 2*q6*q7;
    h_phc_coeffs_Hx[179]=p5*p8 - p6*p7 - p5*q8 + p6*q7 + p7*q6 - p8*q5 + q5*q8 - q6*q7;
    h_phc_coeffs_Hx[180]=q6;
    h_phc_coeffs_Hx[181]=p6 - q6;
    h_phc_coeffs_Hx[182] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[183]=q6 - q8;
    h_phc_coeffs_Hx[184]=p6 - p8 - q6 + q8;
    h_phc_coeffs_Hx[185] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[186]=q6*q7;
    h_phc_coeffs_Hx[187]=p6*q7 + p7*q6 - 2*q6*q7;
    h_phc_coeffs_Hx[188]=(p6 - q6)*(p7 - q7);
    h_phc_coeffs_Hx[189]=q6*q8 + 1;
    h_phc_coeffs_Hx[190]=p6*q8 + p8*q6 - 2*q6*q8;
    h_phc_coeffs_Hx[191]=(p6 - q6)*(p8 - q8);
    h_phc_coeffs_Hx[192]=q6*q8 - 1;
    h_phc_coeffs_Hx[193]=p6*q8 + p8*q6 - 2*q6*q8;
    h_phc_coeffs_Hx[194]=(p6 - q6)*(p8 - q8);
    h_phc_coeffs_Hx[195]=q7;
    h_phc_coeffs_Hx[196]=p7 - q7;
    h_phc_coeffs_Hx[197] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[198]=q8;
    h_phc_coeffs_Hx[199]=p8 - q8;
    h_phc_coeffs_Hx[200] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[201]=q9;
    h_phc_coeffs_Hx[202]=p9 - q9;
    h_phc_coeffs_Hx[203] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[204]=q9 - q11;
    h_phc_coeffs_Hx[205]=p9 - p11 - q9 + q11;
    h_phc_coeffs_Hx[206] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[207]=q9*q11 + 1;
    h_phc_coeffs_Hx[208]=p9*q11 + p11*q9 - 2*q9*q11;
    h_phc_coeffs_Hx[209]=(p9 - q9)*(p11 - q11);
    h_phc_coeffs_Hx[210]=q9*q11 + q10*q12;
    h_phc_coeffs_Hx[211]=p9*q11 + p11*q9 + p10*q12 + p12*q10 - 2*q9*q11 - 2*q10*q12;
    h_phc_coeffs_Hx[212]=p9*p11 + p10*p12 - p9*q11 - p11*q9 - p10*q12 - p12*q10 + q9*q11 + q10*q12;
    h_phc_coeffs_Hx[213]=q9*q11 - 1;
    h_phc_coeffs_Hx[214]=p9*q11 + p11*q9 - 2*q9*q11;
    h_phc_coeffs_Hx[215]=(p9 - q9)*(p11 - q11);
    h_phc_coeffs_Hx[216]=q9*q11 - q10*q12;
    h_phc_coeffs_Hx[217]=p9*q11 + p11*q9 - p10*q12 - p12*q10 - 2*q9*q11 + 2*q10*q12;
    h_phc_coeffs_Hx[218]=p9*p11 - p10*p12 - p9*q11 - p11*q9 + p10*q12 + p12*q10 + q9*q11 - q10*q12;
    h_phc_coeffs_Hx[219]=q9*q12;
    h_phc_coeffs_Hx[220]=p9*q12 + p12*q9 - 2*q9*q12;
    h_phc_coeffs_Hx[221]=(p9 - q9)*(p12 - q12);
    h_phc_coeffs_Hx[222]=q9*q12 - q10*q11;
    h_phc_coeffs_Hx[223]=p9*q12 - p10*q11 - p11*q10 + p12*q9 - 2*q9*q12 + 2*q10*q11;
    h_phc_coeffs_Hx[224]=p9*p12 - p10*p11 - p9*q12 + p10*q11 + p11*q10 - p12*q9 + q9*q12 - q10*q11;
    h_phc_coeffs_Hx[225] = MAGMA_C_ONE;
    h_phc_coeffs_Hx[226] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[227] = MAGMA_C_ZERO;

    h_phc_coeffs_Ht[0]=p1 - q1;
    h_phc_coeffs_Ht[1]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[2]=p1 - p3 - q1 + q3;
    h_phc_coeffs_Ht[3]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[4]=p1*q3 + p3*q1 - 2*q1*q3;
    h_phc_coeffs_Ht[5]=2*(p1 - q1)*(p3 - q3);
    h_phc_coeffs_Ht[6]=p1*q3 + p3*q1 + p2*q4 + p4*q2 - 2*q1*q3 - 2*q2*q4;
    h_phc_coeffs_Ht[7]=2*p1*p3 + 2*p2*p4 - 2*p1*q3 - 2*p3*q1 - 2*p2*q4 - 2*p4*q2 + 2*q1*q3 + 2*q2*q4;
    h_phc_coeffs_Ht[8]=p1*q3 + p3*q1 - 2*q1*q3;
    h_phc_coeffs_Ht[9]=2*(p1 - q1)*(p3 - q3);
    h_phc_coeffs_Ht[10]=p1*q3 + p3*q1 - p2*q4 - p4*q2 - 2*q1*q3 + 2*q2*q4;
    h_phc_coeffs_Ht[11]=2*p1*p3 - 2*p2*p4 - 2*p1*q3 - 2*p3*q1 + 2*p2*q4 + 2*p4*q2 + 2*q1*q3 - 2*q2*q4;
    h_phc_coeffs_Ht[12]=p1*q4 + p4*q1 - 2*q1*q4;
    h_phc_coeffs_Ht[13]=2*(p1 - q1)*(p4 - q4);
    h_phc_coeffs_Ht[14]=p1*q4 - p2*q3 - p3*q2 + p4*q1 - 2*q1*q4 + 2*q2*q3;
    h_phc_coeffs_Ht[15]=2*p1*p4 - 2*p2*p3 - 2*p1*q4 + 2*p2*q3 + 2*p3*q2 - 2*p4*q1 + 2*q1*q4 - 2*q2*q3;
    h_phc_coeffs_Ht[16]=p10 - q10;
    h_phc_coeffs_Ht[17]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[18]=p10 - p12 - q10 + q12;
    h_phc_coeffs_Ht[19]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[20]=p10*q11 + p11*q10 - 2*q10*q11;
    h_phc_coeffs_Ht[21]=2*(p10 - q10)*(p11 - q11);
    h_phc_coeffs_Ht[22]=p10*q12 + p12*q10 - 2*q10*q12;
    h_phc_coeffs_Ht[23]=2*(p10 - q10)*(p12 - q12);
    h_phc_coeffs_Ht[24]=p10*q12 + p12*q10 - 2*q10*q12;
    h_phc_coeffs_Ht[25]=2*(p10 - q10)*(p12 - q12);
    h_phc_coeffs_Ht[26]=p11 - q11;
    h_phc_coeffs_Ht[27]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[28]=p12 - q12;
    h_phc_coeffs_Ht[29]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[30]=p13 - q13;
    h_phc_coeffs_Ht[31]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[32]=p13 - p15 - q13 + q15;
    h_phc_coeffs_Ht[33]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[34]=p13*q15 + p15*q13 - 2*q13*q15;
    h_phc_coeffs_Ht[35]=2*(p13 - q13)*(p15 - q15);
    h_phc_coeffs_Ht[36]=p13*q15 + p15*q13 + p14*q16 + p16*q14 - 2*q13*q15 - 2*q14*q16;
    h_phc_coeffs_Ht[37]=2*p13*p15 + 2*p14*p16 - 2*p13*q15 - 2*p15*q13 - 2*p14*q16 - 2*p16*q14 + 2*q13*q15 + 2*q14*q16;
    h_phc_coeffs_Ht[38]=p13*q15 + p15*q13 - 2*q13*q15;
    h_phc_coeffs_Ht[39]=2*(p13 - q13)*(p15 - q15);
    h_phc_coeffs_Ht[40]=p13*q15 + p15*q13 - p14*q16 - p16*q14 - 2*q13*q15 + 2*q14*q16;
    h_phc_coeffs_Ht[41]=2*p13*p15 - 2*p14*p16 - 2*p13*q15 - 2*p15*q13 + 2*p14*q16 + 2*p16*q14 + 2*q13*q15 - 2*q14*q16;
    h_phc_coeffs_Ht[42]=p13*q16 + p16*q13 - 2*q13*q16;
    h_phc_coeffs_Ht[43]=2*(p13 - q13)*(p16 - q16);
    h_phc_coeffs_Ht[44]=p13*q16 - p14*q15 - p15*q14 + p16*q13 - 2*q13*q16 + 2*q14*q15;
    h_phc_coeffs_Ht[45]=2*p13*p16 - 2*p14*p15 - 2*p13*q16 + 2*p14*q15 + 2*p15*q14 - 2*p16*q13 + 2*q13*q16 - 2*q14*q15;
    h_phc_coeffs_Ht[46]=p14 - q14;
    h_phc_coeffs_Ht[47]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[48]=p14 - p16 - q14 + q16;
    h_phc_coeffs_Ht[49]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[50]=p14*q15 + p15*q14 - 2*q14*q15;
    h_phc_coeffs_Ht[51]=2*(p14 - q14)*(p15 - q15);
    h_phc_coeffs_Ht[52]=p14*q16 + p16*q14 - 2*q14*q16;
    h_phc_coeffs_Ht[53]=2*(p14 - q14)*(p16 - q16);
    h_phc_coeffs_Ht[54]=p14*q16 + p16*q14 - 2*q14*q16;
    h_phc_coeffs_Ht[55]=2*(p14 - q14)*(p16 - q16);
    h_phc_coeffs_Ht[56]=p15 - q15;
    h_phc_coeffs_Ht[57]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[58]=p16 - q16;
    h_phc_coeffs_Ht[59]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[60]=p17 - q17;
    h_phc_coeffs_Ht[61]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[62]=p17 - p19 - q17 + q19;
    h_phc_coeffs_Ht[63]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[64]=p17*q19 + p19*q17 - 2*q17*q19;
    h_phc_coeffs_Ht[65]=2*(p17 - q17)*(p19 - q19);
    h_phc_coeffs_Ht[66]=p17*q19 + p19*q17 + p18*q20 + p20*q18 - 2*q17*q19 - 2*q18*q20;
    h_phc_coeffs_Ht[67]=2*p17*p19 + 2*p18*p20 - 2*p17*q19 - 2*p19*q17 - 2*p18*q20 - 2*p20*q18 + 2*q17*q19 + 2*q18*q20;
    h_phc_coeffs_Ht[68]=p17*q19 + p19*q17 - 2*q17*q19;
    h_phc_coeffs_Ht[69]=2*(p17 - q17)*(p19 - q19);
    h_phc_coeffs_Ht[70]=p17*q19 + p19*q17 - p18*q20 - p20*q18 - 2*q17*q19 + 2*q18*q20;
    h_phc_coeffs_Ht[71]=2*p17*p19 - 2*p18*p20 - 2*p17*q19 - 2*p19*q17 + 2*p18*q20 + 2*p20*q18 + 2*q17*q19 - 2*q18*q20;
    h_phc_coeffs_Ht[72]=p17*q20 + p20*q17 - 2*q17*q20;
    h_phc_coeffs_Ht[73]=2*(p17 - q17)*(p20 - q20);
    h_phc_coeffs_Ht[74]=p17*q20 - p18*q19 - p19*q18 + p20*q17 - 2*q17*q20 + 2*q18*q19;
    h_phc_coeffs_Ht[75]=2*p17*p20 - 2*p18*p19 - 2*p17*q20 + 2*p18*q19 + 2*p19*q18 - 2*p20*q17 + 2*q17*q20 - 2*q18*q19;
    h_phc_coeffs_Ht[76]=p18 - q18;
    h_phc_coeffs_Ht[77]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[78]=p18 - p20 - q18 + q20;
    h_phc_coeffs_Ht[79]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[80]=p18*q19 + p19*q18 - 2*q18*q19;
    h_phc_coeffs_Ht[81]=2*(p18 - q18)*(p19 - q19);
    h_phc_coeffs_Ht[82]=p18*q20 + p20*q18 - 2*q18*q20;
    h_phc_coeffs_Ht[83]=2*(p18 - q18)*(p20 - q20);
    h_phc_coeffs_Ht[84]=p18*q20 + p20*q18 - 2*q18*q20;
    h_phc_coeffs_Ht[85]=2*(p18 - q18)*(p20 - q20);
    h_phc_coeffs_Ht[86]=p19 - q19;
    h_phc_coeffs_Ht[87]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[88]=p2 - q2;
    h_phc_coeffs_Ht[89]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[90]=p2 - p4 - q2 + q4;
    h_phc_coeffs_Ht[91]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[92]=p2*q3 + p3*q2 - 2*q2*q3;
    h_phc_coeffs_Ht[93]=2*(p2 - q2)*(p3 - q3);
    h_phc_coeffs_Ht[94]=p2*q4 + p4*q2 - 2*q2*q4;
    h_phc_coeffs_Ht[95]=2*(p2 - q2)*(p4 - q4);
    h_phc_coeffs_Ht[96]=p2*q4 + p4*q2 - 2*q2*q4;
    h_phc_coeffs_Ht[97]=2*(p2 - q2)*(p4 - q4);
    h_phc_coeffs_Ht[98]=p20 - q20;
    h_phc_coeffs_Ht[99]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[100]=p3 - q3;
    h_phc_coeffs_Ht[101]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[102]=p4 - q4;
    h_phc_coeffs_Ht[103]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[104]=p5 - q5;
    h_phc_coeffs_Ht[105]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[106]=p5 - p7 - q5 + q7;
    h_phc_coeffs_Ht[107]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[108]=p5*q7 + p7*q5 - 2*q5*q7;
    h_phc_coeffs_Ht[109]=2*(p5 - q5)*(p7 - q7);
    h_phc_coeffs_Ht[110]=p5*q7 + p7*q5 + p6*q8 + p8*q6 - 2*q5*q7 - 2*q6*q8;
    h_phc_coeffs_Ht[111]=2*p5*p7 + 2*p6*p8 - 2*p5*q7 - 2*p7*q5 - 2*p6*q8 - 2*p8*q6 + 2*q5*q7 + 2*q6*q8;
    h_phc_coeffs_Ht[112]=p5*q7 + p7*q5 - 2*q5*q7;
    h_phc_coeffs_Ht[113]=2*(p5 - q5)*(p7 - q7);
    h_phc_coeffs_Ht[114]=p5*q7 + p7*q5 - p6*q8 - p8*q6 - 2*q5*q7 + 2*q6*q8;
    h_phc_coeffs_Ht[115]=2*p5*p7 - 2*p6*p8 - 2*p5*q7 - 2*p7*q5 + 2*p6*q8 + 2*p8*q6 + 2*q5*q7 - 2*q6*q8;
    h_phc_coeffs_Ht[116]=p5*q8 + p8*q5 - 2*q5*q8;
    h_phc_coeffs_Ht[117]=2*(p5 - q5)*(p8 - q8);
    h_phc_coeffs_Ht[118]=p5*q8 - p6*q7 - p7*q6 + p8*q5 - 2*q5*q8 + 2*q6*q7;
    h_phc_coeffs_Ht[119]=2*p5*p8 - 2*p6*p7 - 2*p5*q8 + 2*p6*q7 + 2*p7*q6 - 2*p8*q5 + 2*q5*q8 - 2*q6*q7;
    h_phc_coeffs_Ht[120]=p6 - q6;
    h_phc_coeffs_Ht[121]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[122]=p6 - p8 - q6 + q8;
    h_phc_coeffs_Ht[123]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[124]=p6*q7 + p7*q6 - 2*q6*q7;
    h_phc_coeffs_Ht[125]=2*(p6 - q6)*(p7 - q7);
    h_phc_coeffs_Ht[126]=p6*q8 + p8*q6 - 2*q6*q8;
    h_phc_coeffs_Ht[127]=2*(p6 - q6)*(p8 - q8);
    h_phc_coeffs_Ht[128]=p6*q8 + p8*q6 - 2*q6*q8;
    h_phc_coeffs_Ht[129]=2*(p6 - q6)*(p8 - q8);
    h_phc_coeffs_Ht[130]=p7 - q7;
    h_phc_coeffs_Ht[131]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[132]=p8 - q8;
    h_phc_coeffs_Ht[133]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[134]=p9 - q9;
    h_phc_coeffs_Ht[135]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[136]=p9 - p11 - q9 + q11;
    h_phc_coeffs_Ht[137]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[138]=p9*q11 + p11*q9 - 2*q9*q11;
    h_phc_coeffs_Ht[139]=2*(p9 - q9)*(p11 - q11);
    h_phc_coeffs_Ht[140]=p9*q11 + p11*q9 + p10*q12 + p12*q10 - 2*q9*q11 - 2*q10*q12;
    h_phc_coeffs_Ht[141]=2*p9*p11 + 2*p10*p12 - 2*p9*q11 - 2*p11*q9 - 2*p10*q12 - 2*p12*q10 + 2*q9*q11 + 2*q10*q12;
    h_phc_coeffs_Ht[142]=p9*q11 + p11*q9 - 2*q9*q11;
    h_phc_coeffs_Ht[143]=2*(p9 - q9)*(p11 - q11);
    h_phc_coeffs_Ht[144]=p9*q11 + p11*q9 - p10*q12 - p12*q10 - 2*q9*q11 + 2*q10*q12;
    h_phc_coeffs_Ht[145]=2*p9*p11 - 2*p10*p12 - 2*p9*q11 - 2*p11*q9 + 2*p10*q12 + 2*p12*q10 + 2*q9*q11 - 2*q10*q12;
    h_phc_coeffs_Ht[146]=p9*q12 + p12*q9 - 2*q9*q12;
    h_phc_coeffs_Ht[147]=2*(p9 - q9)*(p12 - q12);
    h_phc_coeffs_Ht[148]=p9*q12 - p10*q11 - p11*q10 + p12*q9 - 2*q9*q12 + 2*q10*q11;
    h_phc_coeffs_Ht[149]=2*p9*p12 - 2*p10*p11 - 2*p9*q12 + 2*p10*q11 + 2*p11*q10 - 2*p12*q9 + 2*q9*q12 - 2*q10*q11;
    h_phc_coeffs_Ht[150] = MAGMA_C_ONE;
    h_phc_coeffs_Ht[151] = MAGMA_C_ZERO;
  }
} // end of namespace
#endif