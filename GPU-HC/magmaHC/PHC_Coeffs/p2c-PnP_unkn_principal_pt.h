#ifndef P2C_PNP_UNKN_PRINCIPAL_PT_H
#define P2C_PNP_UNKN_PRINCIPAL_PT_H
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

  void p2c_PnP_unkn_principal_pt(
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
    magmaFloatComplex p34 = h_targetParams[33];
    
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
    magmaFloatComplex q34 = h_startParams[33];

    h_phc_coeffs_Hx[0]=q1 - q10;
    h_phc_coeffs_Hx[1]=p1 - p10 - q1 + q10;
    h_phc_coeffs_Hx[2] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[3]=q1 - q13;
    h_phc_coeffs_Hx[4]=p1 - p13 - q1 + q13;
    h_phc_coeffs_Hx[5] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[6]=q1 - q4;
    h_phc_coeffs_Hx[7]=p1 - p4 - q1 + q4;
    h_phc_coeffs_Hx[8] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[9]=q1 - q7;
    h_phc_coeffs_Hx[10]=p1 - p7 - q1 + q7;
    h_phc_coeffs_Hx[11] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[12]=q16 - q19;
    h_phc_coeffs_Hx[13]=p16 - p19 - q16 + q19;
    h_phc_coeffs_Hx[14] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[15]=q16 - q22;
    h_phc_coeffs_Hx[16]=p16 - p22 - q16 + q22;
    h_phc_coeffs_Hx[17] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[18]=q16 - q25;
    h_phc_coeffs_Hx[19]=p16 - p25 - q16 + q25;
    h_phc_coeffs_Hx[20] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[21]=q16 - q28;
    h_phc_coeffs_Hx[22]=p16 - p28 - q16 + q28;
    h_phc_coeffs_Hx[23] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[24]=q17 - q20;
    h_phc_coeffs_Hx[25]=p17 - p20 - q17 + q20;
    h_phc_coeffs_Hx[26] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[27]=q17 - q23;
    h_phc_coeffs_Hx[28]=p17 - p23 - q17 + q23;
    h_phc_coeffs_Hx[29] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[30]=q17 - q26;
    h_phc_coeffs_Hx[31]=p17 - p26 - q17 + q26;
    h_phc_coeffs_Hx[32] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[33]=q17 - q29;
    h_phc_coeffs_Hx[34]=p17 - p29 - q17 + q29;
    h_phc_coeffs_Hx[35] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[36]=q18 - q21 + q16*q33 + q17*q34 - q19*q33 - q20*q34;
    h_phc_coeffs_Hx[37]=p18 - p21 - q18 + q21 + p16*q33 + p33*q16 + p17*q34 + p34*q17 - p19*q33 - p33*q19 - p20*q34 - p34*q20 - 2*q16*q33 - 2*q17*q34 + 2*q19*q33 + 2*q20*q34;
    h_phc_coeffs_Hx[38]=p16*p33 + p17*p34 - p19*p33 - p20*p34 - p16*q33 - p33*q16 - p17*q34 - p34*q17 + p19*q33 + p33*q19 + p20*q34 + p34*q20 + q16*q33 + q17*q34 - q19*q33 - q20*q34;
    h_phc_coeffs_Hx[39]=q18 - q24 + q16*q33 + q17*q34 - q22*q33 - q23*q34;
    h_phc_coeffs_Hx[40]=p18 - p24 - q18 + q24 + p16*q33 + p33*q16 + p17*q34 + p34*q17 - p22*q33 - p33*q22 - p23*q34 - p34*q23 - 2*q16*q33 - 2*q17*q34 + 2*q22*q33 + 2*q23*q34;
    h_phc_coeffs_Hx[41]=p16*p33 + p17*p34 - p22*p33 - p23*p34 - p16*q33 - p33*q16 - p17*q34 - p34*q17 + p22*q33 + p33*q22 + p23*q34 + p34*q23 + q16*q33 + q17*q34 - q22*q33 - q23*q34;
    h_phc_coeffs_Hx[42]=q18 - q27 + q16*q33 + q17*q34 - q25*q33 - q26*q34;
    h_phc_coeffs_Hx[43]=p18 - p27 - q18 + q27 + p16*q33 + p33*q16 + p17*q34 + p34*q17 - p25*q33 - p33*q25 - p26*q34 - p34*q26 - 2*q16*q33 - 2*q17*q34 + 2*q25*q33 + 2*q26*q34;
    h_phc_coeffs_Hx[44]=p16*p33 + p17*p34 - p25*p33 - p26*p34 - p16*q33 - p33*q16 - p17*q34 - p34*q17 + p25*q33 + p33*q25 + p26*q34 + p34*q26 + q16*q33 + q17*q34 - q25*q33 - q26*q34;
    h_phc_coeffs_Hx[45]=q18 - q30 + q16*q33 + q17*q34 - q28*q33 - q29*q34;
    h_phc_coeffs_Hx[46]=p18 - p30 - q18 + q30 + p16*q33 + p33*q16 + p17*q34 + p34*q17 - p28*q33 - p33*q28 - p29*q34 - p34*q29 - 2*q16*q33 - 2*q17*q34 + 2*q28*q33 + 2*q29*q34;
    h_phc_coeffs_Hx[47]=p16*p33 + p17*p34 - p28*p33 - p29*p34 - p16*q33 - p33*q16 - p17*q34 - p34*q17 + p28*q33 + p33*q28 + p29*q34 + p34*q29 + q16*q33 + q17*q34 - q28*q33 - q29*q34;
    h_phc_coeffs_Hx[48]=q1*q1+q2*q2+q3*q3-q10*q10-q11*q11-q12*q12-q16*q16-q17*q17-q18*q18+q25*q25+q26*q26+q27*q27;
    h_phc_coeffs_Hx[49]=2*p1*q1 + 2*p2*q2 + 2*p3*q3 - 2*p10*q10 - 2*p11*q11 - 2*p12*q12 - 2*p16*q16 - 2*p17*q17 - 2*p18*q18 + 2*p25*q25 + 2*p26*q26 + 2*p27*q27 - 2*q1*q1- 2*q2*q2- 2*q3*q3+ 2*q10*q10+ 2*q11*q11+ 2*q12*q12+ 2*q16*q16+ 2*q17*q17+ 2*q18*q18- 2*q25*q25- 2*q26*q26- 2*q27*q27;
    h_phc_coeffs_Hx[50]=2*p10*q10 - 2*p2*q2 - 2*p3*q3 - 2*p1*q1 + 2*p11*q11 + 2*p12*q12 + 2*p16*q16 + 2*p17*q17 + 2*p18*q18 - 2*p25*q25 - 2*p26*q26 - 2*p27*q27 +p1*p1+p2*p2+p3*p3-p10*p10-p11*p11-p12*p12-p16*p16-p17*p17-p18*p18+p25*p25+p26*p26+p27*p27+q1*q1+q2*q2+q3*q3-q10*q10-q11*q11-q12*q12-q16*q16-q17*q17-q18*q18+q25*q25+q26*q26+q27*q27;
    h_phc_coeffs_Hx[51]=q1*q1+q2*q2+q3*q3-q13*q13-q14*q14-q15*q15-q16*q16-q17*q17-q18*q18+q28*q28+q29*q29+q30*q30;
    h_phc_coeffs_Hx[52]=2*p1*q1 + 2*p2*q2 + 2*p3*q3 - 2*p13*q13 - 2*p14*q14 - 2*p15*q15 - 2*p16*q16 - 2*p17*q17 - 2*p18*q18 + 2*p28*q28 + 2*p29*q29 + 2*p30*q30 - 2*q1*q1- 2*q2*q2- 2*q3*q3+ 2*q13*q13+ 2*q14*q14+ 2*q15*q15+ 2*q16*q16+ 2*q17*q17+ 2*q18*q18- 2*q28*q28- 2*q29*q29- 2*q30*q30;
    h_phc_coeffs_Hx[53]=2*p13*q13 - 2*p2*q2 - 2*p3*q3 - 2*p1*q1 + 2*p14*q14 + 2*p15*q15 + 2*p16*q16 + 2*p17*q17 + 2*p18*q18 - 2*p28*q28 - 2*p29*q29 - 2*p30*q30 +p1*p1+p2*p2+p3*p3-p13*p13-p14*p14-p15*p15-p16*p16-p17*p17-p18*p18+p28*p28+p29*p29+p30*p30+q1*q1+q2*q2+q3*q3-q13*q13-q14*q14-q15*q15-q16*q16-q17*q17-q18*q18+q28*q28+q29*q29+q30*q30;
    h_phc_coeffs_Hx[54]=q1*q1+q2*q2+q3*q3-q4*q4-q5*q5-q6*q6-q16*q16-q17*q17-q18*q18+q19*q19+q20*q20+q21*q21;
    h_phc_coeffs_Hx[55]=2*p1*q1 + 2*p2*q2 + 2*p3*q3 - 2*p4*q4 - 2*p5*q5 - 2*p6*q6 - 2*p16*q16 - 2*p17*q17 - 2*p18*q18 + 2*p19*q19 + 2*p20*q20 + 2*p21*q21 - 2*q1*q1- 2*q2*q2- 2*q3*q3+ 2*q4*q4+ 2*q5*q5+ 2*q6*q6+ 2*q16*q16+ 2*q17*q17+ 2*q18*q18- 2*q19*q19- 2*q20*q20- 2*q21*q21;
    h_phc_coeffs_Hx[56]=2*p4*q4 - 2*p2*q2 - 2*p3*q3 - 2*p1*q1 + 2*p5*q5 + 2*p6*q6 + 2*p16*q16 + 2*p17*q17 + 2*p18*q18 - 2*p19*q19 - 2*p20*q20 - 2*p21*q21 +p1*p1+p2*p2+p3*p3-p4*p4-p5*p5-p6*p6-p16*p16-p17*p17-p18*p18+p19*p19+p20*p20+p21*p21+q1*q1+q2*q2+q3*q3-q4*q4-q5*q5-q6*q6-q16*q16-q17*q17-q18*q18+q19*q19+q20*q20+q21*q21;
    h_phc_coeffs_Hx[57]=q1*q1+q2*q2+q3*q3-q7*q7-q8*q8-q9*q9-q16*q16-q17*q17-q18*q18+q22*q22+q23*q23+q24*q24;
    h_phc_coeffs_Hx[58]=2*p1*q1 + 2*p2*q2 + 2*p3*q3 - 2*p7*q7 - 2*p8*q8 - 2*p9*q9 - 2*p16*q16 - 2*p17*q17 - 2*p18*q18 + 2*p22*q22 + 2*p23*q23 + 2*p24*q24 - 2*q1*q1- 2*q2*q2- 2*q3*q3+ 2*q7*q7+ 2*q8*q8+ 2*q9*q9+ 2*q16*q16+ 2*q17*q17+ 2*q18*q18- 2*q22*q22- 2*q23*q23- 2*q24*q24;
    h_phc_coeffs_Hx[59]=2*p7*q7 - 2*p2*q2 - 2*p3*q3 - 2*p1*q1 + 2*p8*q8 + 2*p9*q9 + 2*p16*q16 + 2*p17*q17 + 2*p18*q18 - 2*p22*q22 - 2*p23*q23 - 2*p24*q24 +p1*p1+p2*p2+p3*p3-p7*p7-p8*p8-p9*p9-p16*p16-p17*p17-p18*p18+p22*p22+p23*p23+p24*p24+q1*q1+q2*q2+q3*q3-q7*q7-q8*q8-q9*q9-q16*q16-q17*q17-q18*q18+q22*q22+q23*q23+q24*q24;
    h_phc_coeffs_Hx[60]=q2 - q11;
    h_phc_coeffs_Hx[61]=p2 - p11 - q2 + q11;
    h_phc_coeffs_Hx[62] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[63]=q2 - q14;
    h_phc_coeffs_Hx[64]=p2 - p14 - q2 + q14;
    h_phc_coeffs_Hx[65] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[66]=q2 - q5;
    h_phc_coeffs_Hx[67]=p2 - p5 - q2 + q5;
    h_phc_coeffs_Hx[68] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[69]=q2 - q8;
    h_phc_coeffs_Hx[70]=p2 - p8 - q2 + q8;
    h_phc_coeffs_Hx[71] = MAGMA_C_ZERO;
    h_phc_coeffs_Hx[72]=q3 - q12 + q1*q31 + q2*q32 - q10*q31 - q11*q32;
    h_phc_coeffs_Hx[73]=p3 - p12 - q3 + q12 + p1*q31 + p31*q1 + p2*q32 + p32*q2 - p10*q31 - p31*q10 - p11*q32 - p32*q11 - 2*q1*q31 - 2*q2*q32 + 2*q10*q31 + 2*q11*q32;
    h_phc_coeffs_Hx[74]=p1*p31 + p2*p32 - p10*p31 - p11*p32 - p1*q31 - p31*q1 - p2*q32 - p32*q2 + p10*q31 + p31*q10 + p11*q32 + p32*q11 + q1*q31 + q2*q32 - q10*q31 - q11*q32;
    h_phc_coeffs_Hx[75]=q3 - q15 + q1*q31 + q2*q32 - q13*q31 - q14*q32;
    h_phc_coeffs_Hx[76]=p3 - p15 - q3 + q15 + p1*q31 + p31*q1 + p2*q32 + p32*q2 - p13*q31 - p31*q13 - p14*q32 - p32*q14 - 2*q1*q31 - 2*q2*q32 + 2*q13*q31 + 2*q14*q32;
    h_phc_coeffs_Hx[77]=p1*p31 + p2*p32 - p13*p31 - p14*p32 - p1*q31 - p31*q1 - p2*q32 - p32*q2 + p13*q31 + p31*q13 + p14*q32 + p32*q14 + q1*q31 + q2*q32 - q13*q31 - q14*q32;
    h_phc_coeffs_Hx[78]=q3 - q6 + q1*q31 + q2*q32 - q4*q31 - q5*q32;
    h_phc_coeffs_Hx[79]=p3 - p6 - q3 + q6 + p1*q31 + p31*q1 + p2*q32 + p32*q2 - p4*q31 - p31*q4 - p5*q32 - p32*q5 - 2*q1*q31 - 2*q2*q32 + 2*q4*q31 + 2*q5*q32;
    h_phc_coeffs_Hx[80]=p1*p31 + p2*p32 - p4*p31 - p5*p32 - p1*q31 - p31*q1 - p2*q32 - p32*q2 + p4*q31 + p31*q4 + p5*q32 + p32*q5 + q1*q31 + q2*q32 - q4*q31 - q5*q32;
    h_phc_coeffs_Hx[81]=q3 - q9 + q1*q31 + q2*q32 - q7*q31 - q8*q32;
    h_phc_coeffs_Hx[82]=p3 - p9 - q3 + q9 + p1*q31 + p31*q1 + p2*q32 + p32*q2 - p7*q31 - p31*q7 - p8*q32 - p32*q8 - 2*q1*q31 - 2*q2*q32 + 2*q7*q31 + 2*q8*q32;
    h_phc_coeffs_Hx[83]=p1*p31 + p2*p32 - p7*p31 - p8*p32 - p1*q31 - p31*q1 - p2*q32 - p32*q2 + p7*q31 + p31*q7 + p8*q32 + p32*q8 + q1*q31 + q2*q32 - q7*q31 - q8*q32;
    h_phc_coeffs_Hx[84]=MAGMA_C_ONE;
    h_phc_coeffs_Hx[85]=MAGMA_C_ZERO;
    h_phc_coeffs_Hx[86]=MAGMA_C_ZERO;

    h_phc_coeffs_Ht[0]=p1 - p10 - q1 + q10;
    h_phc_coeffs_Ht[1]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[2]=p1 - p13 - q1 + q13;
    h_phc_coeffs_Ht[3]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[4]=p1 - p4 - q1 + q4;
    h_phc_coeffs_Ht[5]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[6]=p1 - p7 - q1 + q7;
    h_phc_coeffs_Ht[7]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[8]=p16 - p19 - q16 + q19;
    h_phc_coeffs_Ht[9]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[10]=p16 - p22 - q16 + q22;
    h_phc_coeffs_Ht[11]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[12]=p16 - p25 - q16 + q25;
    h_phc_coeffs_Ht[13]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[14]=p16 - p28 - q16 + q28;
    h_phc_coeffs_Ht[15]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[16]=p17 - p20 - q17 + q20;
    h_phc_coeffs_Ht[17]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[18]=p17 - p23 - q17 + q23;
    h_phc_coeffs_Ht[19]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[20]=p17 - p26 - q17 + q26;
    h_phc_coeffs_Ht[21]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[22]=p17 - p29 - q17 + q29;
    h_phc_coeffs_Ht[23]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[24]=p18 - p21 - q18 + q21 + p16*q33 + p33*q16 + p17*q34 + p34*q17 - p19*q33 - p33*q19 - p20*q34 - p34*q20 - 2*q16*q33 - 2*q17*q34 + 2*q19*q33 + 2*q20*q34;
    h_phc_coeffs_Ht[25]=2*p16*p33 + 2*p17*p34 - 2*p19*p33 - 2*p20*p34 - 2*p16*q33 - 2*p33*q16 - 2*p17*q34 - 2*p34*q17 + 2*p19*q33 + 2*p33*q19 + 2*p20*q34 + 2*p34*q20 + 2*q16*q33 + 2*q17*q34 - 2*q19*q33 - 2*q20*q34;
    h_phc_coeffs_Ht[26]=p18 - p24 - q18 + q24 + p16*q33 + p33*q16 + p17*q34 + p34*q17 - p22*q33 - p33*q22 - p23*q34 - p34*q23 - 2*q16*q33 - 2*q17*q34 + 2*q22*q33 + 2*q23*q34;
    h_phc_coeffs_Ht[27]=2*p16*p33 + 2*p17*p34 - 2*p22*p33 - 2*p23*p34 - 2*p16*q33 - 2*p33*q16 - 2*p17*q34 - 2*p34*q17 + 2*p22*q33 + 2*p33*q22 + 2*p23*q34 + 2*p34*q23 + 2*q16*q33 + 2*q17*q34 - 2*q22*q33 - 2*q23*q34;
    h_phc_coeffs_Ht[28]=p18 - p27 - q18 + q27 + p16*q33 + p33*q16 + p17*q34 + p34*q17 - p25*q33 - p33*q25 - p26*q34 - p34*q26 - 2*q16*q33 - 2*q17*q34 + 2*q25*q33 + 2*q26*q34;
    h_phc_coeffs_Ht[29]=2*p16*p33 + 2*p17*p34 - 2*p25*p33 - 2*p26*p34 - 2*p16*q33 - 2*p33*q16 - 2*p17*q34 - 2*p34*q17 + 2*p25*q33 + 2*p33*q25 + 2*p26*q34 + 2*p34*q26 + 2*q16*q33 + 2*q17*q34 - 2*q25*q33 - 2*q26*q34;
    h_phc_coeffs_Ht[30]=p18 - p30 - q18 + q30 + p16*q33 + p33*q16 + p17*q34 + p34*q17 - p28*q33 - p33*q28 - p29*q34 - p34*q29 - 2*q16*q33 - 2*q17*q34 + 2*q28*q33 + 2*q29*q34;
    h_phc_coeffs_Ht[31]=2*p16*p33 + 2*p17*p34 - 2*p28*p33 - 2*p29*p34 - 2*p16*q33 - 2*p33*q16 - 2*p17*q34 - 2*p34*q17 + 2*p28*q33 + 2*p33*q28 + 2*p29*q34 + 2*p34*q29 + 2*q16*q33 + 2*q17*q34 - 2*q28*q33 - 2*q29*q34;
    h_phc_coeffs_Ht[32]=2*q10*q10- 2*q2*q2- 2*q3*q3- 2*q1*q1+ 2*q11*q11+ 2*q12*q12+ 2*q16*q16+ 2*q17*q17+ 2*q18*q18- 2*q25*q25- 2*q26*q26- 2*q27*q27+ 2*p1*q1 + 2*p2*q2 + 2*p3*q3 - 2*p10*q10 - 2*p11*q11 - 2*p12*q12 - 2*p16*q16 - 2*p17*q17 - 2*p18*q18 + 2*p25*q25 + 2*p26*q26 + 2*p27*q27;
    h_phc_coeffs_Ht[33]=2*p1*p1- 4*p1*q1 + 2*p2*p2- 4*p2*q2 + 2*p3*p3- 4*p3*q3 - 2*p10*p10+ 4*p10*q10 - 2*p11*p11+ 4*p11*q11 - 2*p12*p12+ 4*p12*q12 - 2*p16*p16+ 4*p16*q16 - 2*p17*p17+ 4*p17*q17 - 2*p18*p18+ 4*p18*q18 + 2*p25*p25- 4*p25*q25 + 2*p26*p26- 4*p26*q26 + 2*p27*p27- 4*p27*q27 + 2*q1*q1+ 2*q2*q2+ 2*q3*q3- 2*q10*q10- 2*q11*q11- 2*q12*q12- 2*q16*q16- 2*q17*q17- 2*q18*q18+ 2*q25*q25+ 2*q26*q26+ 2*q27*q27;
    h_phc_coeffs_Ht[34]=2*q13*q13- 2*q2*q2- 2*q3*q3- 2*q1*q1+ 2*q14*q14+ 2*q15*q15+ 2*q16*q16+ 2*q17*q17+ 2*q18*q18- 2*q28*q28- 2*q29*q29- 2*q30*q30+ 2*p1*q1 + 2*p2*q2 + 2*p3*q3 - 2*p13*q13 - 2*p14*q14 - 2*p15*q15 - 2*p16*q16 - 2*p17*q17 - 2*p18*q18 + 2*p28*q28 + 2*p29*q29 + 2*p30*q30;
    h_phc_coeffs_Ht[35]=2*p1*p1- 4*p1*q1 + 2*p2*p2- 4*p2*q2 + 2*p3*p3- 4*p3*q3 - 2*p13*p13+ 4*p13*q13 - 2*p14*p14+ 4*p14*q14 - 2*p15*p15+ 4*p15*q15 - 2*p16*p16+ 4*p16*q16 - 2*p17*p17+ 4*p17*q17 - 2*p18*p18+ 4*p18*q18 + 2*p28*p28- 4*p28*q28 + 2*p29*p29- 4*p29*q29 + 2*p30*p30- 4*p30*q30 + 2*q1*q1+ 2*q2*q2+ 2*q3*q3- 2*q13*q13- 2*q14*q14- 2*q15*q15- 2*q16*q16- 2*q17*q17- 2*q18*q18+ 2*q28*q28+ 2*q29*q29+ 2*q30*q30;
    h_phc_coeffs_Ht[36]=2*q4*q4- 2*q2*q2- 2*q3*q3- 2*q1*q1+ 2*q5*q5+ 2*q6*q6+ 2*q16*q16+ 2*q17*q17+ 2*q18*q18- 2*q19*q19- 2*q20*q20- 2*q21*q21+ 2*p1*q1 + 2*p2*q2 + 2*p3*q3 - 2*p4*q4 - 2*p5*q5 - 2*p6*q6 - 2*p16*q16 - 2*p17*q17 - 2*p18*q18 + 2*p19*q19 + 2*p20*q20 + 2*p21*q21;
    h_phc_coeffs_Ht[37]=2*p1*p1- 4*p1*q1 + 2*p2*p2- 4*p2*q2 + 2*p3*p3- 4*p3*q3 - 2*p4*p4+ 4*p4*q4 - 2*p5*p5+ 4*p5*q5 - 2*p6*p6+ 4*p6*q6 - 2*p16*p16+ 4*p16*q16 - 2*p17*p17+ 4*p17*q17 - 2*p18*p18+ 4*p18*q18 + 2*p19*p19- 4*p19*q19 + 2*p20*p20- 4*p20*q20 + 2*p21*p21- 4*p21*q21 + 2*q1*q1+ 2*q2*q2+ 2*q3*q3- 2*q4*q4- 2*q5*q5- 2*q6*q6- 2*q16*q16- 2*q17*q17- 2*q18*q18+ 2*q19*q19+ 2*q20*q20+ 2*q21*q21;
    h_phc_coeffs_Ht[38]=2*q7*q7- 2*q2*q2- 2*q3*q3- 2*q1*q1+ 2*q8*q8+ 2*q9*q9+ 2*q16*q16+ 2*q17*q17+ 2*q18*q18- 2*q22*q22- 2*q23*q23- 2*q24*q24+ 2*p1*q1 + 2*p2*q2 + 2*p3*q3 - 2*p7*q7 - 2*p8*q8 - 2*p9*q9 - 2*p16*q16 - 2*p17*q17 - 2*p18*q18 + 2*p22*q22 + 2*p23*q23 + 2*p24*q24;
    h_phc_coeffs_Ht[39]=2*p1*p1- 4*p1*q1 + 2*p2*p2- 4*p2*q2 + 2*p3*p3- 4*p3*q3 - 2*p7*p7+ 4*p7*q7 - 2*p8*p8+ 4*p8*q8 - 2*p9*p9+ 4*p9*q9 - 2*p16*p16+ 4*p16*q16 - 2*p17*p17+ 4*p17*q17 - 2*p18*p18+ 4*p18*q18 + 2*p22*p22- 4*p22*q22 + 2*p23*p23- 4*p23*q23 + 2*p24*p24- 4*p24*q24 + 2*q1*q1+ 2*q2*q2+ 2*q3*q3- 2*q7*q7- 2*q8*q8- 2*q9*q9- 2*q16*q16- 2*q17*q17- 2*q18*q18+ 2*q22*q22+ 2*q23*q23+ 2*q24*q24;
    h_phc_coeffs_Ht[40]=p2 - p11 - q2 + q11;
    h_phc_coeffs_Ht[41]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[42]=p2 - p14 - q2 + q14;
    h_phc_coeffs_Ht[43]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[44]=p2 - p5 - q2 + q5;
    h_phc_coeffs_Ht[45]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[46]=p2 - p8 - q2 + q8;
    h_phc_coeffs_Ht[47]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[48]=p3 - p12 - q3 + q12 + p1*q31 + p31*q1 + p2*q32 + p32*q2 - p10*q31 - p31*q10 - p11*q32 - p32*q11 - 2*q1*q31 - 2*q2*q32 + 2*q10*q31 + 2*q11*q32;
    h_phc_coeffs_Ht[49]=2*p1*p31 + 2*p2*p32 - 2*p10*p31 - 2*p11*p32 - 2*p1*q31 - 2*p31*q1 - 2*p2*q32 - 2*p32*q2 + 2*p10*q31 + 2*p31*q10 + 2*p11*q32 + 2*p32*q11 + 2*q1*q31 + 2*q2*q32 - 2*q10*q31 - 2*q11*q32;
    h_phc_coeffs_Ht[50]=p3 - p15 - q3 + q15 + p1*q31 + p31*q1 + p2*q32 + p32*q2 - p13*q31 - p31*q13 - p14*q32 - p32*q14 - 2*q1*q31 - 2*q2*q32 + 2*q13*q31 + 2*q14*q32;
    h_phc_coeffs_Ht[51]=2*p1*p31 + 2*p2*p32 - 2*p13*p31 - 2*p14*p32 - 2*p1*q31 - 2*p31*q1 - 2*p2*q32 - 2*p32*q2 + 2*p13*q31 + 2*p31*q13 + 2*p14*q32 + 2*p32*q14 + 2*q1*q31 + 2*q2*q32 - 2*q13*q31 - 2*q14*q32;
    h_phc_coeffs_Ht[52]=p3 - p6 - q3 + q6 + p1*q31 + p31*q1 + p2*q32 + p32*q2 - p4*q31 - p31*q4 - p5*q32 - p32*q5 - 2*q1*q31 - 2*q2*q32 + 2*q4*q31 + 2*q5*q32;
    h_phc_coeffs_Ht[53]=2*p1*p31 + 2*p2*p32 - 2*p4*p31 - 2*p5*p32 - 2*p1*q31 - 2*p31*q1 - 2*p2*q32 - 2*p32*q2 + 2*p4*q31 + 2*p31*q4 + 2*p5*q32 + 2*p32*q5 + 2*q1*q31 + 2*q2*q32 - 2*q4*q31 - 2*q5*q32;
    h_phc_coeffs_Ht[54]=p3 - p9 - q3 + q9 + p1*q31 + p31*q1 + p2*q32 + p32*q2 - p7*q31 - p31*q7 - p8*q32 - p32*q8 - 2*q1*q31 - 2*q2*q32 + 2*q7*q31 + 2*q8*q32;
    h_phc_coeffs_Ht[55]=2*p1*p31 + 2*p2*p32 - 2*p7*p31 - 2*p8*p32 - 2*p1*q31 - 2*p31*q1 - 2*p2*q32 - 2*p32*q2 + 2*p7*q31 + 2*p31*q7 + 2*p8*q32 + 2*p32*q8 + 2*q1*q31 + 2*q2*q32 - 2*q7*q31 - 2*q8*q32;
    h_phc_coeffs_Ht[56]=MAGMA_C_ONE;
    h_phc_coeffs_Ht[57]=MAGMA_C_ZERO;

  }
} // end of namespace

#endif
