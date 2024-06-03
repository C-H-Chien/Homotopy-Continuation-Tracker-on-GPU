#ifndef P2C_UNCALIBRATED_TRIFOCAL_REL_POS_CY_H
#define P2C_UNCALIBRATED_TRIFOCAL_REL_POS_CY_H
// ==========================================================================================================================
//
// Modifications
//    Chiang-Heng Chien  24-05-26:   Initially created. Built on top of generalized three-view cameras minimal problems.
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

  void p2c_uncalibrated_trifocal_rel_pos_CY(
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

    h_phc_coeffs_Hx[0]=q1*q2 + q5*q6;
    h_phc_coeffs_Hx[1]=p1*q2 + p2*q1 + p5*q6 + p6*q5 - 2*q1*q2 - 2*q5*q6;
    h_phc_coeffs_Hx[2]=p1*p2 + p5*p6 - p1*q2 - p2*q1 - p5*q6 - p6*q5 + q1*q2 + q5*q6;
    h_phc_coeffs_Hx[3]=q1*q3 + q5*q7;
    h_phc_coeffs_Hx[4]=p1*q3 + p3*q1 + p5*q7 + p7*q5 - 2*q1*q3 - 2*q5*q7;
    h_phc_coeffs_Hx[5]=p1*p3 + p5*p7 - p1*q3 - p3*q1 - p5*q7 - p7*q5 + q1*q3 + q5*q7;
    h_phc_coeffs_Hx[6]=q1*q4 + q5*q8;
    h_phc_coeffs_Hx[7]=p1*q4 + p4*q1 + p5*q8 + p8*q5 - 2*q1*q4 - 2*q5*q8;
    h_phc_coeffs_Hx[8]=p1*p4 + p5*p8 - p1*q4 - p4*q1 - p5*q8 - p8*q5 + q1*q4 + q5*q8;
    h_phc_coeffs_Hx[9]=q10*q11 + q14*q15;
    h_phc_coeffs_Hx[10]=p10*q11 + p11*q10 + p14*q15 + p15*q14 - 2*q10*q11 - 2*q14*q15;
    h_phc_coeffs_Hx[11]=p10*p11 + p14*p15 - p10*q11 - p11*q10 - p14*q15 - p15*q14 + q10*q11 + q14*q15;
    h_phc_coeffs_Hx[12]=q10*q12 + q14*q16;
    h_phc_coeffs_Hx[13]=p10*q12 + p12*q10 + p14*q16 + p16*q14 - 2*q10*q12 - 2*q14*q16;
    h_phc_coeffs_Hx[14]=p10*p12 + p14*p16 - p10*q12 - p12*q10 - p14*q16 - p16*q14 + q10*q12 + q14*q16;
    h_phc_coeffs_Hx[15]=q10*q10+q14*q14;
    h_phc_coeffs_Hx[16]=2*p10*q10 + 2*p14*q14 - 2*q10*q10- 2*q14*q14;
    h_phc_coeffs_Hx[17]=p10*p10- 2*p14*q14 - 2*p10*q10 +p14*p14+q10*q10+q14*q14;
    h_phc_coeffs_Hx[18]=q11*q12 + q15*q16;
    h_phc_coeffs_Hx[19]=p11*q12 + p12*q11 + p15*q16 + p16*q15 - 2*q11*q12 - 2*q15*q16;
    h_phc_coeffs_Hx[20]=p11*p12 + p15*p16 - p11*q12 - p12*q11 - p15*q16 - p16*q15 + q11*q12 + q15*q16;
    h_phc_coeffs_Hx[21]=q11*q11+q15*q15;
    h_phc_coeffs_Hx[22]=2*p11*q11 + 2*p15*q15 - 2*q11*q11- 2*q15*q15;
    h_phc_coeffs_Hx[23]=p11*p11- 2*p15*q15 - 2*p11*q11 +p15*p15+q11*q11+q15*q15;
    h_phc_coeffs_Hx[24]=q12*q12+q16*q16;
    h_phc_coeffs_Hx[25]=2*p12*q12 + 2*p16*q16 - 2*q12*q12- 2*q16*q16;
    h_phc_coeffs_Hx[26]=p12*p12- 2*p16*q16 - 2*p12*q12 +p16*p16+q12*q12+q16*q16;
    h_phc_coeffs_Hx[27]=q17*q18 + q21*q22;
    h_phc_coeffs_Hx[28]=p17*q18 + p18*q17 + p21*q22 + p22*q21 - 2*q17*q18 - 2*q21*q22;
    h_phc_coeffs_Hx[29]=p17*p18 + p21*p22 - p17*q18 - p18*q17 - p21*q22 - p22*q21 + q17*q18 + q21*q22;
    h_phc_coeffs_Hx[30]=q17*q19 + q21*q23;
    h_phc_coeffs_Hx[31]=p17*q19 + p19*q17 + p21*q23 + p23*q21 - 2*q17*q19 - 2*q21*q23;
    h_phc_coeffs_Hx[32]=p17*p19 + p21*p23 - p17*q19 - p19*q17 - p21*q23 - p23*q21 + q17*q19 + q21*q23;
    h_phc_coeffs_Hx[33]=q17*q20 + q21*q24;
    h_phc_coeffs_Hx[34]=p17*q20 + p20*q17 + p21*q24 + p24*q21 - 2*q17*q20 - 2*q21*q24;
    h_phc_coeffs_Hx[35]=p17*p20 + p21*p24 - p17*q20 - p20*q17 - p21*q24 - p24*q21 + q17*q20 + q21*q24;
    h_phc_coeffs_Hx[36]=q17*q17+q21*q21;
    h_phc_coeffs_Hx[37]=2*p17*q17 + 2*p21*q21 - 2*q17*q17- 2*q21*q21;
    h_phc_coeffs_Hx[38]=p17*p17- 2*p21*q21 - 2*p17*q17 +p21*p21+q17*q17+q21*q21;
    h_phc_coeffs_Hx[39]=q18*q19 + q22*q23;
    h_phc_coeffs_Hx[40]=p18*q19 + p19*q18 + p22*q23 + p23*q22 - 2*q18*q19 - 2*q22*q23;
    h_phc_coeffs_Hx[41]=p18*p19 + p22*p23 - p18*q19 - p19*q18 - p22*q23 - p23*q22 + q18*q19 + q22*q23;
    h_phc_coeffs_Hx[42]=q18*q20 + q22*q24;
    h_phc_coeffs_Hx[43]=p18*q20 + p20*q18 + p22*q24 + p24*q22 - 2*q18*q20 - 2*q22*q24;
    h_phc_coeffs_Hx[44]=p18*p20 + p22*p24 - p18*q20 - p20*q18 - p22*q24 - p24*q22 + q18*q20 + q22*q24;
    h_phc_coeffs_Hx[45]=q18*q18+q22*q22;
    h_phc_coeffs_Hx[46]=2*p18*q18 + 2*p22*q22 - 2*q18*q18- 2*q22*q22;
    h_phc_coeffs_Hx[47]=p18*p18- 2*p22*q22 - 2*p18*q18 +p22*p22+q18*q18+q22*q22;
    h_phc_coeffs_Hx[48]=q19*q20 + q23*q24;
    h_phc_coeffs_Hx[49]=p19*q20 + p20*q19 + p23*q24 + p24*q23 - 2*q19*q20 - 2*q23*q24;
    h_phc_coeffs_Hx[50]=p19*p20 + p23*p24 - p19*q20 - p20*q19 - p23*q24 - p24*q23 + q19*q20 + q23*q24;
    h_phc_coeffs_Hx[51]=q19*q19+q23*q23;
    h_phc_coeffs_Hx[52]=2*p19*q19 + 2*p23*q23 - 2*q19*q19- 2*q23*q23;
    h_phc_coeffs_Hx[53]=p19*p19- 2*p23*q23 - 2*p19*q19 +p23*p23+q19*q19+q23*q23;
    h_phc_coeffs_Hx[54]=q1*q1+q5*q5;
    h_phc_coeffs_Hx[55]=2*p1*q1 + 2*p5*q5 - 2*q1*q1- 2*q5*q5;
    h_phc_coeffs_Hx[56]=p1*p1- 2*p5*q5 - 2*p1*q1 +p5*p5+q1*q1+q5*q5;
    h_phc_coeffs_Hx[57]=q2*q3 + q6*q7;
    h_phc_coeffs_Hx[58]=p2*q3 + p3*q2 + p6*q7 + p7*q6 - 2*q2*q3 - 2*q6*q7;
    h_phc_coeffs_Hx[59]=p2*p3 + p6*p7 - p2*q3 - p3*q2 - p6*q7 - p7*q6 + q2*q3 + q6*q7;
    h_phc_coeffs_Hx[60]=q2*q4 + q6*q8;
    h_phc_coeffs_Hx[61]=p2*q4 + p4*q2 + p6*q8 + p8*q6 - 2*q2*q4 - 2*q6*q8;
    h_phc_coeffs_Hx[62]=p2*p4 + p6*p8 - p2*q4 - p4*q2 - p6*q8 - p8*q6 + q2*q4 + q6*q8;
    h_phc_coeffs_Hx[63]=q20*q20+q24*q24;
    h_phc_coeffs_Hx[64]=2*p20*q20 + 2*p24*q24 - 2*q20*q20- 2*q24*q24;
    h_phc_coeffs_Hx[65]=p20*p20- 2*p24*q24 - 2*p20*q20 +p24*p24+q20*q20+q24*q24;
    h_phc_coeffs_Hx[66]=q2*q2+q6*q6;
    h_phc_coeffs_Hx[67]=2*p2*q2 + 2*p6*q6 - 2*q2*q2- 2*q6*q6;
    h_phc_coeffs_Hx[68]=p2*p2- 2*p6*q6 - 2*p2*q2 +p6*p6+q2*q2+q6*q6;
    h_phc_coeffs_Hx[69]=q3*q4 + q7*q8;
    h_phc_coeffs_Hx[70]=p3*q4 + p4*q3 + p7*q8 + p8*q7 - 2*q3*q4 - 2*q7*q8;
    h_phc_coeffs_Hx[71]=p3*p4 + p7*p8 - p3*q4 - p4*q3 - p7*q8 - p8*q7 + q3*q4 + q7*q8;
    h_phc_coeffs_Hx[72]=q3*q3+q7*q7;
    h_phc_coeffs_Hx[73]=2*p3*q3 + 2*p7*q7 - 2*q3*q3- 2*q7*q7;
    h_phc_coeffs_Hx[74]=p3*p3- 2*p7*q7 - 2*p3*q3 +p7*p7+q3*q3+q7*q7;
    h_phc_coeffs_Hx[75]=q4*q4+q8*q8;
    h_phc_coeffs_Hx[76]=2*p4*q4 + 2*p8*q8 - 2*q4*q4- 2*q8*q8;
    h_phc_coeffs_Hx[77]=p4*p4- 2*p8*q8 - 2*p4*q4 +p8*p8+q4*q4+q8*q8;
    h_phc_coeffs_Hx[78]=q9*q10 + q13*q14;
    h_phc_coeffs_Hx[79]=p9*q10 + p10*q9 + p13*q14 + p14*q13 - 2*q9*q10 - 2*q13*q14;
    h_phc_coeffs_Hx[80]=p9*p10 + p13*p14 - p9*q10 - p10*q9 - p13*q14 - p14*q13 + q9*q10 + q13*q14;
    h_phc_coeffs_Hx[81]=q9*q11 + q13*q15;
    h_phc_coeffs_Hx[82]=p9*q11 + p11*q9 + p13*q15 + p15*q13 - 2*q9*q11 - 2*q13*q15;
    h_phc_coeffs_Hx[83]=p9*p11 + p13*p15 - p9*q11 - p11*q9 - p13*q15 - p15*q13 + q9*q11 + q13*q15;
    h_phc_coeffs_Hx[84]=q9*q12 + q13*q16;
    h_phc_coeffs_Hx[85]=p9*q12 + p12*q9 + p13*q16 + p16*q13 - 2*q9*q12 - 2*q13*q16;
    h_phc_coeffs_Hx[86]=p9*p12 + p13*p16 - p9*q12 - p12*q9 - p13*q16 - p16*q13 + q9*q12 + q13*q16;
    h_phc_coeffs_Hx[87]=q9*q9+q13*q13;
    h_phc_coeffs_Hx[88]=2*p9*q9 + 2*p13*q13 - 2*q9*q9- 2*q13*q13;
    h_phc_coeffs_Hx[89]=p9*p9- 2*p13*q13 - 2*p9*q9 +p13*p13+q9*q9+q13*q13;
    h_phc_coeffs_Hx[90]=MAGMA_C_ONE;
    h_phc_coeffs_Hx[91]=MAGMA_C_ZERO;
    h_phc_coeffs_Hx[92]=MAGMA_C_ZERO;

    h_phc_coeffs_Ht[0]=p1*q2 + p2*q1 + p5*q6 + p6*q5 - 2*q1*q2 - 2*q5*q6;
    h_phc_coeffs_Ht[1]=2*p1*p2 + 2*p5*p6 - 2*p1*q2 - 2*p2*q1 - 2*p5*q6 - 2*p6*q5 + 2*q1*q2 + 2*q5*q6;
    h_phc_coeffs_Ht[2]=p1*q3 + p3*q1 + p5*q7 + p7*q5 - 2*q1*q3 - 2*q5*q7;
    h_phc_coeffs_Ht[3]=2*p1*p3 + 2*p5*p7 - 2*p1*q3 - 2*p3*q1 - 2*p5*q7 - 2*p7*q5 + 2*q1*q3 + 2*q5*q7;
    h_phc_coeffs_Ht[4]=p1*q4 + p4*q1 + p5*q8 + p8*q5 - 2*q1*q4 - 2*q5*q8;
    h_phc_coeffs_Ht[5]=2*p1*p4 + 2*p5*p8 - 2*p1*q4 - 2*p4*q1 - 2*p5*q8 - 2*p8*q5 + 2*q1*q4 + 2*q5*q8;
    h_phc_coeffs_Ht[6]=p10*q11 + p11*q10 + p14*q15 + p15*q14 - 2*q10*q11 - 2*q14*q15;
    h_phc_coeffs_Ht[7]=2*p10*p11 + 2*p14*p15 - 2*p10*q11 - 2*p11*q10 - 2*p14*q15 - 2*p15*q14 + 2*q10*q11 + 2*q14*q15;
    h_phc_coeffs_Ht[8]=p10*q12 + p12*q10 + p14*q16 + p16*q14 - 2*q10*q12 - 2*q14*q16;
    h_phc_coeffs_Ht[9]=2*p10*p12 + 2*p14*p16 - 2*p10*q12 - 2*p12*q10 - 2*p14*q16 - 2*p16*q14 + 2*q10*q12 + 2*q14*q16;
    h_phc_coeffs_Ht[10]=2*p10*q10 - 2*q14*q14- 2*q10*q10+ 2*p14*q14;
    h_phc_coeffs_Ht[11]=2*p10*p10- 4*p10*q10 + 2*p14*p14- 4*p14*q14 + 2*q10*q10+ 2*q14*q14;
    h_phc_coeffs_Ht[12]=p11*q12 + p12*q11 + p15*q16 + p16*q15 - 2*q11*q12 - 2*q15*q16;
    h_phc_coeffs_Ht[13]=2*p11*p12 + 2*p15*p16 - 2*p11*q12 - 2*p12*q11 - 2*p15*q16 - 2*p16*q15 + 2*q11*q12 + 2*q15*q16;
    h_phc_coeffs_Ht[14]=2*p11*q11 - 2*q15*q15- 2*q11*q11+ 2*p15*q15;
    h_phc_coeffs_Ht[15]=2*p11*p11- 4*p11*q11 + 2*p15*p15- 4*p15*q15 + 2*q11*q11+ 2*q15*q15;
    h_phc_coeffs_Ht[16]=2*p12*q12 - 2*q16*q16- 2*q12*q12+ 2*p16*q16;
    h_phc_coeffs_Ht[17]=2*p12*p12- 4*p12*q12 + 2*p16*p16- 4*p16*q16 + 2*q12*q12+ 2*q16*q16;
    h_phc_coeffs_Ht[18]=p17*q18 + p18*q17 + p21*q22 + p22*q21 - 2*q17*q18 - 2*q21*q22;
    h_phc_coeffs_Ht[19]=2*p17*p18 + 2*p21*p22 - 2*p17*q18 - 2*p18*q17 - 2*p21*q22 - 2*p22*q21 + 2*q17*q18 + 2*q21*q22;
    h_phc_coeffs_Ht[20]=p17*q19 + p19*q17 + p21*q23 + p23*q21 - 2*q17*q19 - 2*q21*q23;
    h_phc_coeffs_Ht[21]=2*p17*p19 + 2*p21*p23 - 2*p17*q19 - 2*p19*q17 - 2*p21*q23 - 2*p23*q21 + 2*q17*q19 + 2*q21*q23;
    h_phc_coeffs_Ht[22]=p17*q20 + p20*q17 + p21*q24 + p24*q21 - 2*q17*q20 - 2*q21*q24;
    h_phc_coeffs_Ht[23]=2*p17*p20 + 2*p21*p24 - 2*p17*q20 - 2*p20*q17 - 2*p21*q24 - 2*p24*q21 + 2*q17*q20 + 2*q21*q24;
    h_phc_coeffs_Ht[24]=2*p17*q17 - 2*q21*q21- 2*q17*q17+ 2*p21*q21;
    h_phc_coeffs_Ht[25]=2*p17*p17- 4*p17*q17 + 2*p21*p21- 4*p21*q21 + 2*q17*q17+ 2*q21*q21;
    h_phc_coeffs_Ht[26]=p18*q19 + p19*q18 + p22*q23 + p23*q22 - 2*q18*q19 - 2*q22*q23;
    h_phc_coeffs_Ht[27]=2*p18*p19 + 2*p22*p23 - 2*p18*q19 - 2*p19*q18 - 2*p22*q23 - 2*p23*q22 + 2*q18*q19 + 2*q22*q23;
    h_phc_coeffs_Ht[28]=p18*q20 + p20*q18 + p22*q24 + p24*q22 - 2*q18*q20 - 2*q22*q24;
    h_phc_coeffs_Ht[29]=2*p18*p20 + 2*p22*p24 - 2*p18*q20 - 2*p20*q18 - 2*p22*q24 - 2*p24*q22 + 2*q18*q20 + 2*q22*q24;
    h_phc_coeffs_Ht[30]=2*p18*q18 - 2*q22*q22- 2*q18*q18+ 2*p22*q22;
    h_phc_coeffs_Ht[31]=2*p18*p18- 4*p18*q18 + 2*p22*p22- 4*p22*q22 + 2*q18*q18+ 2*q22*q22;
    h_phc_coeffs_Ht[32]=p19*q20 + p20*q19 + p23*q24 + p24*q23 - 2*q19*q20 - 2*q23*q24;
    h_phc_coeffs_Ht[33]=2*p19*p20 + 2*p23*p24 - 2*p19*q20 - 2*p20*q19 - 2*p23*q24 - 2*p24*q23 + 2*q19*q20 + 2*q23*q24;
    h_phc_coeffs_Ht[34]=2*p19*q19 - 2*q23*q23- 2*q19*q19+ 2*p23*q23;
    h_phc_coeffs_Ht[35]=2*p19*p19- 4*p19*q19 + 2*p23*p23- 4*p23*q23 + 2*q19*q19+ 2*q23*q23;
    h_phc_coeffs_Ht[36]=2*p1*q1 - 2*q5*q5- 2*q1*q1+ 2*p5*q5;
    h_phc_coeffs_Ht[37]=2*p1*p1- 4*p1*q1 + 2*p5*p5- 4*p5*q5 + 2*q1*q1+ 2*q5*q5;
    h_phc_coeffs_Ht[38]=p2*q3 + p3*q2 + p6*q7 + p7*q6 - 2*q2*q3 - 2*q6*q7;
    h_phc_coeffs_Ht[39]=2*p2*p3 + 2*p6*p7 - 2*p2*q3 - 2*p3*q2 - 2*p6*q7 - 2*p7*q6 + 2*q2*q3 + 2*q6*q7;
    h_phc_coeffs_Ht[40]=p2*q4 + p4*q2 + p6*q8 + p8*q6 - 2*q2*q4 - 2*q6*q8;
    h_phc_coeffs_Ht[41]=2*p2*p4 + 2*p6*p8 - 2*p2*q4 - 2*p4*q2 - 2*p6*q8 - 2*p8*q6 + 2*q2*q4 + 2*q6*q8;
    h_phc_coeffs_Ht[42]=2*p20*q20 - 2*q24*q24- 2*q20*q20+ 2*p24*q24;
    h_phc_coeffs_Ht[43]=2*p20*p20- 4*p20*q20 + 2*p24*p24- 4*p24*q24 + 2*q20*q20+ 2*q24*q24;
    h_phc_coeffs_Ht[44]=2*p2*q2 - 2*q6*q6- 2*q2*q2+ 2*p6*q6;
    h_phc_coeffs_Ht[45]=2*p2*p2- 4*p2*q2 + 2*p6*p6- 4*p6*q6 + 2*q2*q2+ 2*q6*q6;
    h_phc_coeffs_Ht[46]=p3*q4 + p4*q3 + p7*q8 + p8*q7 - 2*q3*q4 - 2*q7*q8;
    h_phc_coeffs_Ht[47]=2*p3*p4 + 2*p7*p8 - 2*p3*q4 - 2*p4*q3 - 2*p7*q8 - 2*p8*q7 + 2*q3*q4 + 2*q7*q8;
    h_phc_coeffs_Ht[48]=2*p3*q3 - 2*q7*q7- 2*q3*q3+ 2*p7*q7;
    h_phc_coeffs_Ht[49]=2*p3*p3- 4*p3*q3 + 2*p7*p7- 4*p7*q7 + 2*q3*q3+ 2*q7*q7;
    h_phc_coeffs_Ht[50]=2*p4*q4 - 2*q8*q8- 2*q4*q4+ 2*p8*q8;
    h_phc_coeffs_Ht[51]=2*p4*p4- 4*p4*q4 + 2*p8*p8- 4*p8*q8 + 2*q4*q4+ 2*q8*q8;
    h_phc_coeffs_Ht[52]=p9*q10 + p10*q9 + p13*q14 + p14*q13 - 2*q9*q10 - 2*q13*q14;
    h_phc_coeffs_Ht[53]=2*p9*p10 + 2*p13*p14 - 2*p9*q10 - 2*p10*q9 - 2*p13*q14 - 2*p14*q13 + 2*q9*q10 + 2*q13*q14;
    h_phc_coeffs_Ht[54]=p9*q11 + p11*q9 + p13*q15 + p15*q13 - 2*q9*q11 - 2*q13*q15;
    h_phc_coeffs_Ht[55]=2*p9*p11 + 2*p13*p15 - 2*p9*q11 - 2*p11*q9 - 2*p13*q15 - 2*p15*q13 + 2*q9*q11 + 2*q13*q15;
    h_phc_coeffs_Ht[56]=p9*q12 + p12*q9 + p13*q16 + p16*q13 - 2*q9*q12 - 2*q13*q16;
    h_phc_coeffs_Ht[57]=2*p9*p12 + 2*p13*p16 - 2*p9*q12 - 2*p12*q9 - 2*p13*q16 - 2*p16*q13 + 2*q9*q12 + 2*q13*q16;
    h_phc_coeffs_Ht[58]=2*p9*q9 - 2*q13*q13- 2*q9*q9+ 2*p13*q13;
    h_phc_coeffs_Ht[59]=2*p9*p9- 4*p9*q9 + 2*p13*p13- 4*p13*q13 + 2*q9*q9+ 2*q13*q13;
    h_phc_coeffs_Ht[60]=MAGMA_C_ONE;
    h_phc_coeffs_Ht[61]=MAGMA_C_ZERO;

  }
} // end of namespace

#endif
