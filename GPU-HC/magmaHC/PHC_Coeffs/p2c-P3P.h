#ifndef P2C_P3P_H
#define P2C_P3P_H
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

  void p2c_P3P(
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

    h_phc_coeffs_Hx[0]=q1*q4 + q2*q5 + q3*q6;
    h_phc_coeffs_Hx[1]=p1*q4 + p4*q1 + p2*q5 + p5*q2 + p3*q6 + p6*q3 - 2*q1*q4 - 2*q2*q5 - 2*q3*q6;
    h_phc_coeffs_Hx[2]=p1*p4 + p2*p5 + p3*p6 - p1*q4 - p4*q1 - p2*q5 - p5*q2 - p3*q6 - p6*q3 + q1*q4 + q2*q5 + q3*q6;
    h_phc_coeffs_Hx[3]=q1*q7 + q2*q8 + q3*q9;
    h_phc_coeffs_Hx[4]=p1*q7 + p7*q1 + p2*q8 + p8*q2 + p3*q9 + p9*q3 - 2*q1*q7 - 2*q2*q8 - 2*q3*q9;
    h_phc_coeffs_Hx[5]=p1*p7 + p2*p8 + p3*p9 - p1*q7 - p7*q1 - p2*q8 - p8*q2 - p3*q9 - p9*q3 + q1*q7 + q2*q8 + q3*q9;
    h_phc_coeffs_Hx[6]=q10*q10- 2*q11*q14 - 2*q12*q15 - 2*q10*q13 +q11*q11+q12*q12+q13*q13+q14*q14+q15*q15;
    h_phc_coeffs_Hx[7]=2*p10*q10 + 2*p11*q11 - 2*p10*q13 - 2*p13*q10 + 2*p12*q12 - 2*p11*q14 - 2*p14*q11 + 2*p13*q13 - 2*p12*q15 - 2*p15*q12 + 2*p14*q14 + 2*p15*q15 + 4*q10*q13 + 4*q11*q14 + 4*q12*q15 - 2*q10*q10- 2*q11*q11- 2*q12*q12- 2*q13*q13- 2*q14*q14- 2*q15*q15;
    h_phc_coeffs_Hx[8]=2*p10*q13 - 2*p11*p14 - 2*p12*p15 - 2*p10*q10 - 2*p11*q11 - 2*p10*p13 + 2*p13*q10 - 2*p12*q12 + 2*p11*q14 + 2*p14*q11 - 2*p13*q13 + 2*p12*q15 + 2*p15*q12 - 2*p14*q14 - 2*p15*q15 - 2*q10*q13 - 2*q11*q14 - 2*q12*q15 +p10*p10+p11*p11+p12*p12+p13*p13+p14*p14+p15*p15+q10*q10+q11*q11+q12*q12+q13*q13+q14*q14+q15*q15;
    h_phc_coeffs_Hx[9]=q10*q10- 2*q11*q17 - 2*q12*q18 - 2*q10*q16 +q11*q11+q12*q12+q16*q16+q17*q17+q18*q18;
    h_phc_coeffs_Hx[10]=2*p10*q10 + 2*p11*q11 + 2*p12*q12 - 2*p10*q16 - 2*p16*q10 - 2*p11*q17 - 2*p17*q11 - 2*p12*q18 - 2*p18*q12 + 2*p16*q16 + 2*p17*q17 + 2*p18*q18 + 4*q10*q16 + 4*q11*q17 + 4*q12*q18 - 2*q10*q10- 2*q11*q11- 2*q12*q12- 2*q16*q16- 2*q17*q17- 2*q18*q18;
    h_phc_coeffs_Hx[11]=2*p10*q16 - 2*p11*p17 - 2*p12*p18 - 2*p10*q10 - 2*p11*q11 - 2*p12*q12 - 2*p10*p16 + 2*p16*q10 + 2*p11*q17 + 2*p17*q11 + 2*p12*q18 + 2*p18*q12 - 2*p16*q16 - 2*p17*q17 - 2*p18*q18 - 2*q10*q16 - 2*q11*q17 - 2*q12*q18 +p10*p10+p11*p11+p12*p12+p16*p16+p17*p17+p18*p18+q10*q10+q11*q11+q12*q12+q16*q16+q17*q17+q18*q18;
    h_phc_coeffs_Hx[12]=q13*q13- 2*q14*q17 - 2*q15*q18 - 2*q13*q16 +q14*q14+q15*q15+q16*q16+q17*q17+q18*q18;
    h_phc_coeffs_Hx[13]=2*p13*q13 + 2*p14*q14 - 2*p13*q16 - 2*p16*q13 + 2*p15*q15 - 2*p14*q17 - 2*p17*q14 + 2*p16*q16 - 2*p15*q18 - 2*p18*q15 + 2*p17*q17 + 2*p18*q18 + 4*q13*q16 + 4*q14*q17 + 4*q15*q18 - 2*q13*q13- 2*q14*q14- 2*q15*q15- 2*q16*q16- 2*q17*q17- 2*q18*q18;
    h_phc_coeffs_Hx[14]=2*p13*q16 - 2*p14*p17 - 2*p15*p18 - 2*p13*q13 - 2*p14*q14 - 2*p13*p16 + 2*p16*q13 - 2*p15*q15 + 2*p14*q17 + 2*p17*q14 - 2*p16*q16 + 2*p15*q18 + 2*p18*q15 - 2*p17*q17 - 2*p18*q18 - 2*q13*q16 - 2*q14*q17 - 2*q15*q18 +p13*p13+p14*p14+p15*p15+p16*p16+p17*p17+p18*p18+q13*q13+q14*q14+q15*q15+q16*q16+q17*q17+q18*q18;
    h_phc_coeffs_Hx[15]=q1*q1+q2*q2+q3*q3;
    h_phc_coeffs_Hx[16]=2*p1*q1 + 2*p2*q2 + 2*p3*q3 - 2*q1*q1- 2*q2*q2- 2*q3*q3;
    h_phc_coeffs_Hx[17]=p1*p1- 2*p2*q2 - 2*p3*q3 - 2*p1*q1 +p2*p2+p3*p3+q1*q1+q2*q2+q3*q3;
    h_phc_coeffs_Hx[18]=q4*q7 + q5*q8 + q6*q9;
    h_phc_coeffs_Hx[19]=p4*q7 + p7*q4 + p5*q8 + p8*q5 + p6*q9 + p9*q6 - 2*q4*q7 - 2*q5*q8 - 2*q6*q9;
    h_phc_coeffs_Hx[20]=p4*p7 + p5*p8 + p6*p9 - p4*q7 - p7*q4 - p5*q8 - p8*q5 - p6*q9 - p9*q6 + q4*q7 + q5*q8 + q6*q9;
    h_phc_coeffs_Hx[21]=q4*q4+q5*q5+q6*q6;
    h_phc_coeffs_Hx[22]=2*p4*q4 + 2*p5*q5 + 2*p6*q6 - 2*q4*q4- 2*q5*q5- 2*q6*q6;
    h_phc_coeffs_Hx[23]=p4*p4- 2*p5*q5 - 2*p6*q6 - 2*p4*q4 +p5*p5+p6*p6+q4*q4+q5*q5+q6*q6;
    h_phc_coeffs_Hx[24]=q7*q7+q8*q8+q9*q9;
    h_phc_coeffs_Hx[25]=2*p7*q7 + 2*p8*q8 + 2*p9*q9 - 2*q7*q7- 2*q8*q8- 2*q9*q9;
    h_phc_coeffs_Hx[26]=p7*p7- 2*p8*q8 - 2*p9*q9 - 2*p7*q7 +p8*p8+p9*p9+q7*q7+q8*q8+q9*q9;
    h_phc_coeffs_Hx[27]=MAGMA_C_ONE;
    h_phc_coeffs_Hx[28]=MAGMA_C_ZERO;
    h_phc_coeffs_Hx[29]=MAGMA_C_ZERO;

    h_phc_coeffs_Ht[0]=p1*q4 + p4*q1 + p2*q5 + p5*q2 + p3*q6 + p6*q3 - 2*q1*q4 - 2*q2*q5 - 2*q3*q6;
    h_phc_coeffs_Ht[1]=2*p1*p4 + 2*p2*p5 + 2*p3*p6 - 2*p1*q4 - 2*p4*q1 - 2*p2*q5 - 2*p5*q2 - 2*p3*q6 - 2*p6*q3 + 2*q1*q4 + 2*q2*q5 + 2*q3*q6;
    h_phc_coeffs_Ht[2]=p1*q7 + p7*q1 + p2*q8 + p8*q2 + p3*q9 + p9*q3 - 2*q1*q7 - 2*q2*q8 - 2*q3*q9;
    h_phc_coeffs_Ht[3]=2*p1*p7 + 2*p2*p8 + 2*p3*p9 - 2*p1*q7 - 2*p7*q1 - 2*p2*q8 - 2*p8*q2 - 2*p3*q9 - 2*p9*q3 + 2*q1*q7 + 2*q2*q8 + 2*q3*q9;
    h_phc_coeffs_Ht[4]=2*p10*q10 + 2*p11*q11 - 2*p10*q13 - 2*p13*q10 + 2*p12*q12 - 2*p11*q14 - 2*p14*q11 + 2*p13*q13 - 2*p12*q15 - 2*p15*q12 + 2*p14*q14 + 2*p15*q15 + 4*q10*q13 + 4*q11*q14 + 4*q12*q15 - 2*q10*q10- 2*q11*q11- 2*q12*q12- 2*q13*q13- 2*q14*q14- 2*q15*q15;
    h_phc_coeffs_Ht[5]=2*p10*p10- 4*p10*p13 - 4*p10*q10 + 4*p10*q13 + 2*p11*p11- 4*p11*p14 - 4*p11*q11 + 4*p11*q14 + 2*p12*p12- 4*p12*p15 - 4*p12*q12 + 4*p12*q15 + 2*p13*p13+ 4*p13*q10 - 4*p13*q13 + 2*p14*p14+ 4*p14*q11 - 4*p14*q14 + 2*p15*p15+ 4*p15*q12 - 4*p15*q15 + 2*q10*q10- 4*q10*q13 + 2*q11*q11- 4*q11*q14 + 2*q12*q12- 4*q12*q15 + 2*q13*q13+ 2*q14*q14+ 2*q15*q15;
    h_phc_coeffs_Ht[6]=2*p10*q10 + 2*p11*q11 + 2*p12*q12 - 2*p10*q16 - 2*p16*q10 - 2*p11*q17 - 2*p17*q11 - 2*p12*q18 - 2*p18*q12 + 2*p16*q16 + 2*p17*q17 + 2*p18*q18 + 4*q10*q16 + 4*q11*q17 + 4*q12*q18 - 2*q10*q10- 2*q11*q11- 2*q12*q12- 2*q16*q16- 2*q17*q17- 2*q18*q18;
    h_phc_coeffs_Ht[7]=2*p10*p10- 4*p10*p16 - 4*p10*q10 + 4*p10*q16 + 2*p11*p11- 4*p11*p17 - 4*p11*q11 + 4*p11*q17 + 2*p12*p12- 4*p12*p18 - 4*p12*q12 + 4*p12*q18 + 2*p16*p16+ 4*p16*q10 - 4*p16*q16 + 2*p17*p17+ 4*p17*q11 - 4*p17*q17 + 2*p18*p18+ 4*p18*q12 - 4*p18*q18 + 2*q10*q10- 4*q10*q16 + 2*q11*q11- 4*q11*q17 + 2*q12*q12- 4*q12*q18 + 2*q16*q16+ 2*q17*q17+ 2*q18*q18;
    h_phc_coeffs_Ht[8]=2*p13*q13 + 2*p14*q14 - 2*p13*q16 - 2*p16*q13 + 2*p15*q15 - 2*p14*q17 - 2*p17*q14 + 2*p16*q16 - 2*p15*q18 - 2*p18*q15 + 2*p17*q17 + 2*p18*q18 + 4*q13*q16 + 4*q14*q17 + 4*q15*q18 - 2*q13*q13- 2*q14*q14- 2*q15*q15- 2*q16*q16- 2*q17*q17- 2*q18*q18;
    h_phc_coeffs_Ht[9]=2*p13*p13- 4*p13*p16 - 4*p13*q13 + 4*p13*q16 + 2*p14*p14- 4*p14*p17 - 4*p14*q14 + 4*p14*q17 + 2*p15*p15- 4*p15*p18 - 4*p15*q15 + 4*p15*q18 + 2*p16*p16+ 4*p16*q13 - 4*p16*q16 + 2*p17*p17+ 4*p17*q14 - 4*p17*q17 + 2*p18*p18+ 4*p18*q15 - 4*p18*q18 + 2*q13*q13- 4*q13*q16 + 2*q14*q14- 4*q14*q17 + 2*q15*q15- 4*q15*q18 + 2*q16*q16+ 2*q17*q17+ 2*q18*q18;
    h_phc_coeffs_Ht[10]=2*p1*q1 - 2*q2*q2- 2*q3*q3- 2*q1*q1+ 2*p2*q2 + 2*p3*q3;
    h_phc_coeffs_Ht[11]=2*p1*p1- 4*p1*q1 + 2*p2*p2- 4*p2*q2 + 2*p3*p3- 4*p3*q3 + 2*q1*q1+ 2*q2*q2+ 2*q3*q3;
    h_phc_coeffs_Ht[12]=p4*q7 + p7*q4 + p5*q8 + p8*q5 + p6*q9 + p9*q6 - 2*q4*q7 - 2*q5*q8 - 2*q6*q9;
    h_phc_coeffs_Ht[13]=2*p4*p7 + 2*p5*p8 + 2*p6*p9 - 2*p4*q7 - 2*p7*q4 - 2*p5*q8 - 2*p8*q5 - 2*p6*q9 - 2*p9*q6 + 2*q4*q7 + 2*q5*q8 + 2*q6*q9;
    h_phc_coeffs_Ht[14]=2*p4*q4 - 2*q5*q5- 2*q6*q6- 2*q4*q4+ 2*p5*q5 + 2*p6*q6;
    h_phc_coeffs_Ht[15]=2*p4*p4- 4*p4*q4 + 2*p5*p5- 4*p5*q5 + 2*p6*p6- 4*p6*q6 + 2*q4*q4+ 2*q5*q5+ 2*q6*q6;
    h_phc_coeffs_Ht[16]=2*p7*q7 - 2*q8*q8- 2*q9*q9- 2*q7*q7+ 2*p8*q8 + 2*p9*q9;
    h_phc_coeffs_Ht[17]=2*p7*p7- 4*p7*q7 + 2*p8*p8- 4*p8*q8 + 2*p9*p9- 4*p9*q9 + 2*q7*q7+ 2*q8*q8+ 2*q9*q9;
    h_phc_coeffs_Ht[18]=MAGMA_C_ONE;
    h_phc_coeffs_Ht[19]=MAGMA_C_ZERO;

  }
} // end of namespace

#endif
