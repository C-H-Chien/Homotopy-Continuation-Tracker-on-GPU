#ifndef P2C_P3P_10x10_H
#define P2C_P3P_10x10_H
// ==========================================================================================================================
//
// Modifications
//    Chiang-Heng Chien  25-01-19:   Built from P3P p2c-P3P.h.
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

  void p2c_P3P_10x10(
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

	h_phc_coeffs_Hx[0]=q1;
	h_phc_coeffs_Hx[1]=p1 - q1;
	h_phc_coeffs_Hx[2]=q10;
	h_phc_coeffs_Hx[3]=p10 - q10;
	h_phc_coeffs_Hx[4]=q11;
	h_phc_coeffs_Hx[5]=p11 - q11;
	h_phc_coeffs_Hx[6]=q12;
	h_phc_coeffs_Hx[7]=p12 - q12;
	h_phc_coeffs_Hx[8]=q13;
	h_phc_coeffs_Hx[9]=p13 - q13;
	h_phc_coeffs_Hx[10]=q14;
	h_phc_coeffs_Hx[11]=p14 - q14;
	h_phc_coeffs_Hx[12]=q15;
	h_phc_coeffs_Hx[13]=p15 - q15;
	h_phc_coeffs_Hx[14]=q2;
	h_phc_coeffs_Hx[15]=p2 - q2;
	h_phc_coeffs_Hx[16]=q3;
	h_phc_coeffs_Hx[17]=p3 - q3;
	h_phc_coeffs_Hx[18]=q4;
	h_phc_coeffs_Hx[19]=p4 - q4;
	h_phc_coeffs_Hx[20]=q5;
	h_phc_coeffs_Hx[21]=p5 - q5;
	h_phc_coeffs_Hx[22]=q6;
	h_phc_coeffs_Hx[23]=p6 - q6;
	h_phc_coeffs_Hx[24]=q7;
	h_phc_coeffs_Hx[25]=p7 - q7;
	h_phc_coeffs_Hx[26]=q8;
	h_phc_coeffs_Hx[27]=p8 - q8;
	h_phc_coeffs_Hx[28]=q9;
	h_phc_coeffs_Hx[29]=p9 - q9;
	h_phc_coeffs_Hx[30] = MAGMA_C_ONE;
	h_phc_coeffs_Hx[31] = MAGMA_C_ZERO;

	h_phc_coeffs_Ht[0]=p1 - q1;
	h_phc_coeffs_Ht[1]=p10 - q10;
	h_phc_coeffs_Ht[2]=p11 - q11;
	h_phc_coeffs_Ht[3]=p12 - q12;
	h_phc_coeffs_Ht[4]=p13 - q13;
	h_phc_coeffs_Ht[5]=p14 - q14;
	h_phc_coeffs_Ht[6]=p15 - q15;
	h_phc_coeffs_Ht[7]=p2 - q2;
	h_phc_coeffs_Ht[8]=p3 - q3;
	h_phc_coeffs_Ht[9]=p4 - q4;
	h_phc_coeffs_Ht[10]=p5 - q5;
	h_phc_coeffs_Ht[11]=p6 - q6;
	h_phc_coeffs_Ht[12]=p7 - q7;
	h_phc_coeffs_Ht[13]=p8 - q8;
	h_phc_coeffs_Ht[14]=p9 - q9;
	h_phc_coeffs_Ht[15] = MAGMA_C_ONE;

  }
} // end of namespace

#endif
