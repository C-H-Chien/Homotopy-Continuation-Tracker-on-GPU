#ifndef p2c_5pt_rel_pose_w_depth_recon_h
#define p2c_5pt_rel_pose_w_depth_recon_h
// =======================================================================
//
// Modifications
//    Chien  21-12-31:   Initially Created
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

  void p2c_5pt_rel_pose_w_depth_recon(
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

    h_phc_coeffs_Ht[0]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[1]=2*p2 - 2*q2;
    h_phc_coeffs_Ht[2]=2*q2 - 2*p2;
    h_phc_coeffs_Ht[3]=2*p1 - 2*q1;
    h_phc_coeffs_Ht[4]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[5]=q1 - p1;
    h_phc_coeffs_Ht[6]=p3 - q3;
    h_phc_coeffs_Ht[7]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[8]=2*q1 - 2*p1;
    h_phc_coeffs_Ht[9]=q2 - p2;
    h_phc_coeffs_Ht[10]=p4 - q4;
    h_phc_coeffs_Ht[11]= MAGMA_C_ZERO;
    h_phc_coeffs_Ht[12]=2*p6 - 2*q6;
    h_phc_coeffs_Ht[13]=2*q6 - 2*p6;
    h_phc_coeffs_Ht[14]=2*p5 - 2*q5;
    h_phc_coeffs_Ht[15]=q5 - p5;
    h_phc_coeffs_Ht[16]=p7 - q7;
    h_phc_coeffs_Ht[17]=2*q5 - 2*p5;
    h_phc_coeffs_Ht[18]=q6 - p6;
    h_phc_coeffs_Ht[19]=p8 - q8;
    h_phc_coeffs_Ht[20]=2*p10 - 2*q10;
    h_phc_coeffs_Ht[21]=2*q10 - 2*p10;
    h_phc_coeffs_Ht[22]=2*p9 - 2*q9;
    h_phc_coeffs_Ht[23]=q9 - p9;
    h_phc_coeffs_Ht[24]=p11 - q11;
    h_phc_coeffs_Ht[25]=2*q9 - 2*p9;
    h_phc_coeffs_Ht[26]=q10 - p10;
    h_phc_coeffs_Ht[27]=p12 - q12;
    h_phc_coeffs_Ht[28]=2*p14 - 2*q14;
    h_phc_coeffs_Ht[29]=2*q14 - 2*p14;
    h_phc_coeffs_Ht[30]=2*p13 - 2*q13;
    h_phc_coeffs_Ht[31]=q13 - p13;
    h_phc_coeffs_Ht[32]=p15 - q15;
    h_phc_coeffs_Ht[33]=2*q13 - 2*p13;
    h_phc_coeffs_Ht[34]=q14 - p14;
    h_phc_coeffs_Ht[35]=p16 - q16;
    h_phc_coeffs_Ht[36]=2*p18 - 2*q18;
    h_phc_coeffs_Ht[37]=2*q18 - 2*p18;
    h_phc_coeffs_Ht[38]=2*p17 - 2*q17;
    h_phc_coeffs_Ht[39]=q17 - p17;
    h_phc_coeffs_Ht[40]=p19 - q19;
    h_phc_coeffs_Ht[41]=2*q17 - 2*p17;
    h_phc_coeffs_Ht[42]=q18 - p18;
    h_phc_coeffs_Ht[43]=p20 - q20;

    h_phc_coeffs_H[0]=MAGMA_C_MAKE(-2.0, 0.0);
    h_phc_coeffs_H[1]= MAGMA_C_ZERO;
    h_phc_coeffs_H[2]=2*q2;
    h_phc_coeffs_H[3]=2*p2 - 2*q2;
    h_phc_coeffs_H[4]=-2*q2;
    h_phc_coeffs_H[5]=2*q2 - 2*p2;
    h_phc_coeffs_H[6]=2*q1;
    h_phc_coeffs_H[7]=2*p1 - 2*q1;
    h_phc_coeffs_H[8]=MAGMA_C_MAKE(-1.0, 0.0);
    h_phc_coeffs_H[9]= MAGMA_C_ZERO;
    h_phc_coeffs_H[10]=-q1;
    h_phc_coeffs_H[11]=q1 - p1;
    h_phc_coeffs_H[12]=q3;
    h_phc_coeffs_H[13]=p3 - q3;
    h_phc_coeffs_H[14]=MAGMA_C_MAKE(2.0, 0.0);
    h_phc_coeffs_H[15]= MAGMA_C_ZERO;
    h_phc_coeffs_H[16]=-2*q1;
    h_phc_coeffs_H[17]=2*q1 - 2*p1;
    h_phc_coeffs_H[18]=-q2;
    h_phc_coeffs_H[19]=q2 - p2;
    h_phc_coeffs_H[20]=q4;
    h_phc_coeffs_H[21]=p4 - q4;
    h_phc_coeffs_H[22]=MAGMA_C_ONE;
    h_phc_coeffs_H[23]= MAGMA_C_ZERO;
    h_phc_coeffs_H[24]=2*q6;
    h_phc_coeffs_H[25]=2*p6 - 2*q6;
    h_phc_coeffs_H[26]=-2*q6;
    h_phc_coeffs_H[27]=2*q6 - 2*p6;
    h_phc_coeffs_H[28]=2*q5;
    h_phc_coeffs_H[29]=2*p5 - 2*q5;
    h_phc_coeffs_H[30]=-q5;
    h_phc_coeffs_H[31]=q5 - p5;
    h_phc_coeffs_H[32]=q7;
    h_phc_coeffs_H[33]=p7 - q7;
    h_phc_coeffs_H[34]=-2*q5;
    h_phc_coeffs_H[35]=2*q5 - 2*p5;
    h_phc_coeffs_H[36]=-q6;
    h_phc_coeffs_H[37]=q6 - p6;
    h_phc_coeffs_H[38]=q8;
    h_phc_coeffs_H[39]=p8 - q8;
    h_phc_coeffs_H[40]=2*q10;
    h_phc_coeffs_H[41]=2*p10 - 2*q10;
    h_phc_coeffs_H[42]=-2*q10;
    h_phc_coeffs_H[43]=2*q10 - 2*p10;
    h_phc_coeffs_H[44]=2*q9;
    h_phc_coeffs_H[45]=2*p9 - 2*q9;
    h_phc_coeffs_H[46]=-q9;
    h_phc_coeffs_H[47]=q9 - p9;
    h_phc_coeffs_H[48]=q11;
    h_phc_coeffs_H[49]=p11 - q11;
    h_phc_coeffs_H[50]=-2*q9;
    h_phc_coeffs_H[51]=2*q9 - 2*p9;
    h_phc_coeffs_H[52]=-q10;
    h_phc_coeffs_H[53]=q10 - p10;
    h_phc_coeffs_H[54]=q12;
    h_phc_coeffs_H[55]=p12 - q12;
    h_phc_coeffs_H[56]=2*q14;
    h_phc_coeffs_H[57]=2*p14 - 2*q14;
    h_phc_coeffs_H[58]=-2*q14;
    h_phc_coeffs_H[59]=2*q14 - 2*p14;
    h_phc_coeffs_H[60]=2*q13;
    h_phc_coeffs_H[61]=2*p13 - 2*q13;
    h_phc_coeffs_H[62]=-q13;
    h_phc_coeffs_H[63]=q13 - p13;
    h_phc_coeffs_H[64]=q15;
    h_phc_coeffs_H[65]=p15 - q15;
    h_phc_coeffs_H[66]=-2*q13;
    h_phc_coeffs_H[67]=2*q13 - 2*p13;
    h_phc_coeffs_H[68]=-q14;
    h_phc_coeffs_H[69]=q14 - p14;
    h_phc_coeffs_H[70]=q16;
    h_phc_coeffs_H[71]=p16 - q16;
    h_phc_coeffs_H[72]=2*q18;
    h_phc_coeffs_H[73]=2*p18 - 2*q18;
    h_phc_coeffs_H[74]=-2*q18;
    h_phc_coeffs_H[75]=2*q18 - 2*p18;
    h_phc_coeffs_H[76]=2*q17;
    h_phc_coeffs_H[77]=2*p17 - 2*q17;
    h_phc_coeffs_H[78]=-q17;
    h_phc_coeffs_H[79]=q17 - p17;
    h_phc_coeffs_H[80]=q19;
    h_phc_coeffs_H[81]=p19 - q19;
    h_phc_coeffs_H[82]=-2*q17;
    h_phc_coeffs_H[83]=2*q17 - 2*p17;
    h_phc_coeffs_H[84]=-q18;
    h_phc_coeffs_H[85]=q18 - p18;
    h_phc_coeffs_H[86]=q20;
    h_phc_coeffs_H[87]=p20 - q20;
  }
} // end of namespace

#endif
