#ifndef cpu_eval_HxH_katsura7_h
#define cpu_eval_HxH_katsura7_h
// ============================================================================
// partial derivative evaluations of the katsura7 problem for cpu HC computation
//
// Modifications
//    Chien  21-12-17:   Originally created
//
// ============================================================================
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <chrono>

// -- cuda --
#include <cuda.h>
#include <cuda_runtime.h>

// -- magma --
#include "flops.h"
#include "magma_v2.h"
#include "magma_lapack.h"
#include "magma_internal.h"
#undef max
#undef min

namespace magmaHCWrapper {

  extern "C"
  void cpu_eval_HxH_katsura7(
      magma_int_t s, float t, int N, magmaFloatComplex* s_track,
      const magmaFloatComplex &C0, const magmaFloatComplex &C1,
      magmaFloatComplex* s_startCoefs, magmaFloatComplex* s_targetCoefs,
      magmaFloatComplex* r_cgesvA, magmaFloatComplex* r_cgesvB )
  {
    magmaFloatComplex G1 = C0 - t;
    magmaFloatComplex G2 = G1 * s_startCoefs[0];
    magmaFloatComplex G3 = t * s_targetCoefs[0];
    magmaFloatComplex G4 = G2 + G3;
    magmaFloatComplex G5 = s_track[0] + s_track[0];
    magmaFloatComplex G6 = G4 * G5;
    magmaFloatComplex G7 = G1 * s_startCoefs[1];
    magmaFloatComplex G8 = t * s_targetCoefs[1];
    magmaFloatComplex G9 = G7 + G8;
    magmaFloatComplex G10 = G6 + G9;
    magmaFloatComplex G11 = G1 * s_startCoefs[2];
    magmaFloatComplex G12 = t * s_targetCoefs[2];
    magmaFloatComplex G13 = G11 + G12;
    magmaFloatComplex G14 = s_track[1] * G13;
    magmaFloatComplex G15 = s_track[2] * G13;
    magmaFloatComplex G16 = s_track[3] * G13;
    magmaFloatComplex G17 = s_track[4] * G13;
    magmaFloatComplex G18 = s_track[5] * G13;
    magmaFloatComplex G19 = s_track[6] * G13;
    magmaFloatComplex G20 = s_track[1] + s_track[1];
    magmaFloatComplex G21 = G13 * G20;
    magmaFloatComplex G22 = G13 * s_track[0];
    magmaFloatComplex G23 = G22 + G15;
    magmaFloatComplex G24 = G23 + G9;
    magmaFloatComplex G25 = G4 * G20;
    magmaFloatComplex G26 = G25 + G16;
    magmaFloatComplex G27 = G15 + G17;
    magmaFloatComplex G28 = G16 + G18;
    magmaFloatComplex G29 = G17 + G19;
    magmaFloatComplex G30 = s_track[7] * G13;
    magmaFloatComplex G31 = G18 + G30;
    magmaFloatComplex G32 = s_track[2] + s_track[2];
    magmaFloatComplex G33 = G13 * G32;
    magmaFloatComplex G34 = G13 * s_track[1];
    magmaFloatComplex G35 = G34 + G16;
    magmaFloatComplex G36 = G22 + G17;
    magmaFloatComplex G37 = G36 + G9;
    magmaFloatComplex G38 = G34 + G18;
    magmaFloatComplex G39 = G4 * G32;
    magmaFloatComplex G40 = G39 + G19;
    magmaFloatComplex G41 = G16 + G30;
    magmaFloatComplex G42 = s_track[3] + s_track[3];
    magmaFloatComplex G43 = G13 * G42;
    magmaFloatComplex G44 = G13 * s_track[2];
    magmaFloatComplex G45 = G44 + G17;
    magmaFloatComplex G46 = G22 + G19;
    magmaFloatComplex G47 = G46 + G9;
    magmaFloatComplex G48 = G34 + G30;
    magmaFloatComplex G49 = G4 * G42;
    magmaFloatComplex G50 = s_track[4] + s_track[4];
    magmaFloatComplex G51 = G13 * G50;
    magmaFloatComplex G52 = G13 * s_track[3];
    magmaFloatComplex G53 = G52 + G18;
    magmaFloatComplex G54 = G44 + G19;
    magmaFloatComplex G55 = G22 + G9;
    magmaFloatComplex G56 = s_track[5] + s_track[5];
    magmaFloatComplex G57 = G13 * G56;
    magmaFloatComplex G58 = G13 * s_track[4];
    magmaFloatComplex G59 = G58 + G19;
    magmaFloatComplex G60 = G52 + G30;
    magmaFloatComplex G61 = s_track[6] + s_track[6];
    magmaFloatComplex G62 = G13 * G61;
    magmaFloatComplex G63 = G13 * s_track[5];
    magmaFloatComplex G64 = G63 + G30;
    magmaFloatComplex G65 = s_track[7] + s_track[7];
    magmaFloatComplex G66 = G13 * G65;
    magmaFloatComplex G67 = G13 * s_track[6];
    magmaFloatComplex G68 = s_track[0] * s_track[0];
    magmaFloatComplex G69 = G4 * G68;
    magmaFloatComplex G70 = G9 * s_track[0];
    magmaFloatComplex G71 = G69 + G70;
    magmaFloatComplex G72 = s_track[1] * s_track[1];
    magmaFloatComplex G73 = G13 * G72;
    magmaFloatComplex G74 = G71 + G73;
    magmaFloatComplex G75 = s_track[2] * s_track[2];
    magmaFloatComplex G76 = G13 * G75;
    magmaFloatComplex G77 = G74 + G76;
    magmaFloatComplex G78 = s_track[3] * s_track[3];
    magmaFloatComplex G79 = G13 * G78;
    magmaFloatComplex G80 = G77 + G79;
    magmaFloatComplex G81 = s_track[4] * s_track[4];
    magmaFloatComplex G82 = G13 * G81;
    magmaFloatComplex G83 = G80 + G82;
    magmaFloatComplex G84 = s_track[5] * s_track[5];
    magmaFloatComplex G85 = G13 * G84;
    magmaFloatComplex G86 = G83 + G85;
    magmaFloatComplex G87 = s_track[6] * s_track[6];
    magmaFloatComplex G88 = G13 * G87;
    magmaFloatComplex G89 = G86 + G88;
    magmaFloatComplex G90 = s_track[7] * s_track[7];
    magmaFloatComplex G91 = G13 * G90;
    magmaFloatComplex G92 = G89 + G91;
    magmaFloatComplex G93 = G22 * s_track[1];
    magmaFloatComplex G94 = G34 * s_track[2];
    magmaFloatComplex G95 = G93 + G94;
    magmaFloatComplex G96 = G9 * s_track[1];
    magmaFloatComplex G97 = G95 + G96;
    magmaFloatComplex G98 = G44 * s_track[3];
    magmaFloatComplex G99 = G97 + G98;
    magmaFloatComplex G100 = G52 * s_track[4];
    magmaFloatComplex G101 = G99 + G100;
    magmaFloatComplex G102 = G58 * s_track[5];
    magmaFloatComplex G103 = G101 + G102;
    magmaFloatComplex G104 = G63 * s_track[6];
    magmaFloatComplex G105 = G103 + G104;
    magmaFloatComplex G106 = G67 * s_track[7];
    magmaFloatComplex G107 = G105 + G106;
    magmaFloatComplex G108 = G22 * s_track[2];
    magmaFloatComplex G109 = G4 * G72;
    magmaFloatComplex G110 = G108 + G109;
    magmaFloatComplex G111 = G34 * s_track[3];
    magmaFloatComplex G112 = G110 + G111;
    magmaFloatComplex G113 = G44 * s_track[4];
    magmaFloatComplex G114 = G112 + G113;
    magmaFloatComplex G115 = G9 * s_track[2];
    magmaFloatComplex G116 = G114 + G115;
    magmaFloatComplex G117 = G52 * s_track[5];
    magmaFloatComplex G118 = G116 + G117;
    magmaFloatComplex G119 = G58 * s_track[6];
    magmaFloatComplex G120 = G118 + G119;
    magmaFloatComplex G121 = G63 * s_track[7];
    magmaFloatComplex G122 = G120 + G121;
    magmaFloatComplex G123 = G22 * s_track[3];
    magmaFloatComplex G124 = G123 + G94;
    magmaFloatComplex G125 = G34 * s_track[4];
    magmaFloatComplex G126 = G124 + G125;
    magmaFloatComplex G127 = G44 * s_track[5];
    magmaFloatComplex G128 = G126 + G127;
    magmaFloatComplex G129 = G52 * s_track[6];
    magmaFloatComplex G130 = G128 + G129;
    magmaFloatComplex G131 = G9 * s_track[3];
    magmaFloatComplex G132 = G130 + G131;
    magmaFloatComplex G133 = G58 * s_track[7];
    magmaFloatComplex G134 = G132 + G133;
    magmaFloatComplex G135 = G22 * s_track[4];
    magmaFloatComplex G136 = G135 + G111;
    magmaFloatComplex G137 = G34 * s_track[5];
    magmaFloatComplex G138 = G136 + G137;
    magmaFloatComplex G139 = G4 * G75;
    magmaFloatComplex G140 = G138 + G139;
    magmaFloatComplex G141 = G44 * s_track[6];
    magmaFloatComplex G142 = G140 + G141;
    magmaFloatComplex G143 = G52 * s_track[7];
    magmaFloatComplex G144 = G142 + G143;
    magmaFloatComplex G145 = G9 * s_track[4];
    magmaFloatComplex G146 = G144 + G145;
    magmaFloatComplex G147 = G22 * s_track[5];
    magmaFloatComplex G148 = G147 + G125;
    magmaFloatComplex G149 = G34 * s_track[6];
    magmaFloatComplex G150 = G148 + G149;
    magmaFloatComplex G151 = G150 + G98;
    magmaFloatComplex G152 = G44 * s_track[7];
    magmaFloatComplex G153 = G151 + G152;
    magmaFloatComplex G154 = G9 * s_track[5];
    magmaFloatComplex G155 = G153 + G154;
    magmaFloatComplex G156 = G22 * s_track[6];
    magmaFloatComplex G157 = G156 + G137;
    magmaFloatComplex G158 = G34 * s_track[7];
    magmaFloatComplex G159 = G157 + G158;
    magmaFloatComplex G160 = G159 + G113;
    magmaFloatComplex G161 = G4 * G78;
    magmaFloatComplex G162 = G160 + G161;
    magmaFloatComplex G163 = G9 * s_track[6];
    magmaFloatComplex G164 = G162 + G163;
    magmaFloatComplex G165 = G4 * s_track[0];
    magmaFloatComplex G166 = G165 + G34;
    magmaFloatComplex G167 = G166 + G44;
    magmaFloatComplex G168 = G167 + G52;
    magmaFloatComplex G169 = G168 + G58;
    magmaFloatComplex G170 = G169 + G63;
    magmaFloatComplex G171 = G170 + G67;
    magmaFloatComplex G172 = G13 * s_track[7];
    magmaFloatComplex G173 = G171 + G172;
    magmaFloatComplex G174 = G173 + G9;

    r_cgesvA[0] = G10;
    r_cgesvA[1] = G14;
    r_cgesvA[2] = G15;
    r_cgesvA[3] = G16;
    r_cgesvA[4] = G17;
    r_cgesvA[5] = G18;
    r_cgesvA[6] = G19;
    r_cgesvA[7] = G4;
    r_cgesvB[0] =G92;

    r_cgesvA[8] = G21;
    r_cgesvA[9] = G24;
    r_cgesvA[10] = G26;
    r_cgesvA[11] = G27;
    r_cgesvA[12] = G28;
    r_cgesvA[13] = G29;
    r_cgesvA[14] = G31;
    r_cgesvA[15] = G13;
    r_cgesvB[1] =G107;

    r_cgesvA[16] = G33;
    r_cgesvA[17] = G35;
    r_cgesvA[18] = G37;
    r_cgesvA[19] = G38;
    r_cgesvA[20] = G40;
    r_cgesvA[21] = G41;
    r_cgesvA[22] = G17;
    r_cgesvA[23] = G13;
    r_cgesvB[2] =G122;

    r_cgesvA[24] = G43;
    r_cgesvA[25] = G45;
    r_cgesvA[26] = G38;
    r_cgesvA[27] = G47;
    r_cgesvA[28] = G48;
    r_cgesvA[29] = G44;
    r_cgesvA[30] = G49;
    r_cgesvA[31] = G13;
    r_cgesvB[3] =G134;

    r_cgesvA[32] = G51;
    r_cgesvA[33] = G53;
    r_cgesvA[34] = G54;
    r_cgesvA[35] = G48;
    r_cgesvA[36] = G55;
    r_cgesvA[37] = G34;
    r_cgesvA[38] = G44;
    r_cgesvA[39] = G13;
    r_cgesvB[4] =G146;

    r_cgesvA[40] = G57;
    r_cgesvA[41] = G59;
    r_cgesvA[42] = G60;
    r_cgesvA[43] = G44;
    r_cgesvA[44] = G34;
    r_cgesvA[45] = G55;
    r_cgesvA[46] = G34;
    r_cgesvA[47] = G13;
    r_cgesvB[5] =G155;

    r_cgesvA[48] = G62;
    r_cgesvA[49] = G64;
    r_cgesvA[50] = G58;
    r_cgesvA[51] = G52;
    r_cgesvA[52] = G44;
    r_cgesvA[53] = G34;
    r_cgesvA[54] = G55;
    r_cgesvA[55] = G13;
    r_cgesvB[6] =G164;

    r_cgesvA[56] = G66;
    r_cgesvA[57] = G67;
    r_cgesvA[58] = G63;
    r_cgesvA[59] = G58;
    r_cgesvA[60] = G52;
    r_cgesvA[61] = G44;
    r_cgesvA[62] = G34;
    r_cgesvA[63] = G13;
    r_cgesvB[7] =G174;
  }
}

#endif
