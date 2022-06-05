#ifndef cpu_eval_HxH_katsura8_h
#define cpu_eval_HxH_katsura8_h
// ============================================================================
// partial derivative evaluations of the katsura8 problem for cpu HC computation
//
// Modifications
//    Chien  21-06-17:   Originally created
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
  void cpu_eval_HxH_katsura8(
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
    magmaFloatComplex G20 = s_track[7] * G13;
    magmaFloatComplex G21 = s_track[1] + s_track[1];
    magmaFloatComplex G22 = G13 * G21;
    magmaFloatComplex G23 = G13 * s_track[0];
    magmaFloatComplex G24 = G23 + G15;
    magmaFloatComplex G25 = G24 + G9;
    magmaFloatComplex G26 = G4 * G21;
    magmaFloatComplex G27 = G26 + G16;
    magmaFloatComplex G28 = G15 + G17;
    magmaFloatComplex G29 = G16 + G18;
    magmaFloatComplex G30 = G17 + G19;
    magmaFloatComplex G31 = G18 + G20;
    magmaFloatComplex G32 = s_track[8] * G13;
    magmaFloatComplex G33 = G19 + G32;
    magmaFloatComplex G34 = s_track[2] + s_track[2];
    magmaFloatComplex G35 = G13 * G34;
    magmaFloatComplex G36 = G13 * s_track[1];
    magmaFloatComplex G37 = G36 + G16;
    magmaFloatComplex G38 = G23 + G17;
    magmaFloatComplex G39 = G38 + G9;
    magmaFloatComplex G40 = G36 + G18;
    magmaFloatComplex G41 = G4 * G34;
    magmaFloatComplex G42 = G41 + G19;
    magmaFloatComplex G43 = G16 + G20;
    magmaFloatComplex G44 = G17 + G32;
    magmaFloatComplex G45 = s_track[3] + s_track[3];
    magmaFloatComplex G46 = G13 * G45;
    magmaFloatComplex G47 = G13 * s_track[2];
    magmaFloatComplex G48 = G47 + G17;
    magmaFloatComplex G49 = G23 + G19;
    magmaFloatComplex G50 = G49 + G9;
    magmaFloatComplex G51 = G36 + G20;
    magmaFloatComplex G52 = G47 + G32;
    magmaFloatComplex G53 = G4 * G45;
    magmaFloatComplex G54 = s_track[4] + s_track[4];
    magmaFloatComplex G55 = G13 * G54;
    magmaFloatComplex G56 = G13 * s_track[3];
    magmaFloatComplex G57 = G56 + G18;
    magmaFloatComplex G58 = G47 + G19;
    magmaFloatComplex G59 = G23 + G32;
    magmaFloatComplex G60 = G59 + G9;
    magmaFloatComplex G61 = s_track[5] + s_track[5];
    magmaFloatComplex G62 = G13 * G61;
    magmaFloatComplex G63 = G13 * s_track[4];
    magmaFloatComplex G64 = G63 + G19;
    magmaFloatComplex G65 = G56 + G20;
    magmaFloatComplex G66 = G23 + G9;
    magmaFloatComplex G67 = s_track[6] + s_track[6];
    magmaFloatComplex G68 = G13 * G67;
    magmaFloatComplex G69 = G13 * s_track[5];
    magmaFloatComplex G70 = G69 + G20;
    magmaFloatComplex G71 = G63 + G32;
    magmaFloatComplex G72 = s_track[7] + s_track[7];
    magmaFloatComplex G73 = G13 * G72;
    magmaFloatComplex G74 = G13 * s_track[6];
    magmaFloatComplex G75 = G74 + G32;
    magmaFloatComplex G76 = s_track[8] + s_track[8];
    magmaFloatComplex G77 = G13 * G76;
    magmaFloatComplex G78 = G13 * s_track[7];
    magmaFloatComplex G79 = s_track[0] * s_track[0];
    magmaFloatComplex G80 = G4 * G79;
    magmaFloatComplex G81 = G9 * s_track[0];
    magmaFloatComplex G82 = G80 + G81;
    magmaFloatComplex G83 = s_track[1] * s_track[1];
    magmaFloatComplex G84 = G13 * G83;
    magmaFloatComplex G85 = G82 + G84;
    magmaFloatComplex G86 = s_track[2] * s_track[2];
    magmaFloatComplex G87 = G13 * G86;
    magmaFloatComplex G88 = G85 + G87;
    magmaFloatComplex G89 = s_track[3] * s_track[3];
    magmaFloatComplex G90 = G13 * G89;
    magmaFloatComplex G91 = G88 + G90;
    magmaFloatComplex G92 = s_track[4] * s_track[4];
    magmaFloatComplex G93 = G13 * G92;
    magmaFloatComplex G94 = G91 + G93;
    magmaFloatComplex G95 = s_track[5] * s_track[5];
    magmaFloatComplex G96 = G13 * G95;
    magmaFloatComplex G97 = G94 + G96;
    magmaFloatComplex G98 = s_track[6] * s_track[6];
    magmaFloatComplex G99 = G13 * G98;
    magmaFloatComplex G100 = G97 + G99;
    magmaFloatComplex G101 = s_track[7] * s_track[7];
    magmaFloatComplex G102 = G13 * G101;
    magmaFloatComplex G103 = G100 + G102;
    magmaFloatComplex G104 = s_track[8] * s_track[8];
    magmaFloatComplex G105 = G13 * G104;
    magmaFloatComplex G106 = G103 + G105;
    magmaFloatComplex G107 = G23 * s_track[1];
    magmaFloatComplex G108 = G36 * s_track[2];
    magmaFloatComplex G109 = G107 + G108;
    magmaFloatComplex G110 = G9 * s_track[1];
    magmaFloatComplex G111 = G109 + G110;
    magmaFloatComplex G112 = G47 * s_track[3];
    magmaFloatComplex G113 = G111 + G112;
    magmaFloatComplex G114 = G56 * s_track[4];
    magmaFloatComplex G115 = G113 + G114;
    magmaFloatComplex G116 = G63 * s_track[5];
    magmaFloatComplex G117 = G115 + G116;
    magmaFloatComplex G118 = G69 * s_track[6];
    magmaFloatComplex G119 = G117 + G118;
    magmaFloatComplex G120 = G74 * s_track[7];
    magmaFloatComplex G121 = G119 + G120;
    magmaFloatComplex G122 = G78 * s_track[8];
    magmaFloatComplex G123 = G121 + G122;
    magmaFloatComplex G124 = G23 * s_track[2];
    magmaFloatComplex G125 = G4 * G83;
    magmaFloatComplex G126 = G124 + G125;
    magmaFloatComplex G127 = G36 * s_track[3];
    magmaFloatComplex G128 = G126 + G127;
    magmaFloatComplex G129 = G47 * s_track[4];
    magmaFloatComplex G130 = G128 + G129;
    magmaFloatComplex G131 = G9 * s_track[2];
    magmaFloatComplex G132 = G130 + G131;
    magmaFloatComplex G133 = G56 * s_track[5];
    magmaFloatComplex G134 = G132 + G133;
    magmaFloatComplex G135 = G63 * s_track[6];
    magmaFloatComplex G136 = G134 + G135;
    magmaFloatComplex G137 = G69 * s_track[7];
    magmaFloatComplex G138 = G136 + G137;
    magmaFloatComplex G139 = G74 * s_track[8];
    magmaFloatComplex G140 = G138 + G139;
    magmaFloatComplex G141 = G23 * s_track[3];
    magmaFloatComplex G142 = G141 + G108;
    magmaFloatComplex G143 = G36 * s_track[4];
    magmaFloatComplex G144 = G142 + G143;
    magmaFloatComplex G145 = G47 * s_track[5];
    magmaFloatComplex G146 = G144 + G145;
    magmaFloatComplex G147 = G56 * s_track[6];
    magmaFloatComplex G148 = G146 + G147;
    magmaFloatComplex G149 = G9 * s_track[3];
    magmaFloatComplex G150 = G148 + G149;
    magmaFloatComplex G151 = G63 * s_track[7];
    magmaFloatComplex G152 = G150 + G151;
    magmaFloatComplex G153 = G69 * s_track[8];
    magmaFloatComplex G154 = G152 + G153;
    magmaFloatComplex G155 = G23 * s_track[4];
    magmaFloatComplex G156 = G155 + G127;
    magmaFloatComplex G157 = G36 * s_track[5];
    magmaFloatComplex G158 = G156 + G157;
    magmaFloatComplex G159 = G4 * G86;
    magmaFloatComplex G160 = G158 + G159;
    magmaFloatComplex G161 = G47 * s_track[6];
    magmaFloatComplex G162 = G160 + G161;
    magmaFloatComplex G163 = G56 * s_track[7];
    magmaFloatComplex G164 = G162 + G163;
    magmaFloatComplex G165 = G63 * s_track[8];
    magmaFloatComplex G166 = G164 + G165;
    magmaFloatComplex G167 = G9 * s_track[4];
    magmaFloatComplex G168 = G166 + G167;
    magmaFloatComplex G169 = G23 * s_track[5];
    magmaFloatComplex G170 = G169 + G143;
    magmaFloatComplex G171 = G36 * s_track[6];
    magmaFloatComplex G172 = G170 + G171;
    magmaFloatComplex G173 = G172 + G112;
    magmaFloatComplex G174 = G47 * s_track[7];
    magmaFloatComplex G175 = G173 + G174;
    magmaFloatComplex G176 = G56 * s_track[8];
    magmaFloatComplex G177 = G175 + G176;
    magmaFloatComplex G178 = G9 * s_track[5];
    magmaFloatComplex G179 = G177 + G178;
    magmaFloatComplex G180 = G23 * s_track[6];
    magmaFloatComplex G181 = G180 + G157;
    magmaFloatComplex G182 = G36 * s_track[7];
    magmaFloatComplex G183 = G181 + G182;
    magmaFloatComplex G184 = G183 + G129;
    magmaFloatComplex G185 = G47 * s_track[8];
    magmaFloatComplex G186 = G184 + G185;
    magmaFloatComplex G187 = G4 * G89;
    magmaFloatComplex G188 = G186 + G187;
    magmaFloatComplex G189 = G9 * s_track[6];
    magmaFloatComplex G190 = G188 + G189;
    magmaFloatComplex G191 = G23 * s_track[7];
    magmaFloatComplex G192 = G191 + G171;
    magmaFloatComplex G193 = G36 * s_track[8];
    magmaFloatComplex G194 = G192 + G193;
    magmaFloatComplex G195 = G194 + G145;
    magmaFloatComplex G196 = G195 + G114;
    magmaFloatComplex G197 = G9 * s_track[7];
    magmaFloatComplex G198 = G196 + G197;
    magmaFloatComplex G199 = G4 * s_track[0];
    magmaFloatComplex G200 = G199 + G36;
    magmaFloatComplex G201 = G200 + G47;
    magmaFloatComplex G202 = G201 + G56;
    magmaFloatComplex G203 = G202 + G63;
    magmaFloatComplex G204 = G203 + G69;
    magmaFloatComplex G205 = G204 + G74;
    magmaFloatComplex G206 = G205 + G78;
    magmaFloatComplex G207 = G13 * s_track[8];
    magmaFloatComplex G208 = G206 + G207;
    magmaFloatComplex G209 = G208 + G9;

    r_cgesvA[0] = G10;
    r_cgesvA[1] = G14;
    r_cgesvA[2] = G15;
    r_cgesvA[3] = G16;
    r_cgesvA[4] = G17;
    r_cgesvA[5] = G18;
    r_cgesvA[6] = G19;
    r_cgesvA[7] = G20;
    r_cgesvA[8] = G4;
    r_cgesvB[0] =G106;

    r_cgesvA[9] = G22;
    r_cgesvA[10] = G25;
    r_cgesvA[11] = G27;
    r_cgesvA[12] = G28;
    r_cgesvA[13] = G29;
    r_cgesvA[14] = G30;
    r_cgesvA[15] = G31;
    r_cgesvA[16] = G33;
    r_cgesvA[17] = G13;
    r_cgesvB[1] =G123;

    r_cgesvA[18] = G35;
    r_cgesvA[19] = G37;
    r_cgesvA[20] = G39;
    r_cgesvA[21] = G40;
    r_cgesvA[22] = G42;
    r_cgesvA[23] = G43;
    r_cgesvA[24] = G44;
    r_cgesvA[25] = G18;
    r_cgesvA[26] = G13;
    r_cgesvB[2] =G140;

    r_cgesvA[27] = G46;
    r_cgesvA[28] = G48;
    r_cgesvA[29] = G40;
    r_cgesvA[30] = G50;
    r_cgesvA[31] = G51;
    r_cgesvA[32] = G52;
    r_cgesvA[33] = G53;
    r_cgesvA[34] = G17;
    r_cgesvA[35] = G13;
    r_cgesvB[3] =G154;

    r_cgesvA[36] = G55;
    r_cgesvA[37] = G57;
    r_cgesvA[38] = G58;
    r_cgesvA[39] = G51;
    r_cgesvA[40] = G60;
    r_cgesvA[41] = G36;
    r_cgesvA[42] = G47;
    r_cgesvA[43] = G56;
    r_cgesvA[44] = G13;
    r_cgesvB[4] =G168;

    r_cgesvA[45] = G62;
    r_cgesvA[46] = G64;
    r_cgesvA[47] = G65;
    r_cgesvA[48] = G52;
    r_cgesvA[49] = G36;
    r_cgesvA[50] = G66;
    r_cgesvA[51] = G36;
    r_cgesvA[52] = G47;
    r_cgesvA[53] = G13;
    r_cgesvB[5] =G179;

    r_cgesvA[54] = G68;
    r_cgesvA[55] = G70;
    r_cgesvA[56] = G71;
    r_cgesvA[57] = G56;
    r_cgesvA[58] = G47;
    r_cgesvA[59] = G36;
    r_cgesvA[60] = G66;
    r_cgesvA[61] = G36;
    r_cgesvA[62] = G13;
    r_cgesvB[6] =G190;

    r_cgesvA[63] = G73;
    r_cgesvA[64] = G75;
    r_cgesvA[65] = G69;
    r_cgesvA[66] = G63;
    r_cgesvA[67] = G56;
    r_cgesvA[68] = G47;
    r_cgesvA[69] = G36;
    r_cgesvA[70] = G66;
    r_cgesvA[71] = G13;
    r_cgesvB[7] =G198;

    r_cgesvA[72] = G77;
    r_cgesvA[73] = G78;
    r_cgesvA[74] = G74;
    r_cgesvA[75] = G69;
    r_cgesvA[76] = G63;
    r_cgesvA[77] = G56;
    r_cgesvA[78] = G47;
    r_cgesvA[79] = G36;
    r_cgesvA[80] = G13;
    r_cgesvB[8] =G209;
  }
}

#endif
