#ifndef cpu_eval_HxHt_katsura8_h
#define cpu_eval_HxHt_katsura8_h
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
  void cpu_eval_HxHt_katsura8(
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
    magmaFloatComplex G81 = s_targetCoefs[0] - s_startCoefs[0];
    magmaFloatComplex G82 = G79 * G81;
    magmaFloatComplex G84 = s_targetCoefs[1] - s_startCoefs[1];
    magmaFloatComplex G85 = s_track[0] * G84;
    magmaFloatComplex G86 = G82 + G85;
    magmaFloatComplex G87 = s_track[1] * s_track[1];
    magmaFloatComplex G89 = s_targetCoefs[2] - s_startCoefs[2];
    magmaFloatComplex G90 = G87 * G89;
    magmaFloatComplex G91 = G86 + G90;
    magmaFloatComplex G92 = s_track[2] * s_track[2];
    magmaFloatComplex G93 = G92 * G89;
    magmaFloatComplex G94 = G91 + G93;
    magmaFloatComplex G95 = s_track[3] * s_track[3];
    magmaFloatComplex G96 = G95 * G89;
    magmaFloatComplex G97 = G94 + G96;
    magmaFloatComplex G98 = s_track[4] * s_track[4];
    magmaFloatComplex G99 = G98 * G89;
    magmaFloatComplex G100 = G97 + G99;
    magmaFloatComplex G101 = s_track[5] * s_track[5];
    magmaFloatComplex G102 = G101 * G89;
    magmaFloatComplex G103 = G100 + G102;
    magmaFloatComplex G104 = s_track[6] * s_track[6];
    magmaFloatComplex G105 = G104 * G89;
    magmaFloatComplex G106 = G103 + G105;
    magmaFloatComplex G107 = s_track[7] * s_track[7];
    magmaFloatComplex G108 = G107 * G89;
    magmaFloatComplex G109 = G106 + G108;
    magmaFloatComplex G110 = s_track[8] * s_track[8];
    magmaFloatComplex G111 = G110 * G89;
    magmaFloatComplex G112 = G109 + G111;
    magmaFloatComplex G113 = s_track[0] * G89;
    magmaFloatComplex G114 = s_track[1] * G113;
    magmaFloatComplex G115 = s_track[1] * G89;
    magmaFloatComplex G116 = s_track[2] * G115;
    magmaFloatComplex G117 = G114 + G116;
    magmaFloatComplex G118 = s_track[1] * G84;
    magmaFloatComplex G119 = G117 + G118;
    magmaFloatComplex G120 = s_track[2] * G89;
    magmaFloatComplex G121 = s_track[3] * G120;
    magmaFloatComplex G122 = G119 + G121;
    magmaFloatComplex G123 = s_track[3] * G89;
    magmaFloatComplex G124 = s_track[4] * G123;
    magmaFloatComplex G125 = G122 + G124;
    magmaFloatComplex G126 = s_track[4] * G89;
    magmaFloatComplex G127 = s_track[5] * G126;
    magmaFloatComplex G128 = G125 + G127;
    magmaFloatComplex G129 = s_track[5] * G89;
    magmaFloatComplex G130 = s_track[6] * G129;
    magmaFloatComplex G131 = G128 + G130;
    magmaFloatComplex G132 = s_track[6] * G89;
    magmaFloatComplex G133 = s_track[7] * G132;
    magmaFloatComplex G134 = G131 + G133;
    magmaFloatComplex G135 = s_track[7] * G89;
    magmaFloatComplex G136 = s_track[8] * G135;
    magmaFloatComplex G137 = G134 + G136;
    magmaFloatComplex G138 = s_track[2] * G113;
    magmaFloatComplex G139 = G87 * G81;
    magmaFloatComplex G140 = G138 + G139;
    magmaFloatComplex G141 = s_track[3] * G115;
    magmaFloatComplex G142 = G140 + G141;
    magmaFloatComplex G143 = s_track[4] * G120;
    magmaFloatComplex G144 = G142 + G143;
    magmaFloatComplex G145 = s_track[2] * G84;
    magmaFloatComplex G146 = G144 + G145;
    magmaFloatComplex G147 = s_track[5] * G123;
    magmaFloatComplex G148 = G146 + G147;
    magmaFloatComplex G149 = s_track[6] * G126;
    magmaFloatComplex G150 = G148 + G149;
    magmaFloatComplex G151 = s_track[7] * G129;
    magmaFloatComplex G152 = G150 + G151;
    magmaFloatComplex G153 = s_track[8] * G132;
    magmaFloatComplex G154 = G152 + G153;
    magmaFloatComplex G155 = s_track[3] * G113;
    magmaFloatComplex G156 = G155 + G116;
    magmaFloatComplex G157 = s_track[4] * G115;
    magmaFloatComplex G158 = G156 + G157;
    magmaFloatComplex G159 = s_track[5] * G120;
    magmaFloatComplex G160 = G158 + G159;
    magmaFloatComplex G161 = s_track[6] * G123;
    magmaFloatComplex G162 = G160 + G161;
    magmaFloatComplex G163 = s_track[3] * G84;
    magmaFloatComplex G164 = G162 + G163;
    magmaFloatComplex G165 = s_track[7] * G126;
    magmaFloatComplex G166 = G164 + G165;
    magmaFloatComplex G167 = s_track[8] * G129;
    magmaFloatComplex G168 = G166 + G167;
    magmaFloatComplex G169 = s_track[4] * G113;
    magmaFloatComplex G170 = G169 + G141;
    magmaFloatComplex G171 = s_track[5] * G115;
    magmaFloatComplex G172 = G170 + G171;
    magmaFloatComplex G173 = G92 * G81;
    magmaFloatComplex G174 = G172 + G173;
    magmaFloatComplex G175 = s_track[6] * G120;
    magmaFloatComplex G176 = G174 + G175;
    magmaFloatComplex G177 = s_track[7] * G123;
    magmaFloatComplex G178 = G176 + G177;
    magmaFloatComplex G179 = s_track[8] * G126;
    magmaFloatComplex G180 = G178 + G179;
    magmaFloatComplex G181 = s_track[4] * G84;
    magmaFloatComplex G182 = G180 + G181;
    magmaFloatComplex G183 = s_track[5] * G113;
    magmaFloatComplex G184 = G183 + G157;
    magmaFloatComplex G185 = s_track[6] * G115;
    magmaFloatComplex G186 = G184 + G185;
    magmaFloatComplex G187 = G186 + G121;
    magmaFloatComplex G188 = s_track[7] * G120;
    magmaFloatComplex G189 = G187 + G188;
    magmaFloatComplex G190 = s_track[8] * G123;
    magmaFloatComplex G191 = G189 + G190;
    magmaFloatComplex G192 = s_track[5] * G84;
    magmaFloatComplex G193 = G191 + G192;
    magmaFloatComplex G194 = s_track[6] * G113;
    magmaFloatComplex G195 = G194 + G171;
    magmaFloatComplex G196 = s_track[7] * G115;
    magmaFloatComplex G197 = G195 + G196;
    magmaFloatComplex G198 = G197 + G143;
    magmaFloatComplex G199 = s_track[8] * G120;
    magmaFloatComplex G200 = G198 + G199;
    magmaFloatComplex G201 = G95 * G81;
    magmaFloatComplex G202 = G200 + G201;
    magmaFloatComplex G203 = s_track[6] * G84;
    magmaFloatComplex G204 = G202 + G203;
    magmaFloatComplex G205 = s_track[7] * G113;
    magmaFloatComplex G206 = G205 + G185;
    magmaFloatComplex G207 = s_track[8] * G115;
    magmaFloatComplex G208 = G206 + G207;
    magmaFloatComplex G209 = G208 + G159;
    magmaFloatComplex G210 = G209 + G124;
    magmaFloatComplex G211 = s_track[7] * G84;
    magmaFloatComplex G212 = G210 + G211;
    magmaFloatComplex G213 = s_track[0] * G81;
    magmaFloatComplex G214 = G213 + G115;
    magmaFloatComplex G215 = G214 + G120;
    magmaFloatComplex G216 = G215 + G123;
    magmaFloatComplex G217 = G216 + G126;
    magmaFloatComplex G218 = G217 + G129;
    magmaFloatComplex G219 = G218 + G132;
    magmaFloatComplex G220 = G219 + G135;
    magmaFloatComplex G221 = s_track[8] * G89;
    magmaFloatComplex G222 = G220 + G221;
    magmaFloatComplex G223 = G222 + G84;

    r_cgesvA[0] = G10;
    r_cgesvA[1] = G14;
    r_cgesvA[2] = G15;
    r_cgesvA[3] = G16;
    r_cgesvA[4] = G17;
    r_cgesvA[5] = G18;
    r_cgesvA[6] = G19;
    r_cgesvA[7] = G20;
    r_cgesvA[8] = G4;
    r_cgesvB[0] = -G112;

    r_cgesvA[9] = G22;
    r_cgesvA[10] = G25;
    r_cgesvA[11] = G27;
    r_cgesvA[12] = G28;
    r_cgesvA[13] = G29;
    r_cgesvA[14] = G30;
    r_cgesvA[15] = G31;
    r_cgesvA[16] = G33;
    r_cgesvA[17] = G13;
    r_cgesvB[1] = -G137;

    r_cgesvA[18] = G35;
    r_cgesvA[19] = G37;
    r_cgesvA[20] = G39;
    r_cgesvA[21] = G40;
    r_cgesvA[22] = G42;
    r_cgesvA[23] = G43;
    r_cgesvA[24] = G44;
    r_cgesvA[25] = G18;
    r_cgesvA[26] = G13;
    r_cgesvB[2] = -G154;

    r_cgesvA[27] = G46;
    r_cgesvA[28] = G48;
    r_cgesvA[29] = G40;
    r_cgesvA[30] = G50;
    r_cgesvA[31] = G51;
    r_cgesvA[32] = G52;
    r_cgesvA[33] = G53;
    r_cgesvA[34] = G17;
    r_cgesvA[35] = G13;
    r_cgesvB[3] = -G168;

    r_cgesvA[36] = G55;
    r_cgesvA[37] = G57;
    r_cgesvA[38] = G58;
    r_cgesvA[39] = G51;
    r_cgesvA[40] = G60;
    r_cgesvA[41] = G36;
    r_cgesvA[42] = G47;
    r_cgesvA[43] = G56;
    r_cgesvA[44] = G13;
    r_cgesvB[4] = -G182;

    r_cgesvA[45] = G62;
    r_cgesvA[46] = G64;
    r_cgesvA[47] = G65;
    r_cgesvA[48] = G52;
    r_cgesvA[49] = G36;
    r_cgesvA[50] = G66;
    r_cgesvA[51] = G36;
    r_cgesvA[52] = G47;
    r_cgesvA[53] = G13;
    r_cgesvB[5] = -G193;

    r_cgesvA[54] = G68;
    r_cgesvA[55] = G70;
    r_cgesvA[56] = G71;
    r_cgesvA[57] = G56;
    r_cgesvA[58] = G47;
    r_cgesvA[59] = G36;
    r_cgesvA[60] = G66;
    r_cgesvA[61] = G36;
    r_cgesvA[62] = G13;
    r_cgesvB[6] = -G204;

    r_cgesvA[63] = G73;
    r_cgesvA[64] = G75;
    r_cgesvA[65] = G69;
    r_cgesvA[66] = G63;
    r_cgesvA[67] = G56;
    r_cgesvA[68] = G47;
    r_cgesvA[69] = G36;
    r_cgesvA[70] = G66;
    r_cgesvA[71] = G13;
    r_cgesvB[7] = -G212;

    r_cgesvA[72] = G77;
    r_cgesvA[73] = G78;
    r_cgesvA[74] = G74;
    r_cgesvA[75] = G69;
    r_cgesvA[76] = G63;
    r_cgesvA[77] = G56;
    r_cgesvA[78] = G47;
    r_cgesvA[79] = G36;
    r_cgesvA[80] = G13;
    r_cgesvB[8] = -G223;

  }
}

#endif
