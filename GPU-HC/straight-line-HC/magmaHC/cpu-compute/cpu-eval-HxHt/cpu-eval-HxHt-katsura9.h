#ifndef cpu_eval_HxHt_katsura9_h
#define cpu_eval_HxHt_katsura9_h
// ============================================================================
// partial derivative evaluations of the katsura9 problem for cpu HC computation
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
  void cpu_eval_HxHt_katsura9(
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
    magmaFloatComplex G21 = s_track[8] * G13;
    magmaFloatComplex G22 = s_track[1] + s_track[1];
    magmaFloatComplex G23 = G13 * G22;
    magmaFloatComplex G24 = G13 * s_track[0];
    magmaFloatComplex G25 = G24 + G15;
    magmaFloatComplex G26 = G25 + G9;
    magmaFloatComplex G27 = G4 * G22;
    magmaFloatComplex G28 = G27 + G16;
    magmaFloatComplex G29 = G15 + G17;
    magmaFloatComplex G30 = G16 + G18;
    magmaFloatComplex G31 = G17 + G19;
    magmaFloatComplex G32 = G18 + G20;
    magmaFloatComplex G33 = G19 + G21;
    magmaFloatComplex G34 = s_track[9] * G13;
    magmaFloatComplex G35 = G20 + G34;
    magmaFloatComplex G36 = s_track[2] + s_track[2];
    magmaFloatComplex G37 = G13 * G36;
    magmaFloatComplex G38 = G13 * s_track[1];
    magmaFloatComplex G39 = G38 + G16;
    magmaFloatComplex G40 = G24 + G17;
    magmaFloatComplex G41 = G40 + G9;
    magmaFloatComplex G42 = G38 + G18;
    magmaFloatComplex G43 = G4 * G36;
    magmaFloatComplex G44 = G43 + G19;
    magmaFloatComplex G45 = G16 + G20;
    magmaFloatComplex G46 = G17 + G21;
    magmaFloatComplex G47 = G18 + G34;
    magmaFloatComplex G48 = s_track[3] + s_track[3];
    magmaFloatComplex G49 = G13 * G48;
    magmaFloatComplex G50 = G13 * s_track[2];
    magmaFloatComplex G51 = G50 + G17;
    magmaFloatComplex G52 = G24 + G19;
    magmaFloatComplex G53 = G52 + G9;
    magmaFloatComplex G54 = G38 + G20;
    magmaFloatComplex G55 = G50 + G21;
    magmaFloatComplex G56 = G4 * G48;
    magmaFloatComplex G57 = G56 + G34;
    magmaFloatComplex G58 = s_track[4] + s_track[4];
    magmaFloatComplex G59 = G13 * G58;
    magmaFloatComplex G60 = G13 * s_track[3];
    magmaFloatComplex G61 = G60 + G18;
    magmaFloatComplex G62 = G50 + G19;
    magmaFloatComplex G63 = G24 + G21;
    magmaFloatComplex G64 = G63 + G9;
    magmaFloatComplex G65 = G38 + G34;
    magmaFloatComplex G66 = G4 * G58;
    magmaFloatComplex G67 = s_track[5] + s_track[5];
    magmaFloatComplex G68 = G13 * G67;
    magmaFloatComplex G69 = G13 * s_track[4];
    magmaFloatComplex G70 = G69 + G19;
    magmaFloatComplex G71 = G60 + G20;
    magmaFloatComplex G72 = G24 + G9;
    magmaFloatComplex G73 = s_track[6] + s_track[6];
    magmaFloatComplex G74 = G13 * G73;
    magmaFloatComplex G75 = G13 * s_track[5];
    magmaFloatComplex G76 = G75 + G20;
    magmaFloatComplex G77 = G69 + G21;
    magmaFloatComplex G78 = G60 + G34;
    magmaFloatComplex G79 = s_track[7] + s_track[7];
    magmaFloatComplex G80 = G13 * G79;
    magmaFloatComplex G81 = G13 * s_track[6];
    magmaFloatComplex G82 = G81 + G21;
    magmaFloatComplex G83 = G75 + G34;
    magmaFloatComplex G84 = s_track[8] + s_track[8];
    magmaFloatComplex G85 = G13 * G84;
    magmaFloatComplex G86 = G13 * s_track[7];
    magmaFloatComplex G87 = G86 + G34;
    magmaFloatComplex G88 = s_track[9] + s_track[9];
    magmaFloatComplex G89 = G13 * G88;
    magmaFloatComplex G90 = G13 * s_track[8];
    magmaFloatComplex G91 = s_track[0] * s_track[0];
    magmaFloatComplex G93 = s_targetCoefs[0] - s_startCoefs[0];
    magmaFloatComplex G94 = G91 * G93;
    magmaFloatComplex G96 = s_targetCoefs[1] - s_startCoefs[1];
    magmaFloatComplex G97 = s_track[0] * G96;
    magmaFloatComplex G98 = G94 + G97;
    magmaFloatComplex G99 = s_track[1] * s_track[1];
    magmaFloatComplex G101 = s_targetCoefs[2] - s_startCoefs[2];
    magmaFloatComplex G102 = G99 * G101;
    magmaFloatComplex G103 = G98 + G102;
    magmaFloatComplex G104 = s_track[2] * s_track[2];
    magmaFloatComplex G105 = G104 * G101;
    magmaFloatComplex G106 = G103 + G105;
    magmaFloatComplex G107 = s_track[3] * s_track[3];
    magmaFloatComplex G108 = G107 * G101;
    magmaFloatComplex G109 = G106 + G108;
    magmaFloatComplex G110 = s_track[4] * s_track[4];
    magmaFloatComplex G111 = G110 * G101;
    magmaFloatComplex G112 = G109 + G111;
    magmaFloatComplex G113 = s_track[5] * s_track[5];
    magmaFloatComplex G114 = G113 * G101;
    magmaFloatComplex G115 = G112 + G114;
    magmaFloatComplex G116 = s_track[6] * s_track[6];
    magmaFloatComplex G117 = G116 * G101;
    magmaFloatComplex G118 = G115 + G117;
    magmaFloatComplex G119 = s_track[7] * s_track[7];
    magmaFloatComplex G120 = G119 * G101;
    magmaFloatComplex G121 = G118 + G120;
    magmaFloatComplex G122 = s_track[8] * s_track[8];
    magmaFloatComplex G123 = G122 * G101;
    magmaFloatComplex G124 = G121 + G123;
    magmaFloatComplex G125 = s_track[9] * s_track[9];
    magmaFloatComplex G126 = G125 * G101;
    magmaFloatComplex G127 = G124 + G126;
    magmaFloatComplex G128 = s_track[0] * G101;
    magmaFloatComplex G129 = s_track[1] * G128;
    magmaFloatComplex G130 = s_track[1] * G101;
    magmaFloatComplex G131 = s_track[2] * G130;
    magmaFloatComplex G132 = G129 + G131;
    magmaFloatComplex G133 = s_track[1] * G96;
    magmaFloatComplex G134 = G132 + G133;
    magmaFloatComplex G135 = s_track[2] * G101;
    magmaFloatComplex G136 = s_track[3] * G135;
    magmaFloatComplex G137 = G134 + G136;
    magmaFloatComplex G138 = s_track[3] * G101;
    magmaFloatComplex G139 = s_track[4] * G138;
    magmaFloatComplex G140 = G137 + G139;
    magmaFloatComplex G141 = s_track[4] * G101;
    magmaFloatComplex G142 = s_track[5] * G141;
    magmaFloatComplex G143 = G140 + G142;
    magmaFloatComplex G144 = s_track[5] * G101;
    magmaFloatComplex G145 = s_track[6] * G144;
    magmaFloatComplex G146 = G143 + G145;
    magmaFloatComplex G147 = s_track[6] * G101;
    magmaFloatComplex G148 = s_track[7] * G147;
    magmaFloatComplex G149 = G146 + G148;
    magmaFloatComplex G150 = s_track[7] * G101;
    magmaFloatComplex G151 = s_track[8] * G150;
    magmaFloatComplex G152 = G149 + G151;
    magmaFloatComplex G153 = s_track[8] * G101;
    magmaFloatComplex G154 = s_track[9] * G153;
    magmaFloatComplex G155 = G152 + G154;
    magmaFloatComplex G156 = s_track[2] * G128;
    magmaFloatComplex G157 = G99 * G93;
    magmaFloatComplex G158 = G156 + G157;
    magmaFloatComplex G159 = s_track[3] * G130;
    magmaFloatComplex G160 = G158 + G159;
    magmaFloatComplex G161 = s_track[4] * G135;
    magmaFloatComplex G162 = G160 + G161;
    magmaFloatComplex G163 = s_track[2] * G96;
    magmaFloatComplex G164 = G162 + G163;
    magmaFloatComplex G165 = s_track[5] * G138;
    magmaFloatComplex G166 = G164 + G165;
    magmaFloatComplex G167 = s_track[6] * G141;
    magmaFloatComplex G168 = G166 + G167;
    magmaFloatComplex G169 = s_track[7] * G144;
    magmaFloatComplex G170 = G168 + G169;
    magmaFloatComplex G171 = s_track[8] * G147;
    magmaFloatComplex G172 = G170 + G171;
    magmaFloatComplex G173 = s_track[9] * G150;
    magmaFloatComplex G174 = G172 + G173;
    magmaFloatComplex G175 = s_track[3] * G128;
    magmaFloatComplex G176 = G175 + G131;
    magmaFloatComplex G177 = s_track[4] * G130;
    magmaFloatComplex G178 = G176 + G177;
    magmaFloatComplex G179 = s_track[5] * G135;
    magmaFloatComplex G180 = G178 + G179;
    magmaFloatComplex G181 = s_track[6] * G138;
    magmaFloatComplex G182 = G180 + G181;
    magmaFloatComplex G183 = s_track[3] * G96;
    magmaFloatComplex G184 = G182 + G183;
    magmaFloatComplex G185 = s_track[7] * G141;
    magmaFloatComplex G186 = G184 + G185;
    magmaFloatComplex G187 = s_track[8] * G144;
    magmaFloatComplex G188 = G186 + G187;
    magmaFloatComplex G189 = s_track[9] * G147;
    magmaFloatComplex G190 = G188 + G189;
    magmaFloatComplex G191 = s_track[4] * G128;
    magmaFloatComplex G192 = G191 + G159;
    magmaFloatComplex G193 = s_track[5] * G130;
    magmaFloatComplex G194 = G192 + G193;
    magmaFloatComplex G195 = G104 * G93;
    magmaFloatComplex G196 = G194 + G195;
    magmaFloatComplex G197 = s_track[6] * G135;
    magmaFloatComplex G198 = G196 + G197;
    magmaFloatComplex G199 = s_track[7] * G138;
    magmaFloatComplex G200 = G198 + G199;
    magmaFloatComplex G201 = s_track[8] * G141;
    magmaFloatComplex G202 = G200 + G201;
    magmaFloatComplex G203 = s_track[4] * G96;
    magmaFloatComplex G204 = G202 + G203;
    magmaFloatComplex G205 = s_track[9] * G144;
    magmaFloatComplex G206 = G204 + G205;
    magmaFloatComplex G207 = s_track[5] * G128;
    magmaFloatComplex G208 = G207 + G177;
    magmaFloatComplex G209 = s_track[6] * G130;
    magmaFloatComplex G210 = G208 + G209;
    magmaFloatComplex G211 = G210 + G136;
    magmaFloatComplex G212 = s_track[7] * G135;
    magmaFloatComplex G213 = G211 + G212;
    magmaFloatComplex G214 = s_track[8] * G138;
    magmaFloatComplex G215 = G213 + G214;
    magmaFloatComplex G216 = s_track[9] * G141;
    magmaFloatComplex G217 = G215 + G216;
    magmaFloatComplex G218 = s_track[5] * G96;
    magmaFloatComplex G219 = G217 + G218;
    magmaFloatComplex G220 = s_track[6] * G128;
    magmaFloatComplex G221 = G220 + G193;
    magmaFloatComplex G222 = s_track[7] * G130;
    magmaFloatComplex G223 = G221 + G222;
    magmaFloatComplex G224 = G223 + G161;
    magmaFloatComplex G225 = s_track[8] * G135;
    magmaFloatComplex G226 = G224 + G225;
    magmaFloatComplex G227 = G107 * G93;
    magmaFloatComplex G228 = G226 + G227;
    magmaFloatComplex G229 = s_track[9] * G138;
    magmaFloatComplex G230 = G228 + G229;
    magmaFloatComplex G231 = s_track[6] * G96;
    magmaFloatComplex G232 = G230 + G231;
    magmaFloatComplex G233 = s_track[7] * G128;
    magmaFloatComplex G234 = G233 + G209;
    magmaFloatComplex G235 = s_track[8] * G130;
    magmaFloatComplex G236 = G234 + G235;
    magmaFloatComplex G237 = G236 + G179;
    magmaFloatComplex G238 = s_track[9] * G135;
    magmaFloatComplex G239 = G237 + G238;
    magmaFloatComplex G240 = G239 + G139;
    magmaFloatComplex G241 = s_track[7] * G96;
    magmaFloatComplex G242 = G240 + G241;
    magmaFloatComplex G243 = s_track[8] * G128;
    magmaFloatComplex G244 = G243 + G222;
    magmaFloatComplex G245 = s_track[9] * G130;
    magmaFloatComplex G246 = G244 + G245;
    magmaFloatComplex G247 = G246 + G197;
    magmaFloatComplex G248 = G247 + G165;
    magmaFloatComplex G249 = G110 * G93;
    magmaFloatComplex G250 = G248 + G249;
    magmaFloatComplex G251 = s_track[8] * G96;
    magmaFloatComplex G252 = G250 + G251;
    magmaFloatComplex G253 = s_track[0] * G93;
    magmaFloatComplex G254 = G253 + G130;
    magmaFloatComplex G255 = G254 + G135;
    magmaFloatComplex G256 = G255 + G138;
    magmaFloatComplex G257 = G256 + G141;
    magmaFloatComplex G258 = G257 + G144;
    magmaFloatComplex G259 = G258 + G147;
    magmaFloatComplex G260 = G259 + G150;
    magmaFloatComplex G261 = G260 + G153;
    magmaFloatComplex G262 = s_track[9] * G101;
    magmaFloatComplex G263 = G261 + G262;
    magmaFloatComplex G264 = G263 + G96;

    r_cgesvA[0] = G10;
    r_cgesvA[1] = G14;
    r_cgesvA[2] = G15;
    r_cgesvA[3] = G16;
    r_cgesvA[4] = G17;
    r_cgesvA[5] = G18;
    r_cgesvA[6] = G19;
    r_cgesvA[7] = G20;
    r_cgesvA[8] = G21;
    r_cgesvA[9] = G4;
    r_cgesvB[0] = -G127;

    r_cgesvA[10] = G23;
    r_cgesvA[11] = G26;
    r_cgesvA[12] = G28;
    r_cgesvA[13] = G29;
    r_cgesvA[14] = G30;
    r_cgesvA[15] = G31;
    r_cgesvA[16] = G32;
    r_cgesvA[17] = G33;
    r_cgesvA[18] = G35;
    r_cgesvA[19] = G13;
    r_cgesvB[1] = -G155;

    r_cgesvA[20] = G37;
    r_cgesvA[21] = G39;
    r_cgesvA[22] = G41;
    r_cgesvA[23] = G42;
    r_cgesvA[24] = G44;
    r_cgesvA[25] = G45;
    r_cgesvA[26] = G46;
    r_cgesvA[27] = G47;
    r_cgesvA[28] = G19;
    r_cgesvA[29] = G13;
    r_cgesvB[2] = -G174;

    r_cgesvA[30] = G49;
    r_cgesvA[31] = G51;
    r_cgesvA[32] = G42;
    r_cgesvA[33] = G53;
    r_cgesvA[34] = G54;
    r_cgesvA[35] = G55;
    r_cgesvA[36] = G57;
    r_cgesvA[37] = G17;
    r_cgesvA[38] = G18;
    r_cgesvA[39] = G13;
    r_cgesvB[3] = -G190;

    r_cgesvA[40] = G59;
    r_cgesvA[41] = G61;
    r_cgesvA[42] = G62;
    r_cgesvA[43] = G54;
    r_cgesvA[44] = G64;
    r_cgesvA[45] = G65;
    r_cgesvA[46] = G50;
    r_cgesvA[47] = G60;
    r_cgesvA[48] = G66;
    r_cgesvA[49] = G13;
    r_cgesvB[4] = -G206;

    r_cgesvA[50] = G68;
    r_cgesvA[51] = G70;
    r_cgesvA[52] = G71;
    r_cgesvA[53] = G55;
    r_cgesvA[54] = G65;
    r_cgesvA[55] = G72;
    r_cgesvA[56] = G38;
    r_cgesvA[57] = G50;
    r_cgesvA[58] = G60;
    r_cgesvA[59] = G13;
    r_cgesvB[5] = -G219;

    r_cgesvA[60] = G74;
    r_cgesvA[61] = G76;
    r_cgesvA[62] = G77;
    r_cgesvA[63] = G78;
    r_cgesvA[64] = G50;
    r_cgesvA[65] = G38;
    r_cgesvA[66] = G72;
    r_cgesvA[67] = G38;
    r_cgesvA[68] = G50;
    r_cgesvA[69] = G13;
    r_cgesvB[6] = -G232;

    r_cgesvA[70] = G80;
    r_cgesvA[71] = G82;
    r_cgesvA[72] = G83;
    r_cgesvA[73] = G69;
    r_cgesvA[74] = G60;
    r_cgesvA[75] = G50;
    r_cgesvA[76] = G38;
    r_cgesvA[77] = G72;
    r_cgesvA[78] = G38;
    r_cgesvA[79] = G13;
    r_cgesvB[7] = -G242;

    r_cgesvA[80] = G85;
    r_cgesvA[81] = G87;
    r_cgesvA[82] = G81;
    r_cgesvA[83] = G75;
    r_cgesvA[84] = G69;
    r_cgesvA[85] = G60;
    r_cgesvA[86] = G50;
    r_cgesvA[87] = G38;
    r_cgesvA[88] = G72;
    r_cgesvA[89] = G13;
    r_cgesvB[8] = -G252;

    r_cgesvA[90] = G89;
    r_cgesvA[91] = G90;
    r_cgesvA[92] = G86;
    r_cgesvA[93] = G81;
    r_cgesvA[94] = G75;
    r_cgesvA[95] = G69;
    r_cgesvA[96] = G60;
    r_cgesvA[97] = G50;
    r_cgesvA[98] = G38;
    r_cgesvA[99] = G13;
    r_cgesvB[9] = -G264;
  }
}

#endif
