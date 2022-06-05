#ifndef cpu_eval_HxH_katsura9_h
#define cpu_eval_HxH_katsura9_h
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
  void cpu_eval_HxH_katsura9(
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
    magmaFloatComplex G92 = G4 * G91;
    magmaFloatComplex G93 = G9 * s_track[0];
    magmaFloatComplex G94 = G92 + G93;
    magmaFloatComplex G95 = s_track[1] * s_track[1];
    magmaFloatComplex G96 = G13 * G95;
    magmaFloatComplex G97 = G94 + G96;
    magmaFloatComplex G98 = s_track[2] * s_track[2];
    magmaFloatComplex G99 = G13 * G98;
    magmaFloatComplex G100 = G97 + G99;
    magmaFloatComplex G101 = s_track[3] * s_track[3];
    magmaFloatComplex G102 = G13 * G101;
    magmaFloatComplex G103 = G100 + G102;
    magmaFloatComplex G104 = s_track[4] * s_track[4];
    magmaFloatComplex G105 = G13 * G104;
    magmaFloatComplex G106 = G103 + G105;
    magmaFloatComplex G107 = s_track[5] * s_track[5];
    magmaFloatComplex G108 = G13 * G107;
    magmaFloatComplex G109 = G106 + G108;
    magmaFloatComplex G110 = s_track[6] * s_track[6];
    magmaFloatComplex G111 = G13 * G110;
    magmaFloatComplex G112 = G109 + G111;
    magmaFloatComplex G113 = s_track[7] * s_track[7];
    magmaFloatComplex G114 = G13 * G113;
    magmaFloatComplex G115 = G112 + G114;
    magmaFloatComplex G116 = s_track[8] * s_track[8];
    magmaFloatComplex G117 = G13 * G116;
    magmaFloatComplex G118 = G115 + G117;
    magmaFloatComplex G119 = s_track[9] * s_track[9];
    magmaFloatComplex G120 = G13 * G119;
    magmaFloatComplex G121 = G118 + G120;
    magmaFloatComplex G122 = G24 * s_track[1];
    magmaFloatComplex G123 = G38 * s_track[2];
    magmaFloatComplex G124 = G122 + G123;
    magmaFloatComplex G125 = G9 * s_track[1];
    magmaFloatComplex G126 = G124 + G125;
    magmaFloatComplex G127 = G50 * s_track[3];
    magmaFloatComplex G128 = G126 + G127;
    magmaFloatComplex G129 = G60 * s_track[4];
    magmaFloatComplex G130 = G128 + G129;
    magmaFloatComplex G131 = G69 * s_track[5];
    magmaFloatComplex G132 = G130 + G131;
    magmaFloatComplex G133 = G75 * s_track[6];
    magmaFloatComplex G134 = G132 + G133;
    magmaFloatComplex G135 = G81 * s_track[7];
    magmaFloatComplex G136 = G134 + G135;
    magmaFloatComplex G137 = G86 * s_track[8];
    magmaFloatComplex G138 = G136 + G137;
    magmaFloatComplex G139 = G90 * s_track[9];
    magmaFloatComplex G140 = G138 + G139;
    magmaFloatComplex G141 = G24 * s_track[2];
    magmaFloatComplex G142 = G4 * G95;
    magmaFloatComplex G143 = G141 + G142;
    magmaFloatComplex G144 = G38 * s_track[3];
    magmaFloatComplex G145 = G143 + G144;
    magmaFloatComplex G146 = G50 * s_track[4];
    magmaFloatComplex G147 = G145 + G146;
    magmaFloatComplex G148 = G9 * s_track[2];
    magmaFloatComplex G149 = G147 + G148;
    magmaFloatComplex G150 = G60 * s_track[5];
    magmaFloatComplex G151 = G149 + G150;
    magmaFloatComplex G152 = G69 * s_track[6];
    magmaFloatComplex G153 = G151 + G152;
    magmaFloatComplex G154 = G75 * s_track[7];
    magmaFloatComplex G155 = G153 + G154;
    magmaFloatComplex G156 = G81 * s_track[8];
    magmaFloatComplex G157 = G155 + G156;
    magmaFloatComplex G158 = G86 * s_track[9];
    magmaFloatComplex G159 = G157 + G158;
    magmaFloatComplex G160 = G24 * s_track[3];
    magmaFloatComplex G161 = G160 + G123;
    magmaFloatComplex G162 = G38 * s_track[4];
    magmaFloatComplex G163 = G161 + G162;
    magmaFloatComplex G164 = G50 * s_track[5];
    magmaFloatComplex G165 = G163 + G164;
    magmaFloatComplex G166 = G60 * s_track[6];
    magmaFloatComplex G167 = G165 + G166;
    magmaFloatComplex G168 = G9 * s_track[3];
    magmaFloatComplex G169 = G167 + G168;
    magmaFloatComplex G170 = G69 * s_track[7];
    magmaFloatComplex G171 = G169 + G170;
    magmaFloatComplex G172 = G75 * s_track[8];
    magmaFloatComplex G173 = G171 + G172;
    magmaFloatComplex G174 = G81 * s_track[9];
    magmaFloatComplex G175 = G173 + G174;
    magmaFloatComplex G176 = G24 * s_track[4];
    magmaFloatComplex G177 = G176 + G144;
    magmaFloatComplex G178 = G38 * s_track[5];
    magmaFloatComplex G179 = G177 + G178;
    magmaFloatComplex G180 = G4 * G98;
    magmaFloatComplex G181 = G179 + G180;
    magmaFloatComplex G182 = G50 * s_track[6];
    magmaFloatComplex G183 = G181 + G182;
    magmaFloatComplex G184 = G60 * s_track[7];
    magmaFloatComplex G185 = G183 + G184;
    magmaFloatComplex G186 = G69 * s_track[8];
    magmaFloatComplex G187 = G185 + G186;
    magmaFloatComplex G188 = G9 * s_track[4];
    magmaFloatComplex G189 = G187 + G188;
    magmaFloatComplex G190 = G75 * s_track[9];
    magmaFloatComplex G191 = G189 + G190;
    magmaFloatComplex G192 = G24 * s_track[5];
    magmaFloatComplex G193 = G192 + G162;
    magmaFloatComplex G194 = G38 * s_track[6];
    magmaFloatComplex G195 = G193 + G194;
    magmaFloatComplex G196 = G195 + G127;
    magmaFloatComplex G197 = G50 * s_track[7];
    magmaFloatComplex G198 = G196 + G197;
    magmaFloatComplex G199 = G60 * s_track[8];
    magmaFloatComplex G200 = G198 + G199;
    magmaFloatComplex G201 = G69 * s_track[9];
    magmaFloatComplex G202 = G200 + G201;
    magmaFloatComplex G203 = G9 * s_track[5];
    magmaFloatComplex G204 = G202 + G203;
    magmaFloatComplex G205 = G24 * s_track[6];
    magmaFloatComplex G206 = G205 + G178;
    magmaFloatComplex G207 = G38 * s_track[7];
    magmaFloatComplex G208 = G206 + G207;
    magmaFloatComplex G209 = G208 + G146;
    magmaFloatComplex G210 = G50 * s_track[8];
    magmaFloatComplex G211 = G209 + G210;
    magmaFloatComplex G212 = G4 * G101;
    magmaFloatComplex G213 = G211 + G212;
    magmaFloatComplex G214 = G60 * s_track[9];
    magmaFloatComplex G215 = G213 + G214;
    magmaFloatComplex G216 = G9 * s_track[6];
    magmaFloatComplex G217 = G215 + G216;
    magmaFloatComplex G218 = G24 * s_track[7];
    magmaFloatComplex G219 = G218 + G194;
    magmaFloatComplex G220 = G38 * s_track[8];
    magmaFloatComplex G221 = G219 + G220;
    magmaFloatComplex G222 = G221 + G164;
    magmaFloatComplex G223 = G50 * s_track[9];
    magmaFloatComplex G224 = G222 + G223;
    magmaFloatComplex G225 = G224 + G129;
    magmaFloatComplex G226 = G9 * s_track[7];
    magmaFloatComplex G227 = G225 + G226;
    magmaFloatComplex G228 = G24 * s_track[8];
    magmaFloatComplex G229 = G228 + G207;
    magmaFloatComplex G230 = G38 * s_track[9];
    magmaFloatComplex G231 = G229 + G230;
    magmaFloatComplex G232 = G231 + G182;
    magmaFloatComplex G233 = G232 + G150;
    magmaFloatComplex G234 = G4 * G104;
    magmaFloatComplex G235 = G233 + G234;
    magmaFloatComplex G236 = G9 * s_track[8];
    magmaFloatComplex G237 = G235 + G236;
    magmaFloatComplex G238 = G4 * s_track[0];
    magmaFloatComplex G239 = G238 + G38;
    magmaFloatComplex G240 = G239 + G50;
    magmaFloatComplex G241 = G240 + G60;
    magmaFloatComplex G242 = G241 + G69;
    magmaFloatComplex G243 = G242 + G75;
    magmaFloatComplex G244 = G243 + G81;
    magmaFloatComplex G245 = G244 + G86;
    magmaFloatComplex G246 = G245 + G90;
    magmaFloatComplex G247 = G13 * s_track[9];
    magmaFloatComplex G248 = G246 + G247;
    magmaFloatComplex G249 = G248 + G9;

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
    r_cgesvB[0] = G121;

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
    r_cgesvB[1] = G140;

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
    r_cgesvB[2] = G159;

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
    r_cgesvB[3] = G175;

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
    r_cgesvB[4] = G191;

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
    r_cgesvB[5] = G204;

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
    r_cgesvB[6] = G217;

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
    r_cgesvB[7] = G227;

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
    r_cgesvB[8] = G237;

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
    r_cgesvB[9] = G249;

  }
}

#endif
