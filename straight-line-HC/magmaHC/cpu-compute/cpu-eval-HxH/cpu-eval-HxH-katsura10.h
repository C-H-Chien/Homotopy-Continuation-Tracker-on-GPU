#ifndef cpu_eval_HxH_katsura10_h
#define cpu_eval_HxH_katsura10_h
// ============================================================================
// partial derivative evaluations of the katsura10 problem for cpu HC computation
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
  void cpu_eval_HxH_katsura10(
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
    magmaFloatComplex G22 = s_track[9] * G13;
    magmaFloatComplex G23 = s_track[1] + s_track[1];
    magmaFloatComplex G24 = G13 * G23;
    magmaFloatComplex G25 = G13 * s_track[0];
    magmaFloatComplex G26 = G25 + G15;
    magmaFloatComplex G27 = G26 + G9;
    magmaFloatComplex G28 = G4 * G23;
    magmaFloatComplex G29 = G28 + G16;
    magmaFloatComplex G30 = G15 + G17;
    magmaFloatComplex G31 = G16 + G18;
    magmaFloatComplex G32 = G17 + G19;
    magmaFloatComplex G33 = G18 + G20;
    magmaFloatComplex G34 = G19 + G21;
    magmaFloatComplex G35 = G20 + G22;
    magmaFloatComplex G36 = s_track[10] * G13;
    magmaFloatComplex G37 = G21 + G36;
    magmaFloatComplex G38 = s_track[2] + s_track[2];
    magmaFloatComplex G39 = G13 * G38;
    magmaFloatComplex G40 = G13 * s_track[1];
    magmaFloatComplex G41 = G40 + G16;
    magmaFloatComplex G42 = G25 + G17;
    magmaFloatComplex G43 = G42 + G9;
    magmaFloatComplex G44 = G40 + G18;
    magmaFloatComplex G45 = G4 * G38;
    magmaFloatComplex G46 = G45 + G19;
    magmaFloatComplex G47 = G16 + G20;
    magmaFloatComplex G48 = G17 + G21;
    magmaFloatComplex G49 = G18 + G22;
    magmaFloatComplex G50 = G19 + G36;
    magmaFloatComplex G51 = s_track[3] + s_track[3];
    magmaFloatComplex G52 = G13 * G51;
    magmaFloatComplex G53 = G13 * s_track[2];
    magmaFloatComplex G54 = G53 + G17;
    magmaFloatComplex G55 = G25 + G19;
    magmaFloatComplex G56 = G55 + G9;
    magmaFloatComplex G57 = G40 + G20;
    magmaFloatComplex G58 = G53 + G21;
    magmaFloatComplex G59 = G4 * G51;
    magmaFloatComplex G60 = G59 + G22;
    magmaFloatComplex G61 = G17 + G36;
    magmaFloatComplex G62 = s_track[4] + s_track[4];
    magmaFloatComplex G63 = G13 * G62;
    magmaFloatComplex G64 = G13 * s_track[3];
    magmaFloatComplex G65 = G64 + G18;
    magmaFloatComplex G66 = G53 + G19;
    magmaFloatComplex G67 = G25 + G21;
    magmaFloatComplex G68 = G67 + G9;
    magmaFloatComplex G69 = G40 + G22;
    magmaFloatComplex G70 = G53 + G36;
    magmaFloatComplex G71 = G4 * G62;
    magmaFloatComplex G72 = s_track[5] + s_track[5];
    magmaFloatComplex G73 = G13 * G72;
    magmaFloatComplex G74 = G13 * s_track[4];
    magmaFloatComplex G75 = G74 + G19;
    magmaFloatComplex G76 = G64 + G20;
    magmaFloatComplex G77 = G25 + G36;
    magmaFloatComplex G78 = G77 + G9;
    magmaFloatComplex G79 = s_track[6] + s_track[6];
    magmaFloatComplex G80 = G13 * G79;
    magmaFloatComplex G81 = G13 * s_track[5];
    magmaFloatComplex G82 = G81 + G20;
    magmaFloatComplex G83 = G74 + G21;
    magmaFloatComplex G84 = G64 + G22;
    magmaFloatComplex G85 = G25 + G9;
    magmaFloatComplex G86 = s_track[7] + s_track[7];
    magmaFloatComplex G87 = G13 * G86;
    magmaFloatComplex G88 = G13 * s_track[6];
    magmaFloatComplex G89 = G88 + G21;
    magmaFloatComplex G90 = G81 + G22;
    magmaFloatComplex G91 = G74 + G36;
    magmaFloatComplex G92 = s_track[8] + s_track[8];
    magmaFloatComplex G93 = G13 * G92;
    magmaFloatComplex G94 = G13 * s_track[7];
    magmaFloatComplex G95 = G94 + G22;
    magmaFloatComplex G96 = G88 + G36;
    magmaFloatComplex G97 = s_track[9] + s_track[9];
    magmaFloatComplex G98 = G13 * G97;
    magmaFloatComplex G99 = G13 * s_track[8];
    magmaFloatComplex G100 = G99 + G36;
    magmaFloatComplex G101 = s_track[10] + s_track[10];
    magmaFloatComplex G102 = G13 * G101;
    magmaFloatComplex G103 = G13 * s_track[9];
    magmaFloatComplex G104 = s_track[0] * s_track[0];
    magmaFloatComplex G105 = G4 * G104;
    magmaFloatComplex G106 = G9 * s_track[0];
    magmaFloatComplex G107 = G105 + G106;
    magmaFloatComplex G108 = s_track[1] * s_track[1];
    magmaFloatComplex G109 = G13 * G108;
    magmaFloatComplex G110 = G107 + G109;
    magmaFloatComplex G111 = s_track[2] * s_track[2];
    magmaFloatComplex G112 = G13 * G111;
    magmaFloatComplex G113 = G110 + G112;
    magmaFloatComplex G114 = s_track[3] * s_track[3];
    magmaFloatComplex G115 = G13 * G114;
    magmaFloatComplex G116 = G113 + G115;
    magmaFloatComplex G117 = s_track[4] * s_track[4];
    magmaFloatComplex G118 = G13 * G117;
    magmaFloatComplex G119 = G116 + G118;
    magmaFloatComplex G120 = s_track[5] * s_track[5];
    magmaFloatComplex G121 = G13 * G120;
    magmaFloatComplex G122 = G119 + G121;
    magmaFloatComplex G123 = s_track[6] * s_track[6];
    magmaFloatComplex G124 = G13 * G123;
    magmaFloatComplex G125 = G122 + G124;
    magmaFloatComplex G126 = s_track[7] * s_track[7];
    magmaFloatComplex G127 = G13 * G126;
    magmaFloatComplex G128 = G125 + G127;
    magmaFloatComplex G129 = s_track[8] * s_track[8];
    magmaFloatComplex G130 = G13 * G129;
    magmaFloatComplex G131 = G128 + G130;
    magmaFloatComplex G132 = s_track[9] * s_track[9];
    magmaFloatComplex G133 = G13 * G132;
    magmaFloatComplex G134 = G131 + G133;
    magmaFloatComplex G135 = s_track[10] * s_track[10];
    magmaFloatComplex G136 = G13 * G135;
    magmaFloatComplex G137 = G134 + G136;
    magmaFloatComplex G138 = G25 * s_track[1];
    magmaFloatComplex G139 = G40 * s_track[2];
    magmaFloatComplex G140 = G138 + G139;
    magmaFloatComplex G141 = G9 * s_track[1];
    magmaFloatComplex G142 = G140 + G141;
    magmaFloatComplex G143 = G53 * s_track[3];
    magmaFloatComplex G144 = G142 + G143;
    magmaFloatComplex G145 = G64 * s_track[4];
    magmaFloatComplex G146 = G144 + G145;
    magmaFloatComplex G147 = G74 * s_track[5];
    magmaFloatComplex G148 = G146 + G147;
    magmaFloatComplex G149 = G81 * s_track[6];
    magmaFloatComplex G150 = G148 + G149;
    magmaFloatComplex G151 = G88 * s_track[7];
    magmaFloatComplex G152 = G150 + G151;
    magmaFloatComplex G153 = G94 * s_track[8];
    magmaFloatComplex G154 = G152 + G153;
    magmaFloatComplex G155 = G99 * s_track[9];
    magmaFloatComplex G156 = G154 + G155;
    magmaFloatComplex G157 = G103 * s_track[10];
    magmaFloatComplex G158 = G156 + G157;
    magmaFloatComplex G159 = G25 * s_track[2];
    magmaFloatComplex G160 = G4 * G108;
    magmaFloatComplex G161 = G159 + G160;
    magmaFloatComplex G162 = G40 * s_track[3];
    magmaFloatComplex G163 = G161 + G162;
    magmaFloatComplex G164 = G53 * s_track[4];
    magmaFloatComplex G165 = G163 + G164;
    magmaFloatComplex G166 = G9 * s_track[2];
    magmaFloatComplex G167 = G165 + G166;
    magmaFloatComplex G168 = G64 * s_track[5];
    magmaFloatComplex G169 = G167 + G168;
    magmaFloatComplex G170 = G74 * s_track[6];
    magmaFloatComplex G171 = G169 + G170;
    magmaFloatComplex G172 = G81 * s_track[7];
    magmaFloatComplex G173 = G171 + G172;
    magmaFloatComplex G174 = G88 * s_track[8];
    magmaFloatComplex G175 = G173 + G174;
    magmaFloatComplex G176 = G94 * s_track[9];
    magmaFloatComplex G177 = G175 + G176;
    magmaFloatComplex G178 = G99 * s_track[10];
    magmaFloatComplex G179 = G177 + G178;
    magmaFloatComplex G180 = G25 * s_track[3];
    magmaFloatComplex G181 = G180 + G139;
    magmaFloatComplex G182 = G40 * s_track[4];
    magmaFloatComplex G183 = G181 + G182;
    magmaFloatComplex G184 = G53 * s_track[5];
    magmaFloatComplex G185 = G183 + G184;
    magmaFloatComplex G186 = G64 * s_track[6];
    magmaFloatComplex G187 = G185 + G186;
    magmaFloatComplex G188 = G9 * s_track[3];
    magmaFloatComplex G189 = G187 + G188;
    magmaFloatComplex G190 = G74 * s_track[7];
    magmaFloatComplex G191 = G189 + G190;
    magmaFloatComplex G192 = G81 * s_track[8];
    magmaFloatComplex G193 = G191 + G192;
    magmaFloatComplex G194 = G88 * s_track[9];
    magmaFloatComplex G195 = G193 + G194;
    magmaFloatComplex G196 = G94 * s_track[10];
    magmaFloatComplex G197 = G195 + G196;
    magmaFloatComplex G198 = G25 * s_track[4];
    magmaFloatComplex G199 = G198 + G162;
    magmaFloatComplex G200 = G40 * s_track[5];
    magmaFloatComplex G201 = G199 + G200;
    magmaFloatComplex G202 = G4 * G111;
    magmaFloatComplex G203 = G201 + G202;
    magmaFloatComplex G204 = G53 * s_track[6];
    magmaFloatComplex G205 = G203 + G204;
    magmaFloatComplex G206 = G64 * s_track[7];
    magmaFloatComplex G207 = G205 + G206;
    magmaFloatComplex G208 = G74 * s_track[8];
    magmaFloatComplex G209 = G207 + G208;
    magmaFloatComplex G210 = G9 * s_track[4];
    magmaFloatComplex G211 = G209 + G210;
    magmaFloatComplex G212 = G81 * s_track[9];
    magmaFloatComplex G213 = G211 + G212;
    magmaFloatComplex G214 = G88 * s_track[10];
    magmaFloatComplex G215 = G213 + G214;
    magmaFloatComplex G216 = G25 * s_track[5];
    magmaFloatComplex G217 = G216 + G182;
    magmaFloatComplex G218 = G40 * s_track[6];
    magmaFloatComplex G219 = G217 + G218;
    magmaFloatComplex G220 = G219 + G143;
    magmaFloatComplex G221 = G53 * s_track[7];
    magmaFloatComplex G222 = G220 + G221;
    magmaFloatComplex G223 = G64 * s_track[8];
    magmaFloatComplex G224 = G222 + G223;
    magmaFloatComplex G225 = G74 * s_track[9];
    magmaFloatComplex G226 = G224 + G225;
    magmaFloatComplex G227 = G81 * s_track[10];
    magmaFloatComplex G228 = G226 + G227;
    magmaFloatComplex G229 = G9 * s_track[5];
    magmaFloatComplex G230 = G228 + G229;
    magmaFloatComplex G231 = G25 * s_track[6];
    magmaFloatComplex G232 = G231 + G200;
    magmaFloatComplex G233 = G40 * s_track[7];
    magmaFloatComplex G234 = G232 + G233;
    magmaFloatComplex G235 = G234 + G164;
    magmaFloatComplex G236 = G53 * s_track[8];
    magmaFloatComplex G237 = G235 + G236;
    magmaFloatComplex G238 = G4 * G114;
    magmaFloatComplex G239 = G237 + G238;
    magmaFloatComplex G240 = G64 * s_track[9];
    magmaFloatComplex G241 = G239 + G240;
    magmaFloatComplex G242 = G74 * s_track[10];
    magmaFloatComplex G243 = G241 + G242;
    magmaFloatComplex G244 = G9 * s_track[6];
    magmaFloatComplex G245 = G243 + G244;
    magmaFloatComplex G246 = G25 * s_track[7];
    magmaFloatComplex G247 = G246 + G218;
    magmaFloatComplex G248 = G40 * s_track[8];
    magmaFloatComplex G249 = G247 + G248;
    magmaFloatComplex G250 = G249 + G184;
    magmaFloatComplex G251 = G53 * s_track[9];
    magmaFloatComplex G252 = G250 + G251;
    magmaFloatComplex G253 = G252 + G145;
    magmaFloatComplex G254 = G64 * s_track[10];
    magmaFloatComplex G255 = G253 + G254;
    magmaFloatComplex G256 = G9 * s_track[7];
    magmaFloatComplex G257 = G255 + G256;
    magmaFloatComplex G258 = G25 * s_track[8];
    magmaFloatComplex G259 = G258 + G233;
    magmaFloatComplex G260 = G40 * s_track[9];
    magmaFloatComplex G261 = G259 + G260;
    magmaFloatComplex G262 = G261 + G204;
    magmaFloatComplex G263 = G53 * s_track[10];
    magmaFloatComplex G264 = G262 + G263;
    magmaFloatComplex G265 = G264 + G168;
    magmaFloatComplex G266 = G4 * G117;
    magmaFloatComplex G267 = G265 + G266;
    magmaFloatComplex G268 = G9 * s_track[8];
    magmaFloatComplex G269 = G267 + G268;
    magmaFloatComplex G270 = G25 * s_track[9];
    magmaFloatComplex G271 = G270 + G248;
    magmaFloatComplex G272 = G40 * s_track[10];
    magmaFloatComplex G273 = G271 + G272;
    magmaFloatComplex G274 = G273 + G221;
    magmaFloatComplex G275 = G274 + G186;
    magmaFloatComplex G276 = G275 + G147;
    magmaFloatComplex G277 = G9 * s_track[9];
    magmaFloatComplex G278 = G276 + G277;
    magmaFloatComplex G279 = G4 * s_track[0];
    magmaFloatComplex G280 = G279 + G40;
    magmaFloatComplex G281 = G280 + G53;
    magmaFloatComplex G282 = G281 + G64;
    magmaFloatComplex G283 = G282 + G74;
    magmaFloatComplex G284 = G283 + G81;
    magmaFloatComplex G285 = G284 + G88;
    magmaFloatComplex G286 = G285 + G94;
    magmaFloatComplex G287 = G286 + G99;
    magmaFloatComplex G288 = G287 + G103;
    magmaFloatComplex G289 = G13 * s_track[10];
    magmaFloatComplex G290 = G288 + G289;
    magmaFloatComplex G291 = G290 + G9;

    r_cgesvA[0] = G10;
    r_cgesvA[1] = G14;
    r_cgesvA[2] = G15;
    r_cgesvA[3] = G16;
    r_cgesvA[4] = G17;
    r_cgesvA[5] = G18;
    r_cgesvA[6] = G19;
    r_cgesvA[7] = G20;
    r_cgesvA[8] = G21;
    r_cgesvA[9] = G22;
    r_cgesvA[10] = G4;
    r_cgesvB[0] =G137;

    r_cgesvA[11] = G24;
    r_cgesvA[12] = G27;
    r_cgesvA[13] = G29;
    r_cgesvA[14] = G30;
    r_cgesvA[15] = G31;
    r_cgesvA[16] = G32;
    r_cgesvA[17] = G33;
    r_cgesvA[18] = G34;
    r_cgesvA[19] = G35;
    r_cgesvA[20] = G37;
    r_cgesvA[21] = G13;
    r_cgesvB[1] =G158;

    r_cgesvA[22] = G39;
    r_cgesvA[23] = G41;
    r_cgesvA[24] = G43;
    r_cgesvA[25] = G44;
    r_cgesvA[26] = G46;
    r_cgesvA[27] = G47;
    r_cgesvA[28] = G48;
    r_cgesvA[29] = G49;
    r_cgesvA[30] = G50;
    r_cgesvA[31] = G20;
    r_cgesvA[32] = G13;
    r_cgesvB[2] =G179;

    r_cgesvA[33] = G52;
    r_cgesvA[34] = G54;
    r_cgesvA[35] = G44;
    r_cgesvA[36] = G56;
    r_cgesvA[37] = G57;
    r_cgesvA[38] = G58;
    r_cgesvA[39] = G60;
    r_cgesvA[40] = G61;
    r_cgesvA[41] = G18;
    r_cgesvA[42] = G19;
    r_cgesvA[43] = G13;
    r_cgesvB[3] =G197;

    r_cgesvA[44] = G63;
    r_cgesvA[45] = G65;
    r_cgesvA[46] = G66;
    r_cgesvA[47] = G57;
    r_cgesvA[48] = G68;
    r_cgesvA[49] = G69;
    r_cgesvA[50] = G70;
    r_cgesvA[51] = G64;
    r_cgesvA[52] = G71;
    r_cgesvA[53] = G18;
    r_cgesvA[54] = G13;
    r_cgesvB[4] =G215;

    r_cgesvA[55] = G73;
    r_cgesvA[56] = G75;
    r_cgesvA[57] = G76;
    r_cgesvA[58] = G58;
    r_cgesvA[59] = G69;
    r_cgesvA[60] = G78;
    r_cgesvA[61] = G40;
    r_cgesvA[62] = G53;
    r_cgesvA[63] = G64;
    r_cgesvA[64] = G74;
    r_cgesvA[65] = G13;
    r_cgesvB[5] =G230;

    r_cgesvA[66] = G80;
    r_cgesvA[67] = G82;
    r_cgesvA[68] = G83;
    r_cgesvA[69] = G84;
    r_cgesvA[70] = G70;
    r_cgesvA[71] = G40;
    r_cgesvA[72] = G85;
    r_cgesvA[73] = G40;
    r_cgesvA[74] = G53;
    r_cgesvA[75] = G64;
    r_cgesvA[76] = G13;
    r_cgesvB[6] =G245;

    r_cgesvA[77] = G87;
    r_cgesvA[78] = G89;
    r_cgesvA[79] = G90;
    r_cgesvA[80] = G91;
    r_cgesvA[81] = G64;
    r_cgesvA[82] = G53;
    r_cgesvA[83] = G40;
    r_cgesvA[84] = G85;
    r_cgesvA[85] = G40;
    r_cgesvA[86] = G53;
    r_cgesvA[87] = G13;
    r_cgesvB[7] =G257;

    r_cgesvA[88] = G93;
    r_cgesvA[89] = G95;
    r_cgesvA[90] = G96;
    r_cgesvA[91] = G81;
    r_cgesvA[92] = G74;
    r_cgesvA[93] = G64;
    r_cgesvA[94] = G53;
    r_cgesvA[95] = G40;
    r_cgesvA[96] = G85;
    r_cgesvA[97] = G40;
    r_cgesvA[98] = G13;
    r_cgesvB[8] =G269;

    r_cgesvA[99] = G98;
    r_cgesvA[100] = G100;
    r_cgesvA[101] = G94;
    r_cgesvA[102] = G88;
    r_cgesvA[103] = G81;
    r_cgesvA[104] = G74;
    r_cgesvA[105] = G64;
    r_cgesvA[106] = G53;
    r_cgesvA[107] = G40;
    r_cgesvA[108] = G85;
    r_cgesvA[109] = G13;
    r_cgesvB[9] =G278;

    r_cgesvA[110] = G102;
    r_cgesvA[111] = G103;
    r_cgesvA[112] = G99;
    r_cgesvA[113] = G94;
    r_cgesvA[114] = G88;
    r_cgesvA[115] = G81;
    r_cgesvA[116] = G74;
    r_cgesvA[117] = G64;
    r_cgesvA[118] = G53;
    r_cgesvA[119] = G40;
    r_cgesvA[120] = G13;
    r_cgesvB[10] =G291;
  }
}

#endif
