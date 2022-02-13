#ifndef cpu_eval_HxHt_katsura10_h
#define cpu_eval_HxHt_katsura10_h
// ============================================================================
// partial derivative evaluations of the katsura10 problem for cpu HC computation
//
// Modifications
//    Chien  21-12-13:   Originally created
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
  void cpu_eval_HxHt_katsura10(
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
    magmaFloatComplex G106 = s_targetCoefs[0] - s_startCoefs[0];
    magmaFloatComplex G107 = G104 * G106;
    magmaFloatComplex G109 = s_targetCoefs[1] - s_startCoefs[1];
    magmaFloatComplex G110 = s_track[0] * G109;
    magmaFloatComplex G111 = G107 + G110;
    magmaFloatComplex G112 = s_track[1] * s_track[1];
    magmaFloatComplex G114 = s_targetCoefs[2] - s_startCoefs[2];
    magmaFloatComplex G115 = G112 * G114;
    magmaFloatComplex G116 = G111 + G115;
    magmaFloatComplex G117 = s_track[2] * s_track[2];
    magmaFloatComplex G118 = G117 * G114;
    magmaFloatComplex G119 = G116 + G118;
    magmaFloatComplex G120 = s_track[3] * s_track[3];
    magmaFloatComplex G121 = G120 * G114;
    magmaFloatComplex G122 = G119 + G121;
    magmaFloatComplex G123 = s_track[4] * s_track[4];
    magmaFloatComplex G124 = G123 * G114;
    magmaFloatComplex G125 = G122 + G124;
    magmaFloatComplex G126 = s_track[5] * s_track[5];
    magmaFloatComplex G127 = G126 * G114;
    magmaFloatComplex G128 = G125 + G127;
    magmaFloatComplex G129 = s_track[6] * s_track[6];
    magmaFloatComplex G130 = G129 * G114;
    magmaFloatComplex G131 = G128 + G130;
    magmaFloatComplex G132 = s_track[7] * s_track[7];
    magmaFloatComplex G133 = G132 * G114;
    magmaFloatComplex G134 = G131 + G133;
    magmaFloatComplex G135 = s_track[8] * s_track[8];
    magmaFloatComplex G136 = G135 * G114;
    magmaFloatComplex G137 = G134 + G136;
    magmaFloatComplex G138 = s_track[9] * s_track[9];
    magmaFloatComplex G139 = G138 * G114;
    magmaFloatComplex G140 = G137 + G139;
    magmaFloatComplex G141 = s_track[10] * s_track[10];
    magmaFloatComplex G142 = G141 * G114;
    magmaFloatComplex G143 = G140 + G142;
    magmaFloatComplex G144 = s_track[0] * G114;
    magmaFloatComplex G145 = s_track[1] * G144;
    magmaFloatComplex G146 = s_track[1] * G114;
    magmaFloatComplex G147 = s_track[2] * G146;
    magmaFloatComplex G148 = G145 + G147;
    magmaFloatComplex G149 = s_track[1] * G109;
    magmaFloatComplex G150 = G148 + G149;
    magmaFloatComplex G151 = s_track[2] * G114;
    magmaFloatComplex G152 = s_track[3] * G151;
    magmaFloatComplex G153 = G150 + G152;
    magmaFloatComplex G154 = s_track[3] * G114;
    magmaFloatComplex G155 = s_track[4] * G154;
    magmaFloatComplex G156 = G153 + G155;
    magmaFloatComplex G157 = s_track[4] * G114;
    magmaFloatComplex G158 = s_track[5] * G157;
    magmaFloatComplex G159 = G156 + G158;
    magmaFloatComplex G160 = s_track[5] * G114;
    magmaFloatComplex G161 = s_track[6] * G160;
    magmaFloatComplex G162 = G159 + G161;
    magmaFloatComplex G163 = s_track[6] * G114;
    magmaFloatComplex G164 = s_track[7] * G163;
    magmaFloatComplex G165 = G162 + G164;
    magmaFloatComplex G166 = s_track[7] * G114;
    magmaFloatComplex G167 = s_track[8] * G166;
    magmaFloatComplex G168 = G165 + G167;
    magmaFloatComplex G169 = s_track[8] * G114;
    magmaFloatComplex G170 = s_track[9] * G169;
    magmaFloatComplex G171 = G168 + G170;
    magmaFloatComplex G172 = s_track[9] * G114;
    magmaFloatComplex G173 = s_track[10] * G172;
    magmaFloatComplex G174 = G171 + G173;
    magmaFloatComplex G175 = s_track[2] * G144;
    magmaFloatComplex G176 = G112 * G106;
    magmaFloatComplex G177 = G175 + G176;
    magmaFloatComplex G178 = s_track[3] * G146;
    magmaFloatComplex G179 = G177 + G178;
    magmaFloatComplex G180 = s_track[4] * G151;
    magmaFloatComplex G181 = G179 + G180;
    magmaFloatComplex G182 = s_track[2] * G109;
    magmaFloatComplex G183 = G181 + G182;
    magmaFloatComplex G184 = s_track[5] * G154;
    magmaFloatComplex G185 = G183 + G184;
    magmaFloatComplex G186 = s_track[6] * G157;
    magmaFloatComplex G187 = G185 + G186;
    magmaFloatComplex G188 = s_track[7] * G160;
    magmaFloatComplex G189 = G187 + G188;
    magmaFloatComplex G190 = s_track[8] * G163;
    magmaFloatComplex G191 = G189 + G190;
    magmaFloatComplex G192 = s_track[9] * G166;
    magmaFloatComplex G193 = G191 + G192;
    magmaFloatComplex G194 = s_track[10] * G169;
    magmaFloatComplex G195 = G193 + G194;
    magmaFloatComplex G196 = s_track[3] * G144;
    magmaFloatComplex G197 = G196 + G147;
    magmaFloatComplex G198 = s_track[4] * G146;
    magmaFloatComplex G199 = G197 + G198;
    magmaFloatComplex G200 = s_track[5] * G151;
    magmaFloatComplex G201 = G199 + G200;
    magmaFloatComplex G202 = s_track[6] * G154;
    magmaFloatComplex G203 = G201 + G202;
    magmaFloatComplex G204 = s_track[3] * G109;
    magmaFloatComplex G205 = G203 + G204;
    magmaFloatComplex G206 = s_track[7] * G157;
    magmaFloatComplex G207 = G205 + G206;
    magmaFloatComplex G208 = s_track[8] * G160;
    magmaFloatComplex G209 = G207 + G208;
    magmaFloatComplex G210 = s_track[9] * G163;
    magmaFloatComplex G211 = G209 + G210;
    magmaFloatComplex G212 = s_track[10] * G166;
    magmaFloatComplex G213 = G211 + G212;
    magmaFloatComplex G214 = s_track[4] * G144;
    magmaFloatComplex G215 = G214 + G178;
    magmaFloatComplex G216 = s_track[5] * G146;
    magmaFloatComplex G217 = G215 + G216;
    magmaFloatComplex G218 = G117 * G106;
    magmaFloatComplex G219 = G217 + G218;
    magmaFloatComplex G220 = s_track[6] * G151;
    magmaFloatComplex G221 = G219 + G220;
    magmaFloatComplex G222 = s_track[7] * G154;
    magmaFloatComplex G223 = G221 + G222;
    magmaFloatComplex G224 = s_track[8] * G157;
    magmaFloatComplex G225 = G223 + G224;
    magmaFloatComplex G226 = s_track[4] * G109;
    magmaFloatComplex G227 = G225 + G226;
    magmaFloatComplex G228 = s_track[9] * G160;
    magmaFloatComplex G229 = G227 + G228;
    magmaFloatComplex G230 = s_track[10] * G163;
    magmaFloatComplex G231 = G229 + G230;
    magmaFloatComplex G232 = s_track[5] * G144;
    magmaFloatComplex G233 = G232 + G198;
    magmaFloatComplex G234 = s_track[6] * G146;
    magmaFloatComplex G235 = G233 + G234;
    magmaFloatComplex G236 = G235 + G152;
    magmaFloatComplex G237 = s_track[7] * G151;
    magmaFloatComplex G238 = G236 + G237;
    magmaFloatComplex G239 = s_track[8] * G154;
    magmaFloatComplex G240 = G238 + G239;
    magmaFloatComplex G241 = s_track[9] * G157;
    magmaFloatComplex G242 = G240 + G241;
    magmaFloatComplex G243 = s_track[10] * G160;
    magmaFloatComplex G244 = G242 + G243;
    magmaFloatComplex G245 = s_track[5] * G109;
    magmaFloatComplex G246 = G244 + G245;
    magmaFloatComplex G247 = s_track[6] * G144;
    magmaFloatComplex G248 = G247 + G216;
    magmaFloatComplex G249 = s_track[7] * G146;
    magmaFloatComplex G250 = G248 + G249;
    magmaFloatComplex G251 = G250 + G180;
    magmaFloatComplex G252 = s_track[8] * G151;
    magmaFloatComplex G253 = G251 + G252;
    magmaFloatComplex G254 = G120 * G106;
    magmaFloatComplex G255 = G253 + G254;
    magmaFloatComplex G256 = s_track[9] * G154;
    magmaFloatComplex G257 = G255 + G256;
    magmaFloatComplex G258 = s_track[10] * G157;
    magmaFloatComplex G259 = G257 + G258;
    magmaFloatComplex G260 = s_track[6] * G109;
    magmaFloatComplex G261 = G259 + G260;
    magmaFloatComplex G262 = s_track[7] * G144;
    magmaFloatComplex G263 = G262 + G234;
    magmaFloatComplex G264 = s_track[8] * G146;
    magmaFloatComplex G265 = G263 + G264;
    magmaFloatComplex G266 = G265 + G200;
    magmaFloatComplex G267 = s_track[9] * G151;
    magmaFloatComplex G268 = G266 + G267;
    magmaFloatComplex G269 = G268 + G155;
    magmaFloatComplex G270 = s_track[10] * G154;
    magmaFloatComplex G271 = G269 + G270;
    magmaFloatComplex G272 = s_track[7] * G109;
    magmaFloatComplex G273 = G271 + G272;
    magmaFloatComplex G274 = s_track[8] * G144;
    magmaFloatComplex G275 = G274 + G249;
    magmaFloatComplex G276 = s_track[9] * G146;
    magmaFloatComplex G277 = G275 + G276;
    magmaFloatComplex G278 = G277 + G220;
    magmaFloatComplex G279 = s_track[10] * G151;
    magmaFloatComplex G280 = G278 + G279;
    magmaFloatComplex G281 = G280 + G184;
    magmaFloatComplex G282 = G123 * G106;
    magmaFloatComplex G283 = G281 + G282;
    magmaFloatComplex G284 = s_track[8] * G109;
    magmaFloatComplex G285 = G283 + G284;
    magmaFloatComplex G286 = s_track[9] * G144;
    magmaFloatComplex G287 = G286 + G264;
    magmaFloatComplex G288 = s_track[10] * G146;
    magmaFloatComplex G289 = G287 + G288;
    magmaFloatComplex G290 = G289 + G237;
    magmaFloatComplex G291 = G290 + G202;
    magmaFloatComplex G292 = G291 + G158;
    magmaFloatComplex G293 = s_track[9] * G109;
    magmaFloatComplex G294 = G292 + G293;
    magmaFloatComplex G295 = s_track[0] * G106;
    magmaFloatComplex G296 = G295 + G146;
    magmaFloatComplex G297 = G296 + G151;
    magmaFloatComplex G298 = G297 + G154;
    magmaFloatComplex G299 = G298 + G157;
    magmaFloatComplex G300 = G299 + G160;
    magmaFloatComplex G301 = G300 + G163;
    magmaFloatComplex G302 = G301 + G166;
    magmaFloatComplex G303 = G302 + G169;
    magmaFloatComplex G304 = G303 + G172;
    magmaFloatComplex G305 = s_track[10] * G114;
    magmaFloatComplex G306 = G304 + G305;
    magmaFloatComplex G307 = G306 + G109;

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
    r_cgesvB[0] = -G143;

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
    r_cgesvB[1] = -G174;

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
    r_cgesvB[2] = -G195;

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
    r_cgesvB[3] = -G213;

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
    r_cgesvB[4] = -G231;

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
    r_cgesvB[5] = -G246;

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
    r_cgesvB[6] = -G261;

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
    r_cgesvB[7] = -G273;

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
    r_cgesvB[8] = -G285;

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
    r_cgesvB[9] = -G294;

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
    r_cgesvB[10] = -G307;
  }
}

#endif
