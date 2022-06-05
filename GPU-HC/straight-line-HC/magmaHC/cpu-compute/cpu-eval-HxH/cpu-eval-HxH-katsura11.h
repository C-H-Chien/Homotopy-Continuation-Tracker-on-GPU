#ifndef cpu_eval_HxH_katsura11_h
#define cpu_eval_HxH_katsura11_h
// ============================================================================
// partial derivative evaluations of the katsura11 problem for cpu HC computation
//
// Modifications
//    Chien  22-01-21:   Originally created
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
  void cpu_eval_HxH_katsura11(
      magma_int_t s, float t, int N, magmaFloatComplex* s_track,
      const magmaFloatComplex &C0, const magmaFloatComplex &C1,
      magmaFloatComplex* s_startCoefs, magmaFloatComplex* s_targetCoefs,
      magmaFloatComplex* r_cgesvA, magmaFloatComplex* r_cgesvB )
  {
magmaFloatComplex G0 = C1 * t;
magmaFloatComplex G1 = C0 + G0;
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
magmaFloatComplex G23 = s_track[10] * G13;
magmaFloatComplex G24 = s_track[1] + s_track[1];
magmaFloatComplex G25 = G13 * G24;
magmaFloatComplex G26 = G13 * s_track[0];
magmaFloatComplex G27 = G26 + G15;
magmaFloatComplex G28 = G27 + G9;
magmaFloatComplex G29 = G4 * G24;
magmaFloatComplex G30 = G29 + G16;
magmaFloatComplex G31 = G15 + G17;
magmaFloatComplex G32 = G16 + G18;
magmaFloatComplex G33 = G17 + G19;
magmaFloatComplex G34 = G18 + G20;
magmaFloatComplex G35 = G19 + G21;
magmaFloatComplex G36 = G20 + G22;
magmaFloatComplex G37 = G21 + G23;
magmaFloatComplex G38 = s_track[11] * G13;
magmaFloatComplex G39 = G22 + G38;
magmaFloatComplex G40 = s_track[2] + s_track[2];
magmaFloatComplex G41 = G13 * G40;
magmaFloatComplex G42 = G13 * s_track[1];
magmaFloatComplex G43 = G42 + G16;
magmaFloatComplex G44 = G26 + G17;
magmaFloatComplex G45 = G44 + G9;
magmaFloatComplex G46 = G42 + G18;
magmaFloatComplex G47 = G4 * G40;
magmaFloatComplex G48 = G47 + G19;
magmaFloatComplex G49 = G16 + G20;
magmaFloatComplex G50 = G17 + G21;
magmaFloatComplex G51 = G18 + G22;
magmaFloatComplex G52 = G19 + G23;
magmaFloatComplex G53 = G20 + G38;
magmaFloatComplex G54 = s_track[3] + s_track[3];
magmaFloatComplex G55 = G13 * G54;
magmaFloatComplex G56 = G13 * s_track[2];
magmaFloatComplex G57 = G56 + G17;
magmaFloatComplex G58 = G26 + G19;
magmaFloatComplex G59 = G58 + G9;
magmaFloatComplex G60 = G42 + G20;
magmaFloatComplex G61 = G56 + G21;
magmaFloatComplex G62 = G4 * G54;
magmaFloatComplex G63 = G62 + G22;
magmaFloatComplex G64 = G17 + G23;
magmaFloatComplex G65 = G18 + G38;
magmaFloatComplex G66 = s_track[4] + s_track[4];
magmaFloatComplex G67 = G13 * G66;
magmaFloatComplex G68 = G13 * s_track[3];
magmaFloatComplex G69 = G68 + G18;
magmaFloatComplex G70 = G56 + G19;
magmaFloatComplex G71 = G26 + G21;
magmaFloatComplex G72 = G71 + G9;
magmaFloatComplex G73 = G42 + G22;
magmaFloatComplex G74 = G56 + G23;
magmaFloatComplex G75 = G68 + G38;
magmaFloatComplex G76 = G4 * G66;
magmaFloatComplex G77 = s_track[5] + s_track[5];
magmaFloatComplex G78 = G13 * G77;
magmaFloatComplex G79 = G13 * s_track[4];
magmaFloatComplex G80 = G79 + G19;
magmaFloatComplex G81 = G68 + G20;
magmaFloatComplex G82 = G26 + G23;
magmaFloatComplex G83 = G82 + G9;
magmaFloatComplex G84 = G42 + G38;
magmaFloatComplex G85 = G4 * G77;
magmaFloatComplex G86 = s_track[6] + s_track[6];
magmaFloatComplex G87 = G13 * G86;
magmaFloatComplex G88 = G13 * s_track[5];
magmaFloatComplex G89 = G88 + G20;
magmaFloatComplex G90 = G79 + G21;
magmaFloatComplex G91 = G68 + G22;
magmaFloatComplex G92 = G26 + G9;
magmaFloatComplex G93 = s_track[7] + s_track[7];
magmaFloatComplex G94 = G13 * G93;
magmaFloatComplex G95 = G13 * s_track[6];
magmaFloatComplex G96 = G95 + G21;
magmaFloatComplex G97 = G88 + G22;
magmaFloatComplex G98 = G79 + G23;
magmaFloatComplex G99 = s_track[8] + s_track[8];
magmaFloatComplex G100 = G13 * G99;
magmaFloatComplex G101 = G13 * s_track[7];
magmaFloatComplex G102 = G101 + G22;
magmaFloatComplex G103 = G95 + G23;
magmaFloatComplex G104 = G88 + G38;
magmaFloatComplex G105 = s_track[9] + s_track[9];
magmaFloatComplex G106 = G13 * G105;
magmaFloatComplex G107 = G13 * s_track[8];
magmaFloatComplex G108 = G107 + G23;
magmaFloatComplex G109 = G101 + G38;
magmaFloatComplex G110 = s_track[10] + s_track[10];
magmaFloatComplex G111 = G13 * G110;
magmaFloatComplex G112 = G13 * s_track[9];
magmaFloatComplex G113 = G112 + G38;
magmaFloatComplex G114 = s_track[11] + s_track[11];
magmaFloatComplex G115 = G13 * G114;
magmaFloatComplex G116 = G13 * s_track[10];
magmaFloatComplex G117 = s_track[0] * s_track[0];
magmaFloatComplex G118 = G4 * G117;
magmaFloatComplex G119 = G9 * s_track[0];
magmaFloatComplex G120 = G118 + G119;
magmaFloatComplex G121 = s_track[1] * s_track[1];
magmaFloatComplex G122 = G13 * G121;
magmaFloatComplex G123 = G120 + G122;
magmaFloatComplex G124 = s_track[2] * s_track[2];
magmaFloatComplex G125 = G13 * G124;
magmaFloatComplex G126 = G123 + G125;
magmaFloatComplex G127 = s_track[3] * s_track[3];
magmaFloatComplex G128 = G13 * G127;
magmaFloatComplex G129 = G126 + G128;
magmaFloatComplex G130 = s_track[4] * s_track[4];
magmaFloatComplex G131 = G13 * G130;
magmaFloatComplex G132 = G129 + G131;
magmaFloatComplex G133 = s_track[5] * s_track[5];
magmaFloatComplex G134 = G13 * G133;
magmaFloatComplex G135 = G132 + G134;
magmaFloatComplex G136 = s_track[6] * s_track[6];
magmaFloatComplex G137 = G13 * G136;
magmaFloatComplex G138 = G135 + G137;
magmaFloatComplex G139 = s_track[7] * s_track[7];
magmaFloatComplex G140 = G13 * G139;
magmaFloatComplex G141 = G138 + G140;
magmaFloatComplex G142 = s_track[8] * s_track[8];
magmaFloatComplex G143 = G13 * G142;
magmaFloatComplex G144 = G141 + G143;
magmaFloatComplex G145 = s_track[9] * s_track[9];
magmaFloatComplex G146 = G13 * G145;
magmaFloatComplex G147 = G144 + G146;
magmaFloatComplex G148 = s_track[10] * s_track[10];
magmaFloatComplex G149 = G13 * G148;
magmaFloatComplex G150 = G147 + G149;
magmaFloatComplex G151 = s_track[11] * s_track[11];
magmaFloatComplex G152 = G13 * G151;
magmaFloatComplex G153 = G150 + G152;
magmaFloatComplex G154 = G26 * s_track[1];
magmaFloatComplex G155 = G42 * s_track[2];
magmaFloatComplex G156 = G154 + G155;
magmaFloatComplex G157 = G9 * s_track[1];
magmaFloatComplex G158 = G156 + G157;
magmaFloatComplex G159 = G56 * s_track[3];
magmaFloatComplex G160 = G158 + G159;
magmaFloatComplex G161 = G68 * s_track[4];
magmaFloatComplex G162 = G160 + G161;
magmaFloatComplex G163 = G79 * s_track[5];
magmaFloatComplex G164 = G162 + G163;
magmaFloatComplex G165 = G88 * s_track[6];
magmaFloatComplex G166 = G164 + G165;
magmaFloatComplex G167 = G95 * s_track[7];
magmaFloatComplex G168 = G166 + G167;
magmaFloatComplex G169 = G101 * s_track[8];
magmaFloatComplex G170 = G168 + G169;
magmaFloatComplex G171 = G107 * s_track[9];
magmaFloatComplex G172 = G170 + G171;
magmaFloatComplex G173 = G112 * s_track[10];
magmaFloatComplex G174 = G172 + G173;
magmaFloatComplex G175 = G116 * s_track[11];
magmaFloatComplex G176 = G174 + G175;
magmaFloatComplex G177 = G26 * s_track[2];
magmaFloatComplex G178 = G4 * G121;
magmaFloatComplex G179 = G177 + G178;
magmaFloatComplex G180 = G42 * s_track[3];
magmaFloatComplex G181 = G179 + G180;
magmaFloatComplex G182 = G56 * s_track[4];
magmaFloatComplex G183 = G181 + G182;
magmaFloatComplex G184 = G9 * s_track[2];
magmaFloatComplex G185 = G183 + G184;
magmaFloatComplex G186 = G68 * s_track[5];
magmaFloatComplex G187 = G185 + G186;
magmaFloatComplex G188 = G79 * s_track[6];
magmaFloatComplex G189 = G187 + G188;
magmaFloatComplex G190 = G88 * s_track[7];
magmaFloatComplex G191 = G189 + G190;
magmaFloatComplex G192 = G95 * s_track[8];
magmaFloatComplex G193 = G191 + G192;
magmaFloatComplex G194 = G101 * s_track[9];
magmaFloatComplex G195 = G193 + G194;
magmaFloatComplex G196 = G107 * s_track[10];
magmaFloatComplex G197 = G195 + G196;
magmaFloatComplex G198 = G112 * s_track[11];
magmaFloatComplex G199 = G197 + G198;
magmaFloatComplex G200 = G26 * s_track[3];
magmaFloatComplex G201 = G200 + G155;
magmaFloatComplex G202 = G42 * s_track[4];
magmaFloatComplex G203 = G201 + G202;
magmaFloatComplex G204 = G56 * s_track[5];
magmaFloatComplex G205 = G203 + G204;
magmaFloatComplex G206 = G68 * s_track[6];
magmaFloatComplex G207 = G205 + G206;
magmaFloatComplex G208 = G9 * s_track[3];
magmaFloatComplex G209 = G207 + G208;
magmaFloatComplex G210 = G79 * s_track[7];
magmaFloatComplex G211 = G209 + G210;
magmaFloatComplex G212 = G88 * s_track[8];
magmaFloatComplex G213 = G211 + G212;
magmaFloatComplex G214 = G95 * s_track[9];
magmaFloatComplex G215 = G213 + G214;
magmaFloatComplex G216 = G101 * s_track[10];
magmaFloatComplex G217 = G215 + G216;
magmaFloatComplex G218 = G107 * s_track[11];
magmaFloatComplex G219 = G217 + G218;
magmaFloatComplex G220 = G26 * s_track[4];
magmaFloatComplex G221 = G220 + G180;
magmaFloatComplex G222 = G42 * s_track[5];
magmaFloatComplex G223 = G221 + G222;
magmaFloatComplex G224 = G4 * G124;
magmaFloatComplex G225 = G223 + G224;
magmaFloatComplex G226 = G56 * s_track[6];
magmaFloatComplex G227 = G225 + G226;
magmaFloatComplex G228 = G68 * s_track[7];
magmaFloatComplex G229 = G227 + G228;
magmaFloatComplex G230 = G79 * s_track[8];
magmaFloatComplex G231 = G229 + G230;
magmaFloatComplex G232 = G9 * s_track[4];
magmaFloatComplex G233 = G231 + G232;
magmaFloatComplex G234 = G88 * s_track[9];
magmaFloatComplex G235 = G233 + G234;
magmaFloatComplex G236 = G95 * s_track[10];
magmaFloatComplex G237 = G235 + G236;
magmaFloatComplex G238 = G101 * s_track[11];
magmaFloatComplex G239 = G237 + G238;
magmaFloatComplex G240 = G26 * s_track[5];
magmaFloatComplex G241 = G240 + G202;
magmaFloatComplex G242 = G42 * s_track[6];
magmaFloatComplex G243 = G241 + G242;
magmaFloatComplex G244 = G243 + G159;
magmaFloatComplex G245 = G56 * s_track[7];
magmaFloatComplex G246 = G244 + G245;
magmaFloatComplex G247 = G68 * s_track[8];
magmaFloatComplex G248 = G246 + G247;
magmaFloatComplex G249 = G79 * s_track[9];
magmaFloatComplex G250 = G248 + G249;
magmaFloatComplex G251 = G88 * s_track[10];
magmaFloatComplex G252 = G250 + G251;
magmaFloatComplex G253 = G9 * s_track[5];
magmaFloatComplex G254 = G252 + G253;
magmaFloatComplex G255 = G95 * s_track[11];
magmaFloatComplex G256 = G254 + G255;
magmaFloatComplex G257 = G26 * s_track[6];
magmaFloatComplex G258 = G257 + G222;
magmaFloatComplex G259 = G42 * s_track[7];
magmaFloatComplex G260 = G258 + G259;
magmaFloatComplex G261 = G260 + G182;
magmaFloatComplex G262 = G56 * s_track[8];
magmaFloatComplex G263 = G261 + G262;
magmaFloatComplex G264 = G4 * G127;
magmaFloatComplex G265 = G263 + G264;
magmaFloatComplex G266 = G68 * s_track[9];
magmaFloatComplex G267 = G265 + G266;
magmaFloatComplex G268 = G79 * s_track[10];
magmaFloatComplex G269 = G267 + G268;
magmaFloatComplex G270 = G88 * s_track[11];
magmaFloatComplex G271 = G269 + G270;
magmaFloatComplex G272 = G9 * s_track[6];
magmaFloatComplex G273 = G271 + G272;
magmaFloatComplex G274 = G26 * s_track[7];
magmaFloatComplex G275 = G274 + G242;
magmaFloatComplex G276 = G42 * s_track[8];
magmaFloatComplex G277 = G275 + G276;
magmaFloatComplex G278 = G277 + G204;
magmaFloatComplex G279 = G56 * s_track[9];
magmaFloatComplex G280 = G278 + G279;
magmaFloatComplex G281 = G280 + G161;
magmaFloatComplex G282 = G68 * s_track[10];
magmaFloatComplex G283 = G281 + G282;
magmaFloatComplex G284 = G79 * s_track[11];
magmaFloatComplex G285 = G283 + G284;
magmaFloatComplex G286 = G9 * s_track[7];
magmaFloatComplex G287 = G285 + G286;
magmaFloatComplex G288 = G26 * s_track[8];
magmaFloatComplex G289 = G288 + G259;
magmaFloatComplex G290 = G42 * s_track[9];
magmaFloatComplex G291 = G289 + G290;
magmaFloatComplex G292 = G291 + G226;
magmaFloatComplex G293 = G56 * s_track[10];
magmaFloatComplex G294 = G292 + G293;
magmaFloatComplex G295 = G294 + G186;
magmaFloatComplex G296 = G68 * s_track[11];
magmaFloatComplex G297 = G295 + G296;
magmaFloatComplex G298 = G4 * G130;
magmaFloatComplex G299 = G297 + G298;
magmaFloatComplex G300 = G9 * s_track[8];
magmaFloatComplex G301 = G299 + G300;
magmaFloatComplex G302 = G26 * s_track[9];
magmaFloatComplex G303 = G302 + G276;
magmaFloatComplex G304 = G42 * s_track[10];
magmaFloatComplex G305 = G303 + G304;
magmaFloatComplex G306 = G305 + G245;
magmaFloatComplex G307 = G56 * s_track[11];
magmaFloatComplex G308 = G306 + G307;
magmaFloatComplex G309 = G308 + G206;
magmaFloatComplex G310 = G309 + G163;
magmaFloatComplex G311 = G9 * s_track[9];
magmaFloatComplex G312 = G310 + G311;
magmaFloatComplex G313 = G26 * s_track[10];
magmaFloatComplex G314 = G313 + G290;
magmaFloatComplex G315 = G42 * s_track[11];
magmaFloatComplex G316 = G314 + G315;
magmaFloatComplex G317 = G316 + G262;
magmaFloatComplex G318 = G317 + G228;
magmaFloatComplex G319 = G318 + G188;
magmaFloatComplex G320 = G4 * G133;
magmaFloatComplex G321 = G319 + G320;
magmaFloatComplex G322 = G9 * s_track[10];
magmaFloatComplex G323 = G321 + G322;
magmaFloatComplex G324 = G4 * s_track[0];
magmaFloatComplex G325 = G324 + G42;
magmaFloatComplex G326 = G325 + G56;
magmaFloatComplex G327 = G326 + G68;
magmaFloatComplex G328 = G327 + G79;
magmaFloatComplex G329 = G328 + G88;
magmaFloatComplex G330 = G329 + G95;
magmaFloatComplex G331 = G330 + G101;
magmaFloatComplex G332 = G331 + G107;
magmaFloatComplex G333 = G332 + G112;
magmaFloatComplex G334 = G333 + G116;
magmaFloatComplex G335 = G13 * s_track[11];
magmaFloatComplex G336 = G334 + G335;
magmaFloatComplex G337 = G336 + G9;

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
r_cgesvA[10] = G23;
r_cgesvA[11] = G4;
r_cgesvB[0] =G153;

r_cgesvA[12] = G25;
r_cgesvA[13] = G28;
r_cgesvA[14] = G30;
r_cgesvA[15] = G31;
r_cgesvA[16] = G32;
r_cgesvA[17] = G33;
r_cgesvA[18] = G34;
r_cgesvA[19] = G35;
r_cgesvA[20] = G36;
r_cgesvA[21] = G37;
r_cgesvA[22] = G39;
r_cgesvA[23] = G13;
r_cgesvB[1] =G176;

r_cgesvA[24] = G41;
r_cgesvA[25] = G43;
r_cgesvA[26] = G45;
r_cgesvA[27] = G46;
r_cgesvA[28] = G48;
r_cgesvA[29] = G49;
r_cgesvA[30] = G50;
r_cgesvA[31] = G51;
r_cgesvA[32] = G52;
r_cgesvA[33] = G53;
r_cgesvA[34] = G21;
r_cgesvA[35] = G13;
r_cgesvB[2] =G199;

r_cgesvA[36] = G55;
r_cgesvA[37] = G57;
r_cgesvA[38] = G46;
r_cgesvA[39] = G59;
r_cgesvA[40] = G60;
r_cgesvA[41] = G61;
r_cgesvA[42] = G63;
r_cgesvA[43] = G64;
r_cgesvA[44] = G65;
r_cgesvA[45] = G19;
r_cgesvA[46] = G20;
r_cgesvA[47] = G13;
r_cgesvB[3] =G219;

r_cgesvA[48] = G67;
r_cgesvA[49] = G69;
r_cgesvA[50] = G70;
r_cgesvA[51] = G60;
r_cgesvA[52] = G72;
r_cgesvA[53] = G73;
r_cgesvA[54] = G74;
r_cgesvA[55] = G75;
r_cgesvA[56] = G76;
r_cgesvA[57] = G18;
r_cgesvA[58] = G19;
r_cgesvA[59] = G13;
r_cgesvB[4] =G239;

r_cgesvA[60] = G78;
r_cgesvA[61] = G80;
r_cgesvA[62] = G81;
r_cgesvA[63] = G61;
r_cgesvA[64] = G73;
r_cgesvA[65] = G83;
r_cgesvA[66] = G84;
r_cgesvA[67] = G56;
r_cgesvA[68] = G68;
r_cgesvA[69] = G79;
r_cgesvA[70] = G85;
r_cgesvA[71] = G13;
r_cgesvB[5] =G256;

r_cgesvA[72] = G87;
r_cgesvA[73] = G89;
r_cgesvA[74] = G90;
r_cgesvA[75] = G91;
r_cgesvA[76] = G74;
r_cgesvA[77] = G84;
r_cgesvA[78] = G92;
r_cgesvA[79] = G42;
r_cgesvA[80] = G56;
r_cgesvA[81] = G68;
r_cgesvA[82] = G79;
r_cgesvA[83] = G13;
r_cgesvB[6] =G273;

r_cgesvA[84] = G94;
r_cgesvA[85] = G96;
r_cgesvA[86] = G97;
r_cgesvA[87] = G98;
r_cgesvA[88] = G75;
r_cgesvA[89] = G56;
r_cgesvA[90] = G42;
r_cgesvA[91] = G92;
r_cgesvA[92] = G42;
r_cgesvA[93] = G56;
r_cgesvA[94] = G68;
r_cgesvA[95] = G13;
r_cgesvB[7] =G287;

r_cgesvA[96] = G100;
r_cgesvA[97] = G102;
r_cgesvA[98] = G103;
r_cgesvA[99] = G104;
r_cgesvA[100] = G79;
r_cgesvA[101] = G68;
r_cgesvA[102] = G56;
r_cgesvA[103] = G42;
r_cgesvA[104] = G92;
r_cgesvA[105] = G42;
r_cgesvA[106] = G56;
r_cgesvA[107] = G13;
r_cgesvB[8] =G301;

r_cgesvA[108] = G106;
r_cgesvA[109] = G108;
r_cgesvA[110] = G109;
r_cgesvA[111] = G95;
r_cgesvA[112] = G88;
r_cgesvA[113] = G79;
r_cgesvA[114] = G68;
r_cgesvA[115] = G56;
r_cgesvA[116] = G42;
r_cgesvA[117] = G92;
r_cgesvA[118] = G42;
r_cgesvA[119] = G13;
r_cgesvB[9] =G312;

r_cgesvA[120] = G111;
r_cgesvA[121] = G113;
r_cgesvA[122] = G107;
r_cgesvA[123] = G101;
r_cgesvA[124] = G95;
r_cgesvA[125] = G88;
r_cgesvA[126] = G79;
r_cgesvA[127] = G68;
r_cgesvA[128] = G56;
r_cgesvA[129] = G42;
r_cgesvA[130] = G92;
r_cgesvA[131] = G13;
r_cgesvB[10] =G323;

r_cgesvA[132] = G115;
r_cgesvA[133] = G116;
r_cgesvA[134] = G112;
r_cgesvA[135] = G107;
r_cgesvA[136] = G101;
r_cgesvA[137] = G95;
r_cgesvA[138] = G88;
r_cgesvA[139] = G79;
r_cgesvA[140] = G68;
r_cgesvA[141] = G56;
r_cgesvA[142] = G42;
r_cgesvA[143] = G13;
r_cgesvB[11] =G337;

  }
}

#endif
