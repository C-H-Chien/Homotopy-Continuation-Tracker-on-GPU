#ifndef cpu_eval_HxH_3vTrg_h
#define cpu_eval_HxH_3vTrg_h
// ============================================================================
// partial derivative evaluations of the 4-view triangulation problem
//
// Modifications
//    Chien  21-12-29:   Originally created
//
// ============================================================================
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

  extern "C"
  void cpu_eval_HxH_3vTrg(
      magma_int_t s, float t, int N, magmaFloatComplex* s_track,
      const magmaFloatComplex &C0, const magmaFloatComplex &C1, 
      const magmaFloatComplex &C2, const magmaFloatComplex &C3,
      magmaFloatComplex* s_startCoefs, magmaFloatComplex* s_targetCoefs,
      magmaFloatComplex* r_cgesvA, magmaFloatComplex* r_cgesvB )
  {
    magmaFloatComplex G0 = C3 * t;
magmaFloatComplex G1 = C2 + G0;
magmaFloatComplex G2 = G1 * s_startCoefs[0];
magmaFloatComplex G3 = t * s_targetCoefs[0];
magmaFloatComplex G4 = G2 + G3;
magmaFloatComplex G5 = G4 * s_track[6];
magmaFloatComplex G6 = G1 * s_startCoefs[1];
magmaFloatComplex G7 = t * s_targetCoefs[1];
magmaFloatComplex G8 = G6 + G7;
magmaFloatComplex G9 = G8 * s_track[6];
magmaFloatComplex G10 = G1 * s_startCoefs[18];
magmaFloatComplex G11 = t * s_targetCoefs[18];
magmaFloatComplex G12 = G10 + G11;
magmaFloatComplex G13 = G12 * s_track[8];
magmaFloatComplex G14 = G1 * s_startCoefs[19];
magmaFloatComplex G15 = t * s_targetCoefs[19];
magmaFloatComplex G16 = G14 + G15;
magmaFloatComplex G17 = G16 * s_track[8];
magmaFloatComplex G18 = G4 * s_track[2];
magmaFloatComplex G19 = G8 * s_track[3];
magmaFloatComplex G20 = G18 + G19;
magmaFloatComplex G21 = G1 * s_startCoefs[2];
magmaFloatComplex G22 = t * s_targetCoefs[2];
magmaFloatComplex G23 = G21 + G22;
magmaFloatComplex G24 = G20 + G23;
magmaFloatComplex G25 = G12 * s_track[4];
magmaFloatComplex G26 = G16 * s_track[5];
magmaFloatComplex G27 = G25 + G26;
magmaFloatComplex G28 = G1 * s_startCoefs[20];
magmaFloatComplex G29 = t * s_targetCoefs[20];
magmaFloatComplex G30 = G28 + G29;
magmaFloatComplex G31 = G27 + G30;
magmaFloatComplex G32 = G1 * s_startCoefs[3];
magmaFloatComplex G33 = t * s_targetCoefs[3];
magmaFloatComplex G34 = G32 + G33;
magmaFloatComplex G35 = G34 * s_track[6];
magmaFloatComplex G36 = G1 * s_startCoefs[4];
magmaFloatComplex G37 = t * s_targetCoefs[4];
magmaFloatComplex G38 = G36 + G37;
magmaFloatComplex G39 = G38 * s_track[6];
magmaFloatComplex G40 = G1 * s_startCoefs[21];
magmaFloatComplex G41 = t * s_targetCoefs[21];
magmaFloatComplex G42 = G40 + G41;
magmaFloatComplex G43 = G42 * s_track[8];
magmaFloatComplex G44 = G1 * s_startCoefs[22];
magmaFloatComplex G45 = t * s_targetCoefs[22];
magmaFloatComplex G46 = G44 + G45;
magmaFloatComplex G47 = G46 * s_track[8];
magmaFloatComplex G48 = G34 * s_track[2];
magmaFloatComplex G49 = G38 * s_track[3];
magmaFloatComplex G50 = G48 + G49;
magmaFloatComplex G51 = G1 * s_startCoefs[5];
magmaFloatComplex G52 = t * s_targetCoefs[5];
magmaFloatComplex G53 = G51 + G52;
magmaFloatComplex G54 = G50 + G53;
magmaFloatComplex G55 = G42 * s_track[4];
magmaFloatComplex G56 = G46 * s_track[5];
magmaFloatComplex G57 = G55 + G56;
magmaFloatComplex G58 = G1 * s_startCoefs[23];
magmaFloatComplex G59 = t * s_targetCoefs[23];
magmaFloatComplex G60 = G58 + G59;
magmaFloatComplex G61 = G57 + G60;
magmaFloatComplex G62 = G1 * s_startCoefs[9];
magmaFloatComplex G63 = t * s_targetCoefs[9];
magmaFloatComplex G64 = G62 + G63;
magmaFloatComplex G65 = G64 * s_track[7];
magmaFloatComplex G66 = G1 * s_startCoefs[10];
magmaFloatComplex G67 = t * s_targetCoefs[10];
magmaFloatComplex G68 = G66 + G67;
magmaFloatComplex G69 = G68 * s_track[7];
magmaFloatComplex G70 = G4 * s_track[0];
magmaFloatComplex G71 = G34 * s_track[1];
magmaFloatComplex G72 = G70 + G71;
magmaFloatComplex G73 = G1 * s_startCoefs[6];
magmaFloatComplex G74 = t * s_targetCoefs[6];
magmaFloatComplex G75 = G73 + G74;
magmaFloatComplex G76 = G72 + G75;
magmaFloatComplex G77 = G64 * s_track[4];
magmaFloatComplex G78 = G68 * s_track[5];
magmaFloatComplex G79 = G77 + G78;
magmaFloatComplex G80 = G1 * s_startCoefs[11];
magmaFloatComplex G81 = t * s_targetCoefs[11];
magmaFloatComplex G82 = G80 + G81;
magmaFloatComplex G83 = G79 + G82;
magmaFloatComplex G84 = G1 * s_startCoefs[12];
magmaFloatComplex G85 = t * s_targetCoefs[12];
magmaFloatComplex G86 = G84 + G85;
magmaFloatComplex G87 = G86 * s_track[7];
magmaFloatComplex G88 = G1 * s_startCoefs[13];
magmaFloatComplex G89 = t * s_targetCoefs[13];
magmaFloatComplex G90 = G88 + G89;
magmaFloatComplex G91 = G90 * s_track[7];
magmaFloatComplex G92 = G8 * s_track[0];
magmaFloatComplex G93 = G38 * s_track[1];
magmaFloatComplex G94 = G92 + G93;
magmaFloatComplex G95 = G1 * s_startCoefs[7];
magmaFloatComplex G96 = t * s_targetCoefs[7];
magmaFloatComplex G97 = G95 + G96;
magmaFloatComplex G98 = G94 + G97;
magmaFloatComplex G99 = G86 * s_track[4];
magmaFloatComplex G100 = G90 * s_track[5];
magmaFloatComplex G101 = G99 + G100;
magmaFloatComplex G102 = G1 * s_startCoefs[14];
magmaFloatComplex G103 = t * s_targetCoefs[14];
magmaFloatComplex G104 = G102 + G103;
magmaFloatComplex G105 = G101 + G104;
magmaFloatComplex G106 = G64 * s_track[2];
magmaFloatComplex G107 = G86 * s_track[3];
magmaFloatComplex G108 = G106 + G107;
magmaFloatComplex G109 = G1 * s_startCoefs[15];
magmaFloatComplex G110 = t * s_targetCoefs[15];
magmaFloatComplex G111 = G109 + G110;
magmaFloatComplex G112 = G108 + G111;
magmaFloatComplex G113 = G12 * s_track[0];
magmaFloatComplex G114 = G42 * s_track[1];
magmaFloatComplex G115 = G113 + G114;
magmaFloatComplex G116 = G1 * s_startCoefs[24];
magmaFloatComplex G117 = t * s_targetCoefs[24];
magmaFloatComplex G118 = G116 + G117;
magmaFloatComplex G119 = G115 + G118;
magmaFloatComplex G120 = G68 * s_track[2];
magmaFloatComplex G121 = G90 * s_track[3];
magmaFloatComplex G122 = G120 + G121;
magmaFloatComplex G123 = G1 * s_startCoefs[16];
magmaFloatComplex G124 = t * s_targetCoefs[16];
magmaFloatComplex G125 = G123 + G124;
magmaFloatComplex G126 = G122 + G125;
magmaFloatComplex G127 = G16 * s_track[0];
magmaFloatComplex G128 = G46 * s_track[1];
magmaFloatComplex G129 = G127 + G128;
magmaFloatComplex G130 = G1 * s_startCoefs[25];
magmaFloatComplex G131 = t * s_targetCoefs[25];
magmaFloatComplex G132 = G130 + G131;
magmaFloatComplex G133 = G129 + G132;
magmaFloatComplex G134 = C0 * s_track[0];
magmaFloatComplex G135 = s_track[2] * s_track[6];
magmaFloatComplex G136 = G135 * G4;
magmaFloatComplex G137 = G134 + G136;
magmaFloatComplex G138 = s_track[3] * s_track[6];
magmaFloatComplex G139 = G138 * G8;
magmaFloatComplex G140 = G137 + G139;
magmaFloatComplex G141 = s_track[4] * s_track[8];
magmaFloatComplex G142 = G141 * G12;
magmaFloatComplex G143 = G140 + G142;
magmaFloatComplex G144 = s_track[5] * s_track[8];
magmaFloatComplex G145 = G144 * G16;
magmaFloatComplex G146 = G143 + G145;
magmaFloatComplex G147 = s_track[6] * G23;
magmaFloatComplex G148 = G146 + G147;
magmaFloatComplex G149 = s_track[8] * G30;
magmaFloatComplex G150 = G148 + G149;
magmaFloatComplex G151 = G1 * s_startCoefs[27];
magmaFloatComplex G152 = t * s_targetCoefs[27];
magmaFloatComplex G153 = G151 + G152;
magmaFloatComplex G154 = C0 * G153;
magmaFloatComplex G155 = C3 * G154;
magmaFloatComplex G156 = G150 + G155;
magmaFloatComplex G157 = C0 * s_track[1];
magmaFloatComplex G158 = G135 * G34;
magmaFloatComplex G159 = G157 + G158;
magmaFloatComplex G160 = G138 * G38;
magmaFloatComplex G161 = G159 + G160;
magmaFloatComplex G162 = G141 * G42;
magmaFloatComplex G163 = G161 + G162;
magmaFloatComplex G164 = G144 * G46;
magmaFloatComplex G165 = G163 + G164;
magmaFloatComplex G166 = s_track[6] * G53;
magmaFloatComplex G167 = G165 + G166;
magmaFloatComplex G168 = s_track[8] * G60;
magmaFloatComplex G169 = G167 + G168;
magmaFloatComplex G170 = G1 * s_startCoefs[28];
magmaFloatComplex G171 = t * s_targetCoefs[28];
magmaFloatComplex G172 = G170 + G171;
magmaFloatComplex G173 = C0 * G172;
magmaFloatComplex G174 = C3 * G173;
magmaFloatComplex G175 = G169 + G174;
magmaFloatComplex G176 = s_track[0] * s_track[6];
magmaFloatComplex G177 = G176 * G4;
magmaFloatComplex G178 = s_track[1] * s_track[6];
magmaFloatComplex G179 = G178 * G34;
magmaFloatComplex G180 = G177 + G179;
magmaFloatComplex G181 = C0 * s_track[2];
magmaFloatComplex G182 = G180 + G181;
magmaFloatComplex G183 = s_track[4] * s_track[7];
magmaFloatComplex G184 = G183 * G64;
magmaFloatComplex G185 = G182 + G184;
magmaFloatComplex G186 = s_track[5] * s_track[7];
magmaFloatComplex G187 = G186 * G68;
magmaFloatComplex G188 = G185 + G187;
magmaFloatComplex G189 = s_track[6] * G75;
magmaFloatComplex G190 = G188 + G189;
magmaFloatComplex G191 = s_track[7] * G82;
magmaFloatComplex G192 = G190 + G191;
magmaFloatComplex G193 = G1 * s_startCoefs[29];
magmaFloatComplex G194 = t * s_targetCoefs[29];
magmaFloatComplex G195 = G193 + G194;
magmaFloatComplex G196 = C0 * G195;
magmaFloatComplex G197 = C3 * G196;
magmaFloatComplex G198 = G192 + G197;
magmaFloatComplex G199 = G176 * G8;
magmaFloatComplex G200 = G178 * G38;
magmaFloatComplex G201 = G199 + G200;
magmaFloatComplex G202 = C0 * s_track[3];
magmaFloatComplex G203 = G201 + G202;
magmaFloatComplex G204 = G183 * G86;
magmaFloatComplex G205 = G203 + G204;
magmaFloatComplex G206 = G186 * G90;
magmaFloatComplex G207 = G205 + G206;
magmaFloatComplex G208 = s_track[6] * G97;
magmaFloatComplex G209 = G207 + G208;
magmaFloatComplex G210 = s_track[7] * G104;
magmaFloatComplex G211 = G209 + G210;
magmaFloatComplex G212 = G1 * s_startCoefs[30];
magmaFloatComplex G213 = t * s_targetCoefs[30];
magmaFloatComplex G214 = G212 + G213;
magmaFloatComplex G215 = C0 * G214;
magmaFloatComplex G216 = C3 * G215;
magmaFloatComplex G217 = G211 + G216;
magmaFloatComplex G218 = s_track[0] * s_track[8];
magmaFloatComplex G219 = G218 * G12;
magmaFloatComplex G220 = s_track[1] * s_track[8];
magmaFloatComplex G221 = G220 * G42;
magmaFloatComplex G222 = G219 + G221;
magmaFloatComplex G223 = s_track[2] * s_track[7];
magmaFloatComplex G224 = G223 * G64;
magmaFloatComplex G225 = G222 + G224;
magmaFloatComplex G226 = s_track[3] * s_track[7];
magmaFloatComplex G227 = G226 * G86;
magmaFloatComplex G228 = G225 + G227;
magmaFloatComplex G229 = C0 * s_track[4];
magmaFloatComplex G230 = G228 + G229;
magmaFloatComplex G231 = s_track[7] * G111;
magmaFloatComplex G232 = G230 + G231;
magmaFloatComplex G233 = s_track[8] * G118;
magmaFloatComplex G234 = G232 + G233;
magmaFloatComplex G235 = G1 * s_startCoefs[31];
magmaFloatComplex G236 = t * s_targetCoefs[31];
magmaFloatComplex G237 = G235 + G236;
magmaFloatComplex G238 = C0 * G237;
magmaFloatComplex G239 = C3 * G238;
magmaFloatComplex G240 = G234 + G239;
magmaFloatComplex G241 = G218 * G16;
magmaFloatComplex G242 = G220 * G46;
magmaFloatComplex G243 = G241 + G242;
magmaFloatComplex G244 = G223 * G68;
magmaFloatComplex G245 = G243 + G244;
magmaFloatComplex G246 = G226 * G90;
magmaFloatComplex G247 = G245 + G246;
magmaFloatComplex G248 = C0 * s_track[5];
magmaFloatComplex G249 = G247 + G248;
magmaFloatComplex G250 = s_track[7] * G125;
magmaFloatComplex G251 = G249 + G250;
magmaFloatComplex G252 = s_track[8] * G132;
magmaFloatComplex G253 = G251 + G252;
magmaFloatComplex G254 = G1 * s_startCoefs[32];
magmaFloatComplex G255 = t * s_targetCoefs[32];
magmaFloatComplex G256 = G254 + G255;
magmaFloatComplex G257 = C0 * G256;
magmaFloatComplex G258 = C3 * G257;
magmaFloatComplex G259 = G253 + G258;
magmaFloatComplex G260 = s_track[0] * s_track[2];
magmaFloatComplex G261 = G260 * G4;
magmaFloatComplex G262 = s_track[0] * s_track[3];
magmaFloatComplex G263 = G262 * G8;
magmaFloatComplex G264 = G261 + G263;
magmaFloatComplex G265 = s_track[0] * G23;
magmaFloatComplex G266 = G264 + G265;
magmaFloatComplex G267 = s_track[1] * s_track[2];
magmaFloatComplex G268 = G267 * G34;
magmaFloatComplex G269 = G266 + G268;
magmaFloatComplex G270 = s_track[1] * s_track[3];
magmaFloatComplex G271 = G270 * G38;
magmaFloatComplex G272 = G269 + G271;
magmaFloatComplex G273 = s_track[1] * G53;
magmaFloatComplex G274 = G272 + G273;
magmaFloatComplex G275 = s_track[2] * G75;
magmaFloatComplex G276 = G274 + G275;
magmaFloatComplex G277 = s_track[3] * G97;
magmaFloatComplex G278 = G276 + G277;
magmaFloatComplex G279 = G1 * s_startCoefs[8];
magmaFloatComplex G280 = t * s_targetCoefs[8];
magmaFloatComplex G281 = G279 + G280;
magmaFloatComplex G282 = G278 + G281;
magmaFloatComplex G283 = s_track[2] * s_track[4];
magmaFloatComplex G284 = G283 * G64;
magmaFloatComplex G285 = s_track[2] * s_track[5];
magmaFloatComplex G286 = G285 * G68;
magmaFloatComplex G287 = G284 + G286;
magmaFloatComplex G288 = s_track[2] * G82;
magmaFloatComplex G289 = G287 + G288;
magmaFloatComplex G290 = s_track[3] * s_track[4];
magmaFloatComplex G291 = G290 * G86;
magmaFloatComplex G292 = G289 + G291;
magmaFloatComplex G293 = s_track[3] * s_track[5];
magmaFloatComplex G294 = G293 * G90;
magmaFloatComplex G295 = G292 + G294;
magmaFloatComplex G296 = s_track[3] * G104;
magmaFloatComplex G297 = G295 + G296;
magmaFloatComplex G298 = s_track[4] * G111;
magmaFloatComplex G299 = G297 + G298;
magmaFloatComplex G300 = s_track[5] * G125;
magmaFloatComplex G301 = G299 + G300;
magmaFloatComplex G302 = G1 * s_startCoefs[17];
magmaFloatComplex G303 = t * s_targetCoefs[17];
magmaFloatComplex G304 = G302 + G303;
magmaFloatComplex G305 = G301 + G304;
magmaFloatComplex G306 = s_track[0] * s_track[4];
magmaFloatComplex G307 = G306 * G12;
magmaFloatComplex G308 = s_track[0] * s_track[5];
magmaFloatComplex G309 = G308 * G16;
magmaFloatComplex G310 = G307 + G309;
magmaFloatComplex G311 = s_track[0] * G30;
magmaFloatComplex G312 = G310 + G311;
magmaFloatComplex G313 = s_track[1] * s_track[4];
magmaFloatComplex G314 = G313 * G42;
magmaFloatComplex G315 = G312 + G314;
magmaFloatComplex G316 = s_track[1] * s_track[5];
magmaFloatComplex G317 = G316 * G46;
magmaFloatComplex G318 = G315 + G317;
magmaFloatComplex G319 = s_track[1] * G60;
magmaFloatComplex G320 = G318 + G319;
magmaFloatComplex G321 = s_track[4] * G118;
magmaFloatComplex G322 = G320 + G321;
magmaFloatComplex G323 = s_track[5] * G132;
magmaFloatComplex G324 = G322 + G323;
magmaFloatComplex G325 = G1 * s_startCoefs[26];
magmaFloatComplex G326 = t * s_targetCoefs[26];
magmaFloatComplex G327 = G325 + G326;
magmaFloatComplex G328 = G324 + G327;

r_cgesvA[0] = C0;
r_cgesvA[1] = C1;
r_cgesvA[2] = G5;
r_cgesvA[3] = G9;
r_cgesvA[4] = G13;
r_cgesvA[5] = G17;
r_cgesvA[6] = G24;
r_cgesvA[7] = C1;
r_cgesvA[8] = G31;
r_cgesvB[0] =G156;

r_cgesvA[9] = C1;
r_cgesvA[10] = C0;
r_cgesvA[11] = G35;
r_cgesvA[12] = G39;
r_cgesvA[13] = G43;
r_cgesvA[14] = G47;
r_cgesvA[15] = G54;
r_cgesvA[16] = C1;
r_cgesvA[17] = G61;
r_cgesvB[1] =G175;

r_cgesvA[18] = G5;
r_cgesvA[19] = G35;
r_cgesvA[20] = C0;
r_cgesvA[21] = C1;
r_cgesvA[22] = G65;
r_cgesvA[23] = G69;
r_cgesvA[24] = G76;
r_cgesvA[25] = G83;
r_cgesvA[26] = C1;
r_cgesvB[2] =G198;

r_cgesvA[27] = G9;
r_cgesvA[28] = G39;
r_cgesvA[29] = C1;
r_cgesvA[30] = C0;
r_cgesvA[31] = G87;
r_cgesvA[32] = G91;
r_cgesvA[33] = G98;
r_cgesvA[34] = G105;
r_cgesvA[35] = C1;
r_cgesvB[3] =G217;

r_cgesvA[36] = G13;
r_cgesvA[37] = G43;
r_cgesvA[38] = G65;
r_cgesvA[39] = G87;
r_cgesvA[40] = C0;
r_cgesvA[41] = C1;
r_cgesvA[42] = C1;
r_cgesvA[43] = G112;
r_cgesvA[44] = G119;
r_cgesvB[4] =G240;

r_cgesvA[45] = G17;
r_cgesvA[46] = G47;
r_cgesvA[47] = G69;
r_cgesvA[48] = G91;
r_cgesvA[49] = C1;
r_cgesvA[50] = C0;
r_cgesvA[51] = C1;
r_cgesvA[52] = G126;
r_cgesvA[53] = G133;
r_cgesvB[5] =G259;

r_cgesvA[54] = G24;
r_cgesvA[55] = G54;
r_cgesvA[56] = G76;
r_cgesvA[57] = G98;
r_cgesvA[58] = C1;
r_cgesvA[59] = C1;
r_cgesvA[60] = C1;
r_cgesvA[61] = C1;
r_cgesvA[62] = C1;
r_cgesvB[6] =G282;

r_cgesvA[63] = C1;
r_cgesvA[64] = C1;
r_cgesvA[65] = G83;
r_cgesvA[66] = G105;
r_cgesvA[67] = G112;
r_cgesvA[68] = G126;
r_cgesvA[69] = C1;
r_cgesvA[70] = C1;
r_cgesvA[71] = C1;
r_cgesvB[7] =G305;

r_cgesvA[72] = G31;
r_cgesvA[73] = G61;
r_cgesvA[74] = C1;
r_cgesvA[75] = C1;
r_cgesvA[76] = G119;
r_cgesvA[77] = G133;
r_cgesvA[78] = C1;
r_cgesvA[79] = C1;
r_cgesvA[80] = C1;
r_cgesvB[8] =G328;

  }
}

#endif
