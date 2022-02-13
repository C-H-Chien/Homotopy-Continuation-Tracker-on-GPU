#ifndef cpu_eval_HxH_cyclic7_h
#define cpu_eval_HxH_cyclic7_h
// ============================================================================
// partial derivative evaluations of the alea6 problem for cpu HC computation
//
// Modifications
//    Chien  21-06-30:   Originally created
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
  void cpu_eval_HxH_cyclic7(
      magma_int_t s, float t, int N, magmaFloatComplex* s_track,
      const magmaFloatComplex &C0, const magmaFloatComplex &C1, const magmaFloatComplex &C2,
      magmaFloatComplex* s_startCoefs, magmaFloatComplex* s_targetCoefs,
      magmaFloatComplex* r_cgesvA, magmaFloatComplex* r_cgesvB )
  {
    magmaFloatComplex G1 = C0 - t;
magmaFloatComplex G2 = G1 * s_startCoefs[0];
magmaFloatComplex G3 = t * s_targetCoefs[0];
magmaFloatComplex G4 = G2 + G3;
magmaFloatComplex G5 = s_track[1] * G4;
magmaFloatComplex G6 = s_track[6] * G4;
magmaFloatComplex G7 = G5 + G6;
magmaFloatComplex G8 = s_track[2] * G5;
magmaFloatComplex G9 = s_track[6] * G5;
magmaFloatComplex G10 = G8 + G9;
magmaFloatComplex G11 = s_track[5] * G4;
magmaFloatComplex G12 = s_track[6] * G11;
magmaFloatComplex G13 = G10 + G12;
magmaFloatComplex G14 = s_track[3] * G8;
magmaFloatComplex G15 = s_track[6] * G8;
magmaFloatComplex G16 = G14 + G15;
magmaFloatComplex G17 = s_track[5] * G5;
magmaFloatComplex G18 = s_track[6] * G17;
magmaFloatComplex G19 = G16 + G18;
magmaFloatComplex G20 = s_track[4] * G4;
magmaFloatComplex G21 = s_track[5] * G20;
magmaFloatComplex G22 = s_track[6] * G21;
magmaFloatComplex G23 = G19 + G22;
magmaFloatComplex G24 = s_track[4] * G14;
magmaFloatComplex G25 = s_track[6] * G14;
magmaFloatComplex G26 = G24 + G25;
magmaFloatComplex G27 = s_track[5] * G8;
magmaFloatComplex G28 = s_track[6] * G27;
magmaFloatComplex G29 = G26 + G28;
magmaFloatComplex G30 = s_track[4] * G5;
magmaFloatComplex G31 = s_track[5] * G30;
magmaFloatComplex G32 = s_track[6] * G31;
magmaFloatComplex G33 = G29 + G32;
magmaFloatComplex G34 = s_track[3] * G4;
magmaFloatComplex G35 = s_track[4] * G34;
magmaFloatComplex G36 = s_track[5] * G35;
magmaFloatComplex G37 = s_track[6] * G36;
magmaFloatComplex G38 = G33 + G37;
magmaFloatComplex G39 = s_track[5] * G24;
magmaFloatComplex G40 = s_track[6] * G24;
magmaFloatComplex G41 = G39 + G40;
magmaFloatComplex G42 = s_track[5] * G14;
magmaFloatComplex G43 = s_track[6] * G42;
magmaFloatComplex G44 = G41 + G43;
magmaFloatComplex G45 = s_track[4] * G8;
magmaFloatComplex G46 = s_track[5] * G45;
magmaFloatComplex G47 = s_track[6] * G46;
magmaFloatComplex G48 = G44 + G47;
magmaFloatComplex G49 = s_track[3] * G5;
magmaFloatComplex G50 = s_track[4] * G49;
magmaFloatComplex G51 = s_track[5] * G50;
magmaFloatComplex G52 = s_track[6] * G51;
magmaFloatComplex G53 = G48 + G52;
magmaFloatComplex G54 = s_track[2] * G4;
magmaFloatComplex G55 = s_track[3] * G54;
magmaFloatComplex G56 = s_track[4] * G55;
magmaFloatComplex G57 = s_track[5] * G56;
magmaFloatComplex G58 = s_track[6] * G57;
magmaFloatComplex G59 = G53 + G58;
magmaFloatComplex G60 = s_track[6] * G39;
magmaFloatComplex G61 = G4 * s_track[0];
magmaFloatComplex G62 = G61 + G54;
magmaFloatComplex G63 = s_track[2] * G61;
magmaFloatComplex G64 = s_track[6] * G61;
magmaFloatComplex G65 = G63 + G64;
magmaFloatComplex G66 = G65 + G55;
magmaFloatComplex G67 = s_track[3] * G63;
magmaFloatComplex G68 = s_track[6] * G63;
magmaFloatComplex G69 = G67 + G68;
magmaFloatComplex G70 = s_track[5] * G61;
magmaFloatComplex G71 = s_track[6] * G70;
magmaFloatComplex G72 = G69 + G71;
magmaFloatComplex G73 = G72 + G56;
magmaFloatComplex G74 = s_track[4] * G67;
magmaFloatComplex G75 = s_track[6] * G67;
magmaFloatComplex G76 = G74 + G75;
magmaFloatComplex G77 = s_track[5] * G63;
magmaFloatComplex G78 = s_track[6] * G77;
magmaFloatComplex G79 = G76 + G78;
magmaFloatComplex G80 = s_track[4] * G61;
magmaFloatComplex G81 = s_track[5] * G80;
magmaFloatComplex G82 = s_track[6] * G81;
magmaFloatComplex G83 = G79 + G82;
magmaFloatComplex G84 = G83 + G57;
magmaFloatComplex G85 = s_track[5] * G74;
magmaFloatComplex G86 = s_track[6] * G74;
magmaFloatComplex G87 = G85 + G86;
magmaFloatComplex G88 = s_track[5] * G67;
magmaFloatComplex G89 = s_track[6] * G88;
magmaFloatComplex G90 = G87 + G89;
magmaFloatComplex G91 = s_track[4] * G63;
magmaFloatComplex G92 = s_track[5] * G91;
magmaFloatComplex G93 = s_track[6] * G92;
magmaFloatComplex G94 = G90 + G93;
magmaFloatComplex G95 = s_track[3] * G61;
magmaFloatComplex G96 = s_track[4] * G95;
magmaFloatComplex G97 = s_track[5] * G96;
magmaFloatComplex G98 = s_track[6] * G97;
magmaFloatComplex G99 = G94 + G98;
magmaFloatComplex G100 = G99 + G58;
magmaFloatComplex G101 = s_track[6] * G85;
magmaFloatComplex G102 = G4 * s_track[1];
magmaFloatComplex G103 = G102 + G34;
magmaFloatComplex G104 = G61 * s_track[1];
magmaFloatComplex G105 = s_track[3] * G102;
magmaFloatComplex G106 = G104 + G105;
magmaFloatComplex G107 = G106 + G35;
magmaFloatComplex G108 = s_track[3] * G104;
magmaFloatComplex G109 = s_track[6] * G104;
magmaFloatComplex G110 = G108 + G109;
magmaFloatComplex G111 = s_track[4] * G105;
magmaFloatComplex G112 = G110 + G111;
magmaFloatComplex G113 = G112 + G36;
magmaFloatComplex G114 = s_track[4] * G108;
magmaFloatComplex G115 = s_track[6] * G108;
magmaFloatComplex G116 = G114 + G115;
magmaFloatComplex G117 = s_track[5] * G104;
magmaFloatComplex G118 = s_track[6] * G117;
magmaFloatComplex G119 = G116 + G118;
magmaFloatComplex G120 = s_track[5] * G111;
magmaFloatComplex G121 = G119 + G120;
magmaFloatComplex G122 = G121 + G37;
magmaFloatComplex G123 = s_track[5] * G114;
magmaFloatComplex G124 = s_track[6] * G114;
magmaFloatComplex G125 = G123 + G124;
magmaFloatComplex G126 = s_track[5] * G108;
magmaFloatComplex G127 = s_track[6] * G126;
magmaFloatComplex G128 = G125 + G127;
magmaFloatComplex G129 = s_track[4] * G104;
magmaFloatComplex G130 = s_track[5] * G129;
magmaFloatComplex G131 = s_track[6] * G130;
magmaFloatComplex G132 = G128 + G131;
magmaFloatComplex G133 = G132 + G98;
magmaFloatComplex G134 = s_track[6] * G120;
magmaFloatComplex G135 = G133 + G134;
magmaFloatComplex G136 = s_track[6] * G123;
magmaFloatComplex G137 = G4 * s_track[2];
magmaFloatComplex G138 = G137 + G20;
magmaFloatComplex G139 = G102 * s_track[2];
magmaFloatComplex G140 = s_track[4] * G137;
magmaFloatComplex G141 = G139 + G140;
magmaFloatComplex G142 = G141 + G21;
magmaFloatComplex G143 = G104 * s_track[2];
magmaFloatComplex G144 = s_track[4] * G139;
magmaFloatComplex G145 = G143 + G144;
magmaFloatComplex G146 = s_track[5] * G140;
magmaFloatComplex G147 = G145 + G146;
magmaFloatComplex G148 = G147 + G22;
magmaFloatComplex G149 = s_track[4] * G143;
magmaFloatComplex G150 = s_track[6] * G143;
magmaFloatComplex G151 = G149 + G150;
magmaFloatComplex G152 = G151 + G82;
magmaFloatComplex G153 = s_track[5] * G144;
magmaFloatComplex G154 = G152 + G153;
magmaFloatComplex G155 = s_track[6] * G146;
magmaFloatComplex G156 = G154 + G155;
magmaFloatComplex G157 = s_track[5] * G149;
magmaFloatComplex G158 = s_track[6] * G149;
magmaFloatComplex G159 = G157 + G158;
magmaFloatComplex G160 = s_track[5] * G143;
magmaFloatComplex G161 = s_track[6] * G160;
magmaFloatComplex G162 = G159 + G161;
magmaFloatComplex G163 = G162 + G131;
magmaFloatComplex G164 = G61 * s_track[2];
magmaFloatComplex G165 = s_track[4] * G164;
magmaFloatComplex G166 = s_track[5] * G165;
magmaFloatComplex G167 = s_track[6] * G166;
magmaFloatComplex G168 = G163 + G167;
magmaFloatComplex G169 = s_track[6] * G153;
magmaFloatComplex G170 = G168 + G169;
magmaFloatComplex G171 = s_track[6] * G157;
magmaFloatComplex G172 = G4 * s_track[3];
magmaFloatComplex G173 = G172 + G11;
magmaFloatComplex G174 = G137 * s_track[3];
magmaFloatComplex G175 = s_track[5] * G172;
magmaFloatComplex G176 = G174 + G175;
magmaFloatComplex G177 = G176 + G12;
magmaFloatComplex G178 = G139 * s_track[3];
magmaFloatComplex G179 = G71 + G178;
magmaFloatComplex G180 = s_track[5] * G174;
magmaFloatComplex G181 = G179 + G180;
magmaFloatComplex G182 = s_track[6] * G175;
magmaFloatComplex G183 = G181 + G182;
magmaFloatComplex G184 = G143 * s_track[3];
magmaFloatComplex G185 = G184 + G118;
magmaFloatComplex G186 = G61 * s_track[3];
magmaFloatComplex G187 = s_track[5] * G186;
magmaFloatComplex G188 = s_track[6] * G187;
magmaFloatComplex G189 = G185 + G188;
magmaFloatComplex G190 = s_track[5] * G178;
magmaFloatComplex G191 = G189 + G190;
magmaFloatComplex G192 = s_track[6] * G180;
magmaFloatComplex G193 = G191 + G192;
magmaFloatComplex G194 = s_track[5] * G184;
magmaFloatComplex G195 = s_track[6] * G184;
magmaFloatComplex G196 = G194 + G195;
magmaFloatComplex G197 = G196 + G161;
magmaFloatComplex G198 = G104 * s_track[3];
magmaFloatComplex G199 = s_track[5] * G198;
magmaFloatComplex G200 = s_track[6] * G199;
magmaFloatComplex G201 = G197 + G200;
magmaFloatComplex G202 = G164 * s_track[3];
magmaFloatComplex G203 = s_track[5] * G202;
magmaFloatComplex G204 = s_track[6] * G203;
magmaFloatComplex G205 = G201 + G204;
magmaFloatComplex G206 = s_track[6] * G190;
magmaFloatComplex G207 = G205 + G206;
magmaFloatComplex G208 = s_track[6] * G194;
magmaFloatComplex G209 = G4 * s_track[4];
magmaFloatComplex G210 = G209 + G6;
magmaFloatComplex G211 = G172 * s_track[4];
magmaFloatComplex G212 = G64 + G211;
magmaFloatComplex G213 = s_track[6] * G209;
magmaFloatComplex G214 = G212 + G213;
magmaFloatComplex G215 = G61 * s_track[4];
magmaFloatComplex G216 = s_track[6] * G215;
magmaFloatComplex G217 = G109 + G216;
magmaFloatComplex G218 = G174 * s_track[4];
magmaFloatComplex G219 = G217 + G218;
magmaFloatComplex G220 = s_track[6] * G211;
magmaFloatComplex G221 = G219 + G220;
magmaFloatComplex G222 = G104 * s_track[4];
magmaFloatComplex G223 = s_track[6] * G222;
magmaFloatComplex G224 = G150 + G223;
magmaFloatComplex G225 = G186 * s_track[4];
magmaFloatComplex G226 = s_track[6] * G225;
magmaFloatComplex G227 = G224 + G226;
magmaFloatComplex G228 = G178 * s_track[4];
magmaFloatComplex G229 = G227 + G228;
magmaFloatComplex G230 = s_track[6] * G218;
magmaFloatComplex G231 = G229 + G230;
magmaFloatComplex G232 = G184 * s_track[4];
magmaFloatComplex G233 = G232 + G195;
magmaFloatComplex G234 = G143 * s_track[4];
magmaFloatComplex G235 = s_track[6] * G234;
magmaFloatComplex G236 = G233 + G235;
magmaFloatComplex G237 = G198 * s_track[4];
magmaFloatComplex G238 = s_track[6] * G237;
magmaFloatComplex G239 = G236 + G238;
magmaFloatComplex G240 = G202 * s_track[4];
magmaFloatComplex G241 = s_track[6] * G240;
magmaFloatComplex G242 = G239 + G241;
magmaFloatComplex G243 = s_track[6] * G228;
magmaFloatComplex G244 = G242 + G243;
magmaFloatComplex G245 = s_track[6] * G232;
magmaFloatComplex G246 = G4 * s_track[5];
magmaFloatComplex G247 = G61 + G246;
magmaFloatComplex G248 = G61 * s_track[5];
magmaFloatComplex G249 = G104 + G248;
magmaFloatComplex G250 = G209 * s_track[5];
magmaFloatComplex G251 = G249 + G250;
magmaFloatComplex G252 = G104 * s_track[5];
magmaFloatComplex G253 = G143 + G252;
magmaFloatComplex G254 = G215 * s_track[5];
magmaFloatComplex G255 = G253 + G254;
magmaFloatComplex G256 = G211 * s_track[5];
magmaFloatComplex G257 = G255 + G256;
magmaFloatComplex G258 = G143 * s_track[5];
magmaFloatComplex G259 = G184 + G258;
magmaFloatComplex G260 = G222 * s_track[5];
magmaFloatComplex G261 = G259 + G260;
magmaFloatComplex G262 = G225 * s_track[5];
magmaFloatComplex G263 = G261 + G262;
magmaFloatComplex G264 = G218 * s_track[5];
magmaFloatComplex G265 = G263 + G264;
magmaFloatComplex G266 = G184 * s_track[5];
magmaFloatComplex G267 = G232 + G266;
magmaFloatComplex G268 = G234 * s_track[5];
magmaFloatComplex G269 = G267 + G268;
magmaFloatComplex G270 = G237 * s_track[5];
magmaFloatComplex G271 = G269 + G270;
magmaFloatComplex G272 = G240 * s_track[5];
magmaFloatComplex G273 = G271 + G272;
magmaFloatComplex G274 = G228 * s_track[5];
magmaFloatComplex G275 = G273 + G274;
magmaFloatComplex G276 = G232 * s_track[5];
magmaFloatComplex G277 = G61 + G102;
magmaFloatComplex G278 = G277 + G137;
magmaFloatComplex G279 = G278 + G172;
magmaFloatComplex G280 = G279 + G209;
magmaFloatComplex G281 = G280 + G246;
magmaFloatComplex G282 = G4 * s_track[6];
magmaFloatComplex G283 = G281 + G282;
magmaFloatComplex G284 = G61 * s_track[6];
magmaFloatComplex G285 = G104 + G284;
magmaFloatComplex G286 = G285 + G139;
magmaFloatComplex G287 = G286 + G174;
magmaFloatComplex G288 = G287 + G211;
magmaFloatComplex G289 = G288 + G250;
magmaFloatComplex G290 = G246 * s_track[6];
magmaFloatComplex G291 = G289 + G290;
magmaFloatComplex G292 = G104 * s_track[6];
magmaFloatComplex G293 = G143 + G292;
magmaFloatComplex G294 = G248 * s_track[6];
magmaFloatComplex G295 = G293 + G294;
magmaFloatComplex G296 = G295 + G178;
magmaFloatComplex G297 = G296 + G218;
magmaFloatComplex G298 = G297 + G256;
magmaFloatComplex G299 = G250 * s_track[6];
magmaFloatComplex G300 = G298 + G299;
magmaFloatComplex G301 = G143 * s_track[6];
magmaFloatComplex G302 = G184 + G301;
magmaFloatComplex G303 = G252 * s_track[6];
magmaFloatComplex G304 = G302 + G303;
magmaFloatComplex G305 = G254 * s_track[6];
magmaFloatComplex G306 = G304 + G305;
magmaFloatComplex G307 = G306 + G228;
magmaFloatComplex G308 = G307 + G264;
magmaFloatComplex G309 = G256 * s_track[6];
magmaFloatComplex G310 = G308 + G309;
magmaFloatComplex G311 = G184 * s_track[6];
magmaFloatComplex G312 = G232 + G311;
magmaFloatComplex G313 = G258 * s_track[6];
magmaFloatComplex G314 = G312 + G313;
magmaFloatComplex G315 = G260 * s_track[6];
magmaFloatComplex G316 = G314 + G315;
magmaFloatComplex G317 = G262 * s_track[6];
magmaFloatComplex G318 = G316 + G317;
magmaFloatComplex G319 = G318 + G274;
magmaFloatComplex G320 = G264 * s_track[6];
magmaFloatComplex G321 = G319 + G320;
magmaFloatComplex G322 = G232 * s_track[6];
magmaFloatComplex G323 = G276 + G322;
magmaFloatComplex G324 = G266 * s_track[6];
magmaFloatComplex G325 = G323 + G324;
magmaFloatComplex G326 = G268 * s_track[6];
magmaFloatComplex G327 = G325 + G326;
magmaFloatComplex G328 = G270 * s_track[6];
magmaFloatComplex G329 = G327 + G328;
magmaFloatComplex G330 = G272 * s_track[6];
magmaFloatComplex G331 = G329 + G330;
magmaFloatComplex G332 = G274 * s_track[6];
magmaFloatComplex G333 = G331 + G332;
magmaFloatComplex G334 = G276 * s_track[6];
magmaFloatComplex G335 = G1 * s_startCoefs[1];
magmaFloatComplex G336 = t * s_targetCoefs[1];
magmaFloatComplex G337 = G335 + G336;
magmaFloatComplex G338 = G334 + G337;

r_cgesvA[0] = G4;
r_cgesvA[1] = G7;
r_cgesvA[2] = G13;
r_cgesvA[3] = G23;
r_cgesvA[4] = G38;
r_cgesvA[5] = G59;
r_cgesvA[6] = G60;
r_cgesvB[0] =G283;

r_cgesvA[7] = G4;
r_cgesvA[8] = G62;
r_cgesvA[9] = G66;
r_cgesvA[10] = G73;
r_cgesvA[11] = G84;
r_cgesvA[12] = G100;
r_cgesvA[13] = G101;
r_cgesvB[1] =G291;

r_cgesvA[14] = G4;
r_cgesvA[15] = G103;
r_cgesvA[16] = G107;
r_cgesvA[17] = G113;
r_cgesvA[18] = G122;
r_cgesvA[19] = G135;
r_cgesvA[20] = G136;
r_cgesvB[2] =G300;

r_cgesvA[21] = G4;
r_cgesvA[22] = G138;
r_cgesvA[23] = G142;
r_cgesvA[24] = G148;
r_cgesvA[25] = G156;
r_cgesvA[26] = G170;
r_cgesvA[27] = G171;
r_cgesvB[3] =G310;

r_cgesvA[28] = G4;
r_cgesvA[29] = G173;
r_cgesvA[30] = G177;
r_cgesvA[31] = G183;
r_cgesvA[32] = G193;
r_cgesvA[33] = G207;
r_cgesvA[34] = G208;
r_cgesvB[4] =G321;

r_cgesvA[35] = G4;
r_cgesvA[36] = G210;
r_cgesvA[37] = G214;
r_cgesvA[38] = G221;
r_cgesvA[39] = G231;
r_cgesvA[40] = G244;
r_cgesvA[41] = G245;
r_cgesvB[5] =G333;

r_cgesvA[42] = G4;
r_cgesvA[43] = G247;
r_cgesvA[44] = G251;
r_cgesvA[45] = G257;
r_cgesvA[46] = G265;
r_cgesvA[47] = G275;
r_cgesvA[48] = G276;
r_cgesvB[6] =G338;

  }
}

#endif
