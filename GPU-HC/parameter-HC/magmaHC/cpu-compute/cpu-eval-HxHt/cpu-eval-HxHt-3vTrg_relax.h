#ifndef cpu_eval_HxHt_3vTrg_relax_h
#define cpu_eval_HxHt_3vTrg_relax_h
// ============================================================================
// partial derivative evaluations of the relaxed 3-view triangulation problem
//
// Modifications
//    Chien  22-01-02:   Originally created
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
  void cpu_eval_HxHt_3vTrg_relax(
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
magmaFloatComplex G6 = G1 * s_startCoefs[3];
magmaFloatComplex G7 = t * s_targetCoefs[3];
magmaFloatComplex G8 = G6 + G7;
magmaFloatComplex G9 = G8 * s_track[6];
magmaFloatComplex G10 = G4 * s_track[2];
magmaFloatComplex G11 = G8 * s_track[3];
magmaFloatComplex G12 = G10 + G11;
magmaFloatComplex G13 = G1 * s_startCoefs[6];
magmaFloatComplex G14 = t * s_targetCoefs[6];
magmaFloatComplex G15 = G13 + G14;
magmaFloatComplex G16 = G12 + G15;
magmaFloatComplex G17 = G1 * s_startCoefs[1];
magmaFloatComplex G18 = t * s_targetCoefs[1];
magmaFloatComplex G19 = G17 + G18;
magmaFloatComplex G20 = G19 * s_track[6];
magmaFloatComplex G21 = G1 * s_startCoefs[4];
magmaFloatComplex G22 = t * s_targetCoefs[4];
magmaFloatComplex G23 = G21 + G22;
magmaFloatComplex G24 = G23 * s_track[6];
magmaFloatComplex G25 = G19 * s_track[2];
magmaFloatComplex G26 = G23 * s_track[3];
magmaFloatComplex G27 = G25 + G26;
magmaFloatComplex G28 = G1 * s_startCoefs[7];
magmaFloatComplex G29 = t * s_targetCoefs[7];
magmaFloatComplex G30 = G28 + G29;
magmaFloatComplex G31 = G27 + G30;
magmaFloatComplex G32 = G1 * s_startCoefs[9];
magmaFloatComplex G33 = t * s_targetCoefs[9];
magmaFloatComplex G34 = G32 + G33;
magmaFloatComplex G35 = G34 * s_track[7];
magmaFloatComplex G36 = G1 * s_startCoefs[12];
magmaFloatComplex G37 = t * s_targetCoefs[12];
magmaFloatComplex G38 = G36 + G37;
magmaFloatComplex G39 = G38 * s_track[7];
magmaFloatComplex G40 = G4 * s_track[0];
magmaFloatComplex G41 = G19 * s_track[1];
magmaFloatComplex G42 = G40 + G41;
magmaFloatComplex G43 = G1 * s_startCoefs[2];
magmaFloatComplex G44 = t * s_targetCoefs[2];
magmaFloatComplex G45 = G43 + G44;
magmaFloatComplex G46 = G42 + G45;
magmaFloatComplex G47 = G34 * s_track[4];
magmaFloatComplex G48 = G38 * s_track[5];
magmaFloatComplex G49 = G47 + G48;
magmaFloatComplex G50 = G1 * s_startCoefs[15];
magmaFloatComplex G51 = t * s_targetCoefs[15];
magmaFloatComplex G52 = G50 + G51;
magmaFloatComplex G53 = G49 + G52;
magmaFloatComplex G54 = G1 * s_startCoefs[10];
magmaFloatComplex G55 = t * s_targetCoefs[10];
magmaFloatComplex G56 = G54 + G55;
magmaFloatComplex G57 = G56 * s_track[7];
magmaFloatComplex G58 = G1 * s_startCoefs[13];
magmaFloatComplex G59 = t * s_targetCoefs[13];
magmaFloatComplex G60 = G58 + G59;
magmaFloatComplex G61 = G60 * s_track[7];
magmaFloatComplex G62 = G8 * s_track[0];
magmaFloatComplex G63 = G23 * s_track[1];
magmaFloatComplex G64 = G62 + G63;
magmaFloatComplex G65 = G1 * s_startCoefs[5];
magmaFloatComplex G66 = t * s_targetCoefs[5];
magmaFloatComplex G67 = G65 + G66;
magmaFloatComplex G68 = G64 + G67;
magmaFloatComplex G69 = G56 * s_track[4];
magmaFloatComplex G70 = G60 * s_track[5];
magmaFloatComplex G71 = G69 + G70;
magmaFloatComplex G72 = G1 * s_startCoefs[16];
magmaFloatComplex G73 = t * s_targetCoefs[16];
magmaFloatComplex G74 = G72 + G73;
magmaFloatComplex G75 = G71 + G74;
magmaFloatComplex G76 = G34 * s_track[2];
magmaFloatComplex G77 = G56 * s_track[3];
magmaFloatComplex G78 = G76 + G77;
magmaFloatComplex G79 = G1 * s_startCoefs[11];
magmaFloatComplex G80 = t * s_targetCoefs[11];
magmaFloatComplex G81 = G79 + G80;
magmaFloatComplex G82 = G78 + G81;
magmaFloatComplex G83 = G38 * s_track[2];
magmaFloatComplex G84 = G60 * s_track[3];
magmaFloatComplex G85 = G83 + G84;
magmaFloatComplex G86 = G1 * s_startCoefs[14];
magmaFloatComplex G87 = t * s_targetCoefs[14];
magmaFloatComplex G88 = G86 + G87;
magmaFloatComplex G89 = G85 + G88;
magmaFloatComplex G90 = s_track[2] * s_track[6];
magmaFloatComplex G91 = C3 * s_startCoefs[0];
magmaFloatComplex G92 = G91 + s_targetCoefs[0];
magmaFloatComplex G93 = G90 * G92;
magmaFloatComplex G94 = s_track[3] * s_track[6];
magmaFloatComplex G95 = C3 * s_startCoefs[3];
magmaFloatComplex G96 = G95 + s_targetCoefs[3];
magmaFloatComplex G97 = G94 * G96;
magmaFloatComplex G98 = G93 + G97;
magmaFloatComplex G99 = C3 * s_startCoefs[6];
magmaFloatComplex G100 = G99 + s_targetCoefs[6];
magmaFloatComplex G101 = s_track[6] * G100;
magmaFloatComplex G102 = G98 + G101;
magmaFloatComplex G103 = C3 * s_startCoefs[18];
magmaFloatComplex G104 = G103 + s_targetCoefs[18];
magmaFloatComplex G105 = C0 * G104;
magmaFloatComplex G106 = C3 * G105;
magmaFloatComplex G107 = G102 + G106;
magmaFloatComplex G108 = C3 * s_startCoefs[1];
magmaFloatComplex G109 = G108 + s_targetCoefs[1];
magmaFloatComplex G110 = G90 * G109;
magmaFloatComplex G111 = C3 * s_startCoefs[4];
magmaFloatComplex G112 = G111 + s_targetCoefs[4];
magmaFloatComplex G113 = G94 * G112;
magmaFloatComplex G114 = G110 + G113;
magmaFloatComplex G115 = C3 * s_startCoefs[7];
magmaFloatComplex G116 = G115 + s_targetCoefs[7];
magmaFloatComplex G117 = s_track[6] * G116;
magmaFloatComplex G118 = G114 + G117;
magmaFloatComplex G119 = C3 * s_startCoefs[19];
magmaFloatComplex G120 = G119 + s_targetCoefs[19];
magmaFloatComplex G121 = C0 * G120;
magmaFloatComplex G122 = C3 * G121;
magmaFloatComplex G123 = G118 + G122;
magmaFloatComplex G124 = s_track[0] * s_track[6];
magmaFloatComplex G125 = G124 * G92;
magmaFloatComplex G126 = s_track[1] * s_track[6];
magmaFloatComplex G127 = G126 * G109;
magmaFloatComplex G128 = G125 + G127;
magmaFloatComplex G129 = s_track[4] * s_track[7];
magmaFloatComplex G130 = C3 * s_startCoefs[9];
magmaFloatComplex G131 = G130 + s_targetCoefs[9];
magmaFloatComplex G132 = G129 * G131;
magmaFloatComplex G133 = G128 + G132;
magmaFloatComplex G134 = s_track[5] * s_track[7];
magmaFloatComplex G135 = C3 * s_startCoefs[12];
magmaFloatComplex G136 = G135 + s_targetCoefs[12];
magmaFloatComplex G137 = G134 * G136;
magmaFloatComplex G138 = G133 + G137;
magmaFloatComplex G139 = C3 * s_startCoefs[2];
magmaFloatComplex G140 = G139 + s_targetCoefs[2];
magmaFloatComplex G141 = s_track[6] * G140;
magmaFloatComplex G142 = G138 + G141;
magmaFloatComplex G143 = C3 * s_startCoefs[15];
magmaFloatComplex G144 = G143 + s_targetCoefs[15];
magmaFloatComplex G145 = s_track[7] * G144;
magmaFloatComplex G146 = G142 + G145;
magmaFloatComplex G147 = C3 * s_startCoefs[20];
magmaFloatComplex G148 = G147 + s_targetCoefs[20];
magmaFloatComplex G149 = C0 * G148;
magmaFloatComplex G150 = C3 * G149;
magmaFloatComplex G151 = G146 + G150;
magmaFloatComplex G152 = G124 * G96;
magmaFloatComplex G153 = G126 * G112;
magmaFloatComplex G154 = G152 + G153;
magmaFloatComplex G155 = C3 * s_startCoefs[10];
magmaFloatComplex G156 = G155 + s_targetCoefs[10];
magmaFloatComplex G157 = G129 * G156;
magmaFloatComplex G158 = G154 + G157;
magmaFloatComplex G159 = C3 * s_startCoefs[13];
magmaFloatComplex G160 = G159 + s_targetCoefs[13];
magmaFloatComplex G161 = G134 * G160;
magmaFloatComplex G162 = G158 + G161;
magmaFloatComplex G163 = C3 * s_startCoefs[5];
magmaFloatComplex G164 = G163 + s_targetCoefs[5];
magmaFloatComplex G165 = s_track[6] * G164;
magmaFloatComplex G166 = G162 + G165;
magmaFloatComplex G167 = C3 * s_startCoefs[16];
magmaFloatComplex G168 = G167 + s_targetCoefs[16];
magmaFloatComplex G169 = s_track[7] * G168;
magmaFloatComplex G170 = G166 + G169;
magmaFloatComplex G171 = C3 * s_startCoefs[21];
magmaFloatComplex G172 = G171 + s_targetCoefs[21];
magmaFloatComplex G173 = C0 * G172;
magmaFloatComplex G174 = C3 * G173;
magmaFloatComplex G175 = G170 + G174;
magmaFloatComplex G176 = s_track[2] * s_track[7];
magmaFloatComplex G177 = G176 * G131;
magmaFloatComplex G178 = s_track[3] * s_track[7];
magmaFloatComplex G179 = G178 * G156;
magmaFloatComplex G180 = G177 + G179;
magmaFloatComplex G181 = C3 * s_startCoefs[11];
magmaFloatComplex G182 = G181 + s_targetCoefs[11];
magmaFloatComplex G183 = s_track[7] * G182;
magmaFloatComplex G184 = G180 + G183;
magmaFloatComplex G185 = C3 * s_startCoefs[22];
magmaFloatComplex G186 = G185 + s_targetCoefs[22];
magmaFloatComplex G187 = C0 * G186;
magmaFloatComplex G188 = C3 * G187;
magmaFloatComplex G189 = G184 + G188;
magmaFloatComplex G190 = G176 * G136;
magmaFloatComplex G191 = G178 * G160;
magmaFloatComplex G192 = G190 + G191;
magmaFloatComplex G193 = C3 * s_startCoefs[14];
magmaFloatComplex G194 = G193 + s_targetCoefs[14];
magmaFloatComplex G195 = s_track[7] * G194;
magmaFloatComplex G196 = G192 + G195;
magmaFloatComplex G197 = C3 * s_startCoefs[23];
magmaFloatComplex G198 = G197 + s_targetCoefs[23];
magmaFloatComplex G199 = C0 * G198;
magmaFloatComplex G200 = C3 * G199;
magmaFloatComplex G201 = G196 + G200;
magmaFloatComplex G202 = s_track[0] * s_track[2];
magmaFloatComplex G203 = G202 * G92;
magmaFloatComplex G204 = s_track[0] * s_track[3];
magmaFloatComplex G205 = G204 * G96;
magmaFloatComplex G206 = G203 + G205;
magmaFloatComplex G207 = s_track[0] * G100;
magmaFloatComplex G208 = G206 + G207;
magmaFloatComplex G209 = s_track[1] * s_track[2];
magmaFloatComplex G210 = G209 * G109;
magmaFloatComplex G211 = G208 + G210;
magmaFloatComplex G212 = s_track[1] * s_track[3];
magmaFloatComplex G213 = G212 * G112;
magmaFloatComplex G214 = G211 + G213;
magmaFloatComplex G215 = s_track[1] * G116;
magmaFloatComplex G216 = G214 + G215;
magmaFloatComplex G217 = s_track[2] * G140;
magmaFloatComplex G218 = G216 + G217;
magmaFloatComplex G219 = s_track[3] * G164;
magmaFloatComplex G220 = G218 + G219;
magmaFloatComplex G221 = C3 * s_startCoefs[8];
magmaFloatComplex G222 = G221 + s_targetCoefs[8];
magmaFloatComplex G223 = G220 + G222;
magmaFloatComplex G224 = s_track[2] * s_track[4];
magmaFloatComplex G225 = G224 * G131;
magmaFloatComplex G226 = s_track[2] * s_track[5];
magmaFloatComplex G227 = G226 * G136;
magmaFloatComplex G228 = G225 + G227;
magmaFloatComplex G229 = s_track[2] * G144;
magmaFloatComplex G230 = G228 + G229;
magmaFloatComplex G231 = s_track[3] * s_track[4];
magmaFloatComplex G232 = G231 * G156;
magmaFloatComplex G233 = G230 + G232;
magmaFloatComplex G234 = s_track[3] * s_track[5];
magmaFloatComplex G235 = G234 * G160;
magmaFloatComplex G236 = G233 + G235;
magmaFloatComplex G237 = s_track[3] * G168;
magmaFloatComplex G238 = G236 + G237;
magmaFloatComplex G239 = s_track[4] * G182;
magmaFloatComplex G240 = G238 + G239;
magmaFloatComplex G241 = s_track[5] * G194;
magmaFloatComplex G242 = G240 + G241;
magmaFloatComplex G243 = C3 * s_startCoefs[17];
magmaFloatComplex G244 = G243 + s_targetCoefs[17];
magmaFloatComplex G245 = G242 + G244;

r_cgesvA[0] = C0;
r_cgesvA[1] = C1;
r_cgesvA[2] = G5;
r_cgesvA[3] = G9;
r_cgesvA[4] = C1;
r_cgesvA[5] = C1;
r_cgesvA[6] = G16;
r_cgesvA[7] = C1;
r_cgesvB[0] = -G107;

r_cgesvA[8] = C1;
r_cgesvA[9] = C0;
r_cgesvA[10] = G20;
r_cgesvA[11] = G24;
r_cgesvA[12] = C1;
r_cgesvA[13] = C1;
r_cgesvA[14] = G31;
r_cgesvA[15] = C1;
r_cgesvB[1] = -G123;

r_cgesvA[16] = G5;
r_cgesvA[17] = G20;
r_cgesvA[18] = C0;
r_cgesvA[19] = C1;
r_cgesvA[20] = G35;
r_cgesvA[21] = G39;
r_cgesvA[22] = G46;
r_cgesvA[23] = G53;
r_cgesvB[2] = -G151;

r_cgesvA[24] = G9;
r_cgesvA[25] = G24;
r_cgesvA[26] = C1;
r_cgesvA[27] = C0;
r_cgesvA[28] = G57;
r_cgesvA[29] = G61;
r_cgesvA[30] = G68;
r_cgesvA[31] = G75;
r_cgesvB[3] = -G175;

r_cgesvA[32] = C1;
r_cgesvA[33] = C1;
r_cgesvA[34] = G35;
r_cgesvA[35] = G57;
r_cgesvA[36] = C0;
r_cgesvA[37] = C1;
r_cgesvA[38] = C1;
r_cgesvA[39] = G82;
r_cgesvB[4] = -G189;

r_cgesvA[40] = C1;
r_cgesvA[41] = C1;
r_cgesvA[42] = G39;
r_cgesvA[43] = G61;
r_cgesvA[44] = C1;
r_cgesvA[45] = C0;
r_cgesvA[46] = C1;
r_cgesvA[47] = G89;
r_cgesvB[5] = -G201;

r_cgesvA[48] = G16;
r_cgesvA[49] = G31;
r_cgesvA[50] = G46;
r_cgesvA[51] = G68;
r_cgesvA[52] = C1;
r_cgesvA[53] = C1;
r_cgesvA[54] = C1;
r_cgesvA[55] = C1;
r_cgesvB[6] = -G223;

r_cgesvA[56] = C1;
r_cgesvA[57] = C1;
r_cgesvA[58] = G53;
r_cgesvA[59] = G75;
r_cgesvA[60] = G82;
r_cgesvA[61] = G89;
r_cgesvA[62] = C1;
r_cgesvA[63] = C1;
r_cgesvB[7] = -G245;

  }
}

#endif
