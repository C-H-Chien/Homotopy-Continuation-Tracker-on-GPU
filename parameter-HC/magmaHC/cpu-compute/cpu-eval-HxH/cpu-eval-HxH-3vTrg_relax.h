#ifndef cpu_eval_HxH_3vTrg_relax_h
#define cpu_eval_HxH_3vTrg_relax_h
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
  void cpu_eval_HxH_3vTrg_relax(
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
magmaFloatComplex G90 = C0 * s_track[0];
magmaFloatComplex G91 = s_track[2] * s_track[6];
magmaFloatComplex G92 = G91 * G4;
magmaFloatComplex G93 = G90 + G92;
magmaFloatComplex G94 = s_track[3] * s_track[6];
magmaFloatComplex G95 = G94 * G8;
magmaFloatComplex G96 = G93 + G95;
magmaFloatComplex G97 = s_track[6] * G15;
magmaFloatComplex G98 = G96 + G97;
magmaFloatComplex G99 = G1 * s_startCoefs[18];
magmaFloatComplex G100 = t * s_targetCoefs[18];
magmaFloatComplex G101 = G99 + G100;
magmaFloatComplex G102 = C0 * G101;
magmaFloatComplex G103 = C3 * G102;
magmaFloatComplex G104 = G98 + G103;
magmaFloatComplex G105 = C0 * s_track[1];
magmaFloatComplex G106 = G91 * G19;
magmaFloatComplex G107 = G105 + G106;
magmaFloatComplex G108 = G94 * G23;
magmaFloatComplex G109 = G107 + G108;
magmaFloatComplex G110 = s_track[6] * G30;
magmaFloatComplex G111 = G109 + G110;
magmaFloatComplex G112 = G1 * s_startCoefs[19];
magmaFloatComplex G113 = t * s_targetCoefs[19];
magmaFloatComplex G114 = G112 + G113;
magmaFloatComplex G115 = C0 * G114;
magmaFloatComplex G116 = C3 * G115;
magmaFloatComplex G117 = G111 + G116;
magmaFloatComplex G118 = s_track[0] * s_track[6];
magmaFloatComplex G119 = G118 * G4;
magmaFloatComplex G120 = s_track[1] * s_track[6];
magmaFloatComplex G121 = G120 * G19;
magmaFloatComplex G122 = G119 + G121;
magmaFloatComplex G123 = C0 * s_track[2];
magmaFloatComplex G124 = G122 + G123;
magmaFloatComplex G125 = s_track[4] * s_track[7];
magmaFloatComplex G126 = G125 * G34;
magmaFloatComplex G127 = G124 + G126;
magmaFloatComplex G128 = s_track[5] * s_track[7];
magmaFloatComplex G129 = G128 * G38;
magmaFloatComplex G130 = G127 + G129;
magmaFloatComplex G131 = s_track[6] * G45;
magmaFloatComplex G132 = G130 + G131;
magmaFloatComplex G133 = s_track[7] * G52;
magmaFloatComplex G134 = G132 + G133;
magmaFloatComplex G135 = G1 * s_startCoefs[20];
magmaFloatComplex G136 = t * s_targetCoefs[20];
magmaFloatComplex G137 = G135 + G136;
magmaFloatComplex G138 = C0 * G137;
magmaFloatComplex G139 = C3 * G138;
magmaFloatComplex G140 = G134 + G139;
magmaFloatComplex G141 = G118 * G8;
magmaFloatComplex G142 = G120 * G23;
magmaFloatComplex G143 = G141 + G142;
magmaFloatComplex G144 = C0 * s_track[3];
magmaFloatComplex G145 = G143 + G144;
magmaFloatComplex G146 = G125 * G56;
magmaFloatComplex G147 = G145 + G146;
magmaFloatComplex G148 = G128 * G60;
magmaFloatComplex G149 = G147 + G148;
magmaFloatComplex G150 = s_track[6] * G67;
magmaFloatComplex G151 = G149 + G150;
magmaFloatComplex G152 = s_track[7] * G74;
magmaFloatComplex G153 = G151 + G152;
magmaFloatComplex G154 = G1 * s_startCoefs[21];
magmaFloatComplex G155 = t * s_targetCoefs[21];
magmaFloatComplex G156 = G154 + G155;
magmaFloatComplex G157 = C0 * G156;
magmaFloatComplex G158 = C3 * G157;
magmaFloatComplex G159 = G153 + G158;
magmaFloatComplex G160 = s_track[2] * s_track[7];
magmaFloatComplex G161 = G160 * G34;
magmaFloatComplex G162 = s_track[3] * s_track[7];
magmaFloatComplex G163 = G162 * G56;
magmaFloatComplex G164 = G161 + G163;
magmaFloatComplex G165 = C0 * s_track[4];
magmaFloatComplex G166 = G164 + G165;
magmaFloatComplex G167 = s_track[7] * G81;
magmaFloatComplex G168 = G166 + G167;
magmaFloatComplex G169 = G1 * s_startCoefs[22];
magmaFloatComplex G170 = t * s_targetCoefs[22];
magmaFloatComplex G171 = G169 + G170;
magmaFloatComplex G172 = C0 * G171;
magmaFloatComplex G173 = C3 * G172;
magmaFloatComplex G174 = G168 + G173;
magmaFloatComplex G175 = G160 * G38;
magmaFloatComplex G176 = G162 * G60;
magmaFloatComplex G177 = G175 + G176;
magmaFloatComplex G178 = C0 * s_track[5];
magmaFloatComplex G179 = G177 + G178;
magmaFloatComplex G180 = s_track[7] * G88;
magmaFloatComplex G181 = G179 + G180;
magmaFloatComplex G182 = G1 * s_startCoefs[23];
magmaFloatComplex G183 = t * s_targetCoefs[23];
magmaFloatComplex G184 = G182 + G183;
magmaFloatComplex G185 = C0 * G184;
magmaFloatComplex G186 = C3 * G185;
magmaFloatComplex G187 = G181 + G186;
magmaFloatComplex G188 = s_track[0] * s_track[2];
magmaFloatComplex G189 = G188 * G4;
magmaFloatComplex G190 = s_track[0] * s_track[3];
magmaFloatComplex G191 = G190 * G8;
magmaFloatComplex G192 = G189 + G191;
magmaFloatComplex G193 = s_track[0] * G15;
magmaFloatComplex G194 = G192 + G193;
magmaFloatComplex G195 = s_track[1] * s_track[2];
magmaFloatComplex G196 = G195 * G19;
magmaFloatComplex G197 = G194 + G196;
magmaFloatComplex G198 = s_track[1] * s_track[3];
magmaFloatComplex G199 = G198 * G23;
magmaFloatComplex G200 = G197 + G199;
magmaFloatComplex G201 = s_track[1] * G30;
magmaFloatComplex G202 = G200 + G201;
magmaFloatComplex G203 = s_track[2] * G45;
magmaFloatComplex G204 = G202 + G203;
magmaFloatComplex G205 = s_track[3] * G67;
magmaFloatComplex G206 = G204 + G205;
magmaFloatComplex G207 = G1 * s_startCoefs[8];
magmaFloatComplex G208 = t * s_targetCoefs[8];
magmaFloatComplex G209 = G207 + G208;
magmaFloatComplex G210 = G206 + G209;
magmaFloatComplex G211 = s_track[2] * s_track[4];
magmaFloatComplex G212 = G211 * G34;
magmaFloatComplex G213 = s_track[2] * s_track[5];
magmaFloatComplex G214 = G213 * G38;
magmaFloatComplex G215 = G212 + G214;
magmaFloatComplex G216 = s_track[2] * G52;
magmaFloatComplex G217 = G215 + G216;
magmaFloatComplex G218 = s_track[3] * s_track[4];
magmaFloatComplex G219 = G218 * G56;
magmaFloatComplex G220 = G217 + G219;
magmaFloatComplex G221 = s_track[3] * s_track[5];
magmaFloatComplex G222 = G221 * G60;
magmaFloatComplex G223 = G220 + G222;
magmaFloatComplex G224 = s_track[3] * G74;
magmaFloatComplex G225 = G223 + G224;
magmaFloatComplex G226 = s_track[4] * G81;
magmaFloatComplex G227 = G225 + G226;
magmaFloatComplex G228 = s_track[5] * G88;
magmaFloatComplex G229 = G227 + G228;
magmaFloatComplex G230 = G1 * s_startCoefs[17];
magmaFloatComplex G231 = t * s_targetCoefs[17];
magmaFloatComplex G232 = G230 + G231;
magmaFloatComplex G233 = G229 + G232;

r_cgesvA[0] = C0;
r_cgesvA[1] = C1;
r_cgesvA[2] = G5;
r_cgesvA[3] = G9;
r_cgesvA[4] = C1;
r_cgesvA[5] = C1;
r_cgesvA[6] = G16;
r_cgesvA[7] = C1;
r_cgesvB[0] =G104;

r_cgesvA[8] = C1;
r_cgesvA[9] = C0;
r_cgesvA[10] = G20;
r_cgesvA[11] = G24;
r_cgesvA[12] = C1;
r_cgesvA[13] = C1;
r_cgesvA[14] = G31;
r_cgesvA[15] = C1;
r_cgesvB[1] =G117;

r_cgesvA[16] = G5;
r_cgesvA[17] = G20;
r_cgesvA[18] = C0;
r_cgesvA[19] = C1;
r_cgesvA[20] = G35;
r_cgesvA[21] = G39;
r_cgesvA[22] = G46;
r_cgesvA[23] = G53;
r_cgesvB[2] =G140;

r_cgesvA[24] = G9;
r_cgesvA[25] = G24;
r_cgesvA[26] = C1;
r_cgesvA[27] = C0;
r_cgesvA[28] = G57;
r_cgesvA[29] = G61;
r_cgesvA[30] = G68;
r_cgesvA[31] = G75;
r_cgesvB[3] =G159;

r_cgesvA[32] = C1;
r_cgesvA[33] = C1;
r_cgesvA[34] = G35;
r_cgesvA[35] = G57;
r_cgesvA[36] = C0;
r_cgesvA[37] = C1;
r_cgesvA[38] = C1;
r_cgesvA[39] = G82;
r_cgesvB[4] =G174;

r_cgesvA[40] = C1;
r_cgesvA[41] = C1;
r_cgesvA[42] = G39;
r_cgesvA[43] = G61;
r_cgesvA[44] = C1;
r_cgesvA[45] = C0;
r_cgesvA[46] = C1;
r_cgesvA[47] = G89;
r_cgesvB[5] =G187;

r_cgesvA[48] = G16;
r_cgesvA[49] = G31;
r_cgesvA[50] = G46;
r_cgesvA[51] = G68;
r_cgesvA[52] = C1;
r_cgesvA[53] = C1;
r_cgesvA[54] = C1;
r_cgesvA[55] = C1;
r_cgesvB[6] =G210;

r_cgesvA[56] = C1;
r_cgesvA[57] = C1;
r_cgesvA[58] = G53;
r_cgesvA[59] = G75;
r_cgesvA[60] = G82;
r_cgesvA[61] = G89;
r_cgesvA[62] = C1;
r_cgesvA[63] = C1;
r_cgesvB[7] =G233;

  }
}

#endif
