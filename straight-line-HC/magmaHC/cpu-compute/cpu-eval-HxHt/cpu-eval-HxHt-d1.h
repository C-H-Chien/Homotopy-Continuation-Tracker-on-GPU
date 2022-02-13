#ifndef cpu_eval_HxHt_d1_h
#define cpu_eval_HxHt_d1_h
// ============================================================================
// partial derivative evaluations of the d1 problem for cpu HC computation
//
// Modifications
//    Chien  21-09-14:   Originally created
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
  void cpu_eval_HxHt_d1(
      magma_int_t s, float t, int N, magmaFloatComplex* s_track,
      const magmaFloatComplex &C0, const magmaFloatComplex &C1, const magmaFloatComplex &C2,
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
    magmaFloatComplex G7 = G1 * s_startCoefs[2];
    magmaFloatComplex G8 = t * s_targetCoefs[2];
    magmaFloatComplex G9 = G7 + G8;
    magmaFloatComplex G10 = s_track[3] * G9;
    magmaFloatComplex G11 = G1 * s_startCoefs[3];
    magmaFloatComplex G12 = t * s_targetCoefs[3];
    magmaFloatComplex G13 = G11 + G12;
    magmaFloatComplex G14 = s_track[5] * G13;
    magmaFloatComplex G15 = G10 + G14;
    magmaFloatComplex G16 = s_track[7] * G4;
    magmaFloatComplex G17 = G15 + G16;
    magmaFloatComplex G18 = s_track[9] * G4;
    magmaFloatComplex G19 = s_track[1] + s_track[1];
    magmaFloatComplex G20 = G4 * G19;
    magmaFloatComplex G21 = s_track[3] * G4;
    magmaFloatComplex G22 = s_track[8] * G21;
    magmaFloatComplex G23 = s_track[5] * G4;
    magmaFloatComplex G24 = s_track[8] * G23;
    magmaFloatComplex G25 = G22 + G24;
    magmaFloatComplex G26 = s_track[8] * G16;
    magmaFloatComplex G27 = G25 + G26;
    magmaFloatComplex G28 = s_track[2] + s_track[2];
    magmaFloatComplex G29 = G4 * G28;
    magmaFloatComplex G30 = s_track[8] * G4;
    magmaFloatComplex G31 = G1 * s_startCoefs[1];
    magmaFloatComplex G32 = t * s_targetCoefs[1];
    magmaFloatComplex G33 = G31 + G32;
    magmaFloatComplex G34 = s_track[9] * G33;
    magmaFloatComplex G35 = s_track[10] * G34;
    magmaFloatComplex G36 = s_track[3] + s_track[3];
    magmaFloatComplex G37 = G4 * G36;
    magmaFloatComplex G38 = G9 * s_track[0];
    magmaFloatComplex G39 = G9 * s_track[1];
    magmaFloatComplex G40 = G4 * s_track[1];
    magmaFloatComplex G41 = s_track[8] * G40;
    magmaFloatComplex G42 = s_track[11] * G4;
    magmaFloatComplex G43 = s_track[4] + s_track[4];
    magmaFloatComplex G44 = G4 * G43;
    magmaFloatComplex G45 = s_track[5] + s_track[5];
    magmaFloatComplex G46 = G4 * G45;
    magmaFloatComplex G47 = G13 * s_track[0];
    magmaFloatComplex G48 = G13 * s_track[1];
    magmaFloatComplex G49 = s_track[6] + s_track[6];
    magmaFloatComplex G50 = G4 * G49;
    magmaFloatComplex G51 = s_track[7] + s_track[7];
    magmaFloatComplex G52 = G4 * G51;
    magmaFloatComplex G53 = G4 * s_track[0];
    magmaFloatComplex G54 = s_track[8] + s_track[8];
    magmaFloatComplex G55 = G4 * G54;
    magmaFloatComplex G56 = G4 * s_track[2];
    magmaFloatComplex G57 = G4 * s_track[4];
    magmaFloatComplex G58 = G56 + G57;
    magmaFloatComplex G59 = G4 * s_track[6];
    magmaFloatComplex G60 = G58 + G59;
    magmaFloatComplex G61 = G40 * s_track[3];
    magmaFloatComplex G62 = G40 * s_track[5];
    magmaFloatComplex G63 = G61 + G62;
    magmaFloatComplex G64 = G40 * s_track[7];
    magmaFloatComplex G65 = G63 + G64;
    magmaFloatComplex G66 = s_track[9] + s_track[9];
    magmaFloatComplex G67 = G4 * G66;
    magmaFloatComplex G68 = G33 * s_track[2];
    magmaFloatComplex G69 = s_track[10] * G68;
    magmaFloatComplex G70 = G33 * s_track[4];
    magmaFloatComplex G71 = s_track[10] * G70;
    magmaFloatComplex G72 = G69 + G71;
    magmaFloatComplex G73 = G33 * s_track[6];
    magmaFloatComplex G74 = s_track[10] * G73;
    magmaFloatComplex G75 = G72 + G74;
    magmaFloatComplex G76 = s_track[10] + s_track[10];
    magmaFloatComplex G77 = G4 * G76;
    magmaFloatComplex G78 = G68 * s_track[9];
    magmaFloatComplex G79 = G70 * s_track[9];
    magmaFloatComplex G80 = G78 + G79;
    magmaFloatComplex G81 = G73 * s_track[9];
    magmaFloatComplex G82 = G80 + G81;
    magmaFloatComplex G83 = s_track[11] + s_track[11];
    magmaFloatComplex G84 = G4 * G83;
    magmaFloatComplex G85 = G4 * s_track[3];
    magmaFloatComplex G86 = G4 * s_track[5];
    magmaFloatComplex G87 = G85 + G86;
    magmaFloatComplex G88 = G4 * s_track[7];
    magmaFloatComplex G89 = G87 + G88;
    magmaFloatComplex G90 = s_track[0] * s_track[0];
    magmaFloatComplex G91 = C1 * s_startCoefs[0];
    magmaFloatComplex G92 = G91 + s_targetCoefs[0];
    magmaFloatComplex G93 = G90 * G92;
    magmaFloatComplex G94 = s_track[1] * s_track[1];
    magmaFloatComplex G95 = G94 * G92;
    magmaFloatComplex G96 = G93 + G95;
    magmaFloatComplex G97 = C1 * s_startCoefs[1];
    magmaFloatComplex G98 = G97 + s_targetCoefs[1];
    magmaFloatComplex G99 = G96 + G98;
    magmaFloatComplex G100 = s_track[2] * s_track[2];
    magmaFloatComplex G101 = G100 * G92;
    magmaFloatComplex G102 = s_track[3] * s_track[3];
    magmaFloatComplex G103 = G102 * G92;
    magmaFloatComplex G104 = G101 + G103;
    magmaFloatComplex G105 = G104 + G98;
    magmaFloatComplex G106 = s_track[4] * s_track[4];
    magmaFloatComplex G107 = G106 * G92;
    magmaFloatComplex G108 = s_track[5] * s_track[5];
    magmaFloatComplex G109 = G108 * G92;
    magmaFloatComplex G110 = G107 + G109;
    magmaFloatComplex G111 = G110 + G98;
    magmaFloatComplex G112 = s_track[6] * s_track[6];
    magmaFloatComplex G113 = G112 * G92;
    magmaFloatComplex G114 = s_track[7] * s_track[7];
    magmaFloatComplex G115 = G114 * G92;
    magmaFloatComplex G116 = G113 + G115;
    magmaFloatComplex G117 = G116 + G98;
    magmaFloatComplex G118 = s_track[8] * s_track[8];
    magmaFloatComplex G119 = G118 * G92;
    magmaFloatComplex G120 = s_track[9] * s_track[9];
    magmaFloatComplex G121 = G120 * G92;
    magmaFloatComplex G122 = G119 + G121;
    magmaFloatComplex G123 = G122 + G98;
    magmaFloatComplex G124 = s_track[10] * s_track[10];
    magmaFloatComplex G125 = G124 * G92;
    magmaFloatComplex G126 = s_track[11] * s_track[11];
    magmaFloatComplex G127 = G126 * G92;
    magmaFloatComplex G128 = G125 + G127;
    magmaFloatComplex G129 = G128 + G98;
    magmaFloatComplex G130 = C1 * s_startCoefs[2];
    magmaFloatComplex G131 = G130 + s_targetCoefs[2];
    magmaFloatComplex G132 = s_track[2] * G131;
    magmaFloatComplex G133 = C1 * s_startCoefs[3];
    magmaFloatComplex G134 = G133 + s_targetCoefs[3];
    magmaFloatComplex G135 = s_track[4] * G134;
    magmaFloatComplex G136 = G132 + G135;
    magmaFloatComplex G137 = s_track[6] * G92;
    magmaFloatComplex G138 = G136 + G137;
    magmaFloatComplex G139 = C1 * s_startCoefs[4];
    magmaFloatComplex G140 = G139 + s_targetCoefs[4];
    magmaFloatComplex G141 = G138 + G140;
    magmaFloatComplex G142 = s_track[0] * G131;
    magmaFloatComplex G143 = s_track[3] * G142;
    magmaFloatComplex G144 = s_track[0] * G134;
    magmaFloatComplex G145 = s_track[5] * G144;
    magmaFloatComplex G146 = G143 + G145;
    magmaFloatComplex G147 = s_track[0] * G92;
    magmaFloatComplex G148 = s_track[7] * G147;
    magmaFloatComplex G149 = G146 + G148;
    magmaFloatComplex G150 = C1 * s_startCoefs[5];
    magmaFloatComplex G151 = G150 + s_targetCoefs[5];
    magmaFloatComplex G152 = G149 + G151;
    magmaFloatComplex G153 = s_track[1] * G131;
    magmaFloatComplex G154 = s_track[3] * G153;
    magmaFloatComplex G155 = s_track[1] * G134;
    magmaFloatComplex G156 = s_track[5] * G155;
    magmaFloatComplex G157 = G154 + G156;
    magmaFloatComplex G158 = s_track[1] * G92;
    magmaFloatComplex G159 = s_track[7] * G158;
    magmaFloatComplex G160 = G157 + G159;
    magmaFloatComplex G161 = C1 * s_startCoefs[6];
    magmaFloatComplex G162 = G161 + s_targetCoefs[6];
    magmaFloatComplex G163 = G160 + G162;
    magmaFloatComplex G164 = s_track[2] * G92;
    magmaFloatComplex G165 = s_track[8] * G164;
    magmaFloatComplex G166 = s_track[4] * G92;
    magmaFloatComplex G167 = s_track[8] * G166;
    magmaFloatComplex G168 = G165 + G167;
    magmaFloatComplex G169 = s_track[8] * G137;
    magmaFloatComplex G170 = G168 + G169;
    magmaFloatComplex G171 = C1 * s_startCoefs[7];
    magmaFloatComplex G172 = G171 + s_targetCoefs[7];
    magmaFloatComplex G173 = G170 + G172;
    magmaFloatComplex G174 = s_track[9] * G147;
    magmaFloatComplex G175 = s_track[3] * G158;
    magmaFloatComplex G176 = s_track[8] * G175;
    magmaFloatComplex G177 = G174 + G176;
    magmaFloatComplex G178 = s_track[5] * G158;
    magmaFloatComplex G179 = s_track[8] * G178;
    magmaFloatComplex G180 = G177 + G179;
    magmaFloatComplex G181 = s_track[8] * G159;
    magmaFloatComplex G182 = G180 + G181;
    magmaFloatComplex G183 = C1 * s_startCoefs[8];
    magmaFloatComplex G184 = G183 + s_targetCoefs[8];
    magmaFloatComplex G185 = G182 + G184;
    magmaFloatComplex G186 = s_track[2] * G98;
    magmaFloatComplex G187 = s_track[9] * G186;
    magmaFloatComplex G188 = s_track[10] * G187;
    magmaFloatComplex G189 = s_track[3] * G92;
    magmaFloatComplex G190 = s_track[11] * G189;
    magmaFloatComplex G191 = G188 + G190;
    magmaFloatComplex G192 = s_track[4] * G98;
    magmaFloatComplex G193 = s_track[9] * G192;
    magmaFloatComplex G194 = s_track[10] * G193;
    magmaFloatComplex G195 = G191 + G194;
    magmaFloatComplex G196 = s_track[5] * G92;
    magmaFloatComplex G197 = s_track[11] * G196;
    magmaFloatComplex G198 = G195 + G197;
    magmaFloatComplex G199 = s_track[6] * G98;
    magmaFloatComplex G200 = s_track[9] * G199;
    magmaFloatComplex G201 = s_track[10] * G200;
    magmaFloatComplex G202 = G198 + G201;
    magmaFloatComplex G203 = s_track[7] * G92;
    magmaFloatComplex G204 = s_track[11] * G203;
    magmaFloatComplex G205 = G202 + G204;
    magmaFloatComplex G206 = C1 * s_startCoefs[9];
    magmaFloatComplex G207 = G206 + s_targetCoefs[9];
    magmaFloatComplex G208 = G205 + G207;

    r_cgesvA[0] = G6;
    r_cgesvA[1] = C2;
    r_cgesvA[2] = C2;
    r_cgesvA[3] = C2;
    r_cgesvA[4] = C2;
    r_cgesvA[5] = C2;
    r_cgesvA[6] = C2;
    r_cgesvA[7] = G17;
    r_cgesvA[8] = C2;
    r_cgesvA[9] = C2;
    r_cgesvA[10] = G18;
    r_cgesvA[11] = C2;
    r_cgesvB[0] = -G99;

    r_cgesvA[12] = G20;
    r_cgesvA[13] = C2;
    r_cgesvA[14] = C2;
    r_cgesvA[15] = C2;
    r_cgesvA[16] = C2;
    r_cgesvA[17] = C2;
    r_cgesvA[18] = C2;
    r_cgesvA[19] = C2;
    r_cgesvA[20] = G17;
    r_cgesvA[21] = C2;
    r_cgesvA[22] = G27;
    r_cgesvA[23] = C2;
    r_cgesvB[1] = -G105;

    r_cgesvA[24] = C2;
    r_cgesvA[25] = G29;
    r_cgesvA[26] = C2;
    r_cgesvA[27] = C2;
    r_cgesvA[28] = C2;
    r_cgesvA[29] = C2;
    r_cgesvA[30] = G9;
    r_cgesvA[31] = C2;
    r_cgesvA[32] = C2;
    r_cgesvA[33] = G30;
    r_cgesvA[34] = C2;
    r_cgesvA[35] = G35;
    r_cgesvB[2] = -G111;

    r_cgesvA[36] = C2;
    r_cgesvA[37] = G37;
    r_cgesvA[38] = C2;
    r_cgesvA[39] = C2;
    r_cgesvA[40] = C2;
    r_cgesvA[41] = C2;
    r_cgesvA[42] = C2;
    r_cgesvA[43] = G38;
    r_cgesvA[44] = G39;
    r_cgesvA[45] = C2;
    r_cgesvA[46] = G41;
    r_cgesvA[47] = G42;
    r_cgesvB[3] = -G117;

    r_cgesvA[48] = C2;
    r_cgesvA[49] = C2;
    r_cgesvA[50] = G44;
    r_cgesvA[51] = C2;
    r_cgesvA[52] = C2;
    r_cgesvA[53] = C2;
    r_cgesvA[54] = G13;
    r_cgesvA[55] = C2;
    r_cgesvA[56] = C2;
    r_cgesvA[57] = G30;
    r_cgesvA[58] = C2;
    r_cgesvA[59] = G35;
    r_cgesvB[4] = -G123;

    r_cgesvA[60] = C2;
    r_cgesvA[61] = C2;
    r_cgesvA[62] = G46;
    r_cgesvA[63] = C2;
    r_cgesvA[64] = C2;
    r_cgesvA[65] = C2;
    r_cgesvA[66] = C2;
    r_cgesvA[67] = G47;
    r_cgesvA[68] = G48;
    r_cgesvA[69] = C2;
    r_cgesvA[70] = G41;
    r_cgesvA[71] = G42;
    r_cgesvB[5] = -G129;

    r_cgesvA[72] = C2;
    r_cgesvA[73] = C2;
    r_cgesvA[74] = C2;
    r_cgesvA[75] = G50;
    r_cgesvA[76] = C2;
    r_cgesvA[77] = C2;
    r_cgesvA[78] = G4;
    r_cgesvA[79] = C2;
    r_cgesvA[80] = C2;
    r_cgesvA[81] = G30;
    r_cgesvA[82] = C2;
    r_cgesvA[83] = G35;
    r_cgesvB[6] = -G141;

    r_cgesvA[84] = C2;
    r_cgesvA[85] = C2;
    r_cgesvA[86] = C2;
    r_cgesvA[87] = G52;
    r_cgesvA[88] = C2;
    r_cgesvA[89] = C2;
    r_cgesvA[90] = C2;
    r_cgesvA[91] = G53;
    r_cgesvA[92] = G40;
    r_cgesvA[93] = C2;
    r_cgesvA[94] = G41;
    r_cgesvA[95] = G42;
    r_cgesvB[7] = -G152;

    r_cgesvA[96] = C2;
    r_cgesvA[97] = C2;
    r_cgesvA[98] = C2;
    r_cgesvA[99] = C2;
    r_cgesvA[100] = G55;
    r_cgesvA[101] = C2;
    r_cgesvA[102] = C2;
    r_cgesvA[103] = C2;
    r_cgesvA[104] = C2;
    r_cgesvA[105] = G60;
    r_cgesvA[106] = G65;
    r_cgesvA[107] = C2;
    r_cgesvB[8] = -G163;

    r_cgesvA[108] = C2;
    r_cgesvA[109] = C2;
    r_cgesvA[110] = C2;
    r_cgesvA[111] = C2;
    r_cgesvA[112] = G67;
    r_cgesvA[113] = C2;
    r_cgesvA[114] = C2;
    r_cgesvA[115] = C2;
    r_cgesvA[116] = C2;
    r_cgesvA[117] = C2;
    r_cgesvA[118] = G53;
    r_cgesvA[119] = G75;
    r_cgesvB[9] = -G173;

    r_cgesvA[120] = C2;
    r_cgesvA[121] = C2;
    r_cgesvA[122] = C2;
    r_cgesvA[123] = C2;
    r_cgesvA[124] = C2;
    r_cgesvA[125] = G77;
    r_cgesvA[126] = C2;
    r_cgesvA[127] = C2;
    r_cgesvA[128] = C2;
    r_cgesvA[129] = C2;
    r_cgesvA[130] = C2;
    r_cgesvA[131] = G82;
    r_cgesvB[10] = -G185;

    r_cgesvA[132] = C2;
    r_cgesvA[133] = C2;
    r_cgesvA[134] = C2;
    r_cgesvA[135] = C2;
    r_cgesvA[136] = C2;
    r_cgesvA[137] = G84;
    r_cgesvA[138] = C2;
    r_cgesvA[139] = C2;
    r_cgesvA[140] = C2;
    r_cgesvA[141] = C2;
    r_cgesvA[142] = C2;
    r_cgesvA[143] = G89;
    r_cgesvB[11] = -G208;
  }
}

#endif
