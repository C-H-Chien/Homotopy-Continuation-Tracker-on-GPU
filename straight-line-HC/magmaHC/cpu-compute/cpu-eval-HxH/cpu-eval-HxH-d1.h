#ifndef cpu_eval_HxH_d1_h
#define cpu_eval_HxH_d1_h
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
  void cpu_eval_HxH_d1(
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
    magmaFloatComplex G91 = G4 * G90;
    magmaFloatComplex G92 = s_track[1] * s_track[1];
    magmaFloatComplex G93 = G4 * G92;
    magmaFloatComplex G94 = G91 + G93;
    magmaFloatComplex G95 = G94 + G33;
    magmaFloatComplex G96 = s_track[2] * s_track[2];
    magmaFloatComplex G97 = G4 * G96;
    magmaFloatComplex G98 = s_track[3] * s_track[3];
    magmaFloatComplex G99 = G4 * G98;
    magmaFloatComplex G100 = G97 + G99;
    magmaFloatComplex G101 = G100 + G33;
    magmaFloatComplex G102 = s_track[4] * s_track[4];
    magmaFloatComplex G103 = G4 * G102;
    magmaFloatComplex G104 = s_track[5] * s_track[5];
    magmaFloatComplex G105 = G4 * G104;
    magmaFloatComplex G106 = G103 + G105;
    magmaFloatComplex G107 = G106 + G33;
    magmaFloatComplex G108 = s_track[6] * s_track[6];
    magmaFloatComplex G109 = G4 * G108;
    magmaFloatComplex G110 = s_track[7] * s_track[7];
    magmaFloatComplex G111 = G4 * G110;
    magmaFloatComplex G112 = G109 + G111;
    magmaFloatComplex G113 = G112 + G33;
    magmaFloatComplex G114 = s_track[8] * s_track[8];
    magmaFloatComplex G115 = G4 * G114;
    magmaFloatComplex G116 = s_track[9] * s_track[9];
    magmaFloatComplex G117 = G4 * G116;
    magmaFloatComplex G118 = G115 + G117;
    magmaFloatComplex G119 = G118 + G33;
    magmaFloatComplex G120 = s_track[10] * s_track[10];
    magmaFloatComplex G121 = G4 * G120;
    magmaFloatComplex G122 = s_track[11] * s_track[11];
    magmaFloatComplex G123 = G4 * G122;
    magmaFloatComplex G124 = G121 + G123;
    magmaFloatComplex G125 = G124 + G33;
    magmaFloatComplex G126 = G9 * s_track[2];
    magmaFloatComplex G127 = G13 * s_track[4];
    magmaFloatComplex G128 = G126 + G127;
    magmaFloatComplex G129 = G128 + G59;
    magmaFloatComplex G130 = G1 * s_startCoefs[4];
    magmaFloatComplex G131 = t * s_targetCoefs[4];
    magmaFloatComplex G132 = G130 + G131;
    magmaFloatComplex G133 = G129 + G132;
    magmaFloatComplex G134 = G38 * s_track[3];
    magmaFloatComplex G135 = G47 * s_track[5];
    magmaFloatComplex G136 = G134 + G135;
    magmaFloatComplex G137 = G53 * s_track[7];
    magmaFloatComplex G138 = G136 + G137;
    magmaFloatComplex G139 = G1 * s_startCoefs[5];
    magmaFloatComplex G140 = t * s_targetCoefs[5];
    magmaFloatComplex G141 = G139 + G140;
    magmaFloatComplex G142 = G138 + G141;
    magmaFloatComplex G143 = G39 * s_track[3];
    magmaFloatComplex G144 = G48 * s_track[5];
    magmaFloatComplex G145 = G143 + G144;
    magmaFloatComplex G146 = G145 + G64;
    magmaFloatComplex G147 = G1 * s_startCoefs[6];
    magmaFloatComplex G148 = t * s_targetCoefs[6];
    magmaFloatComplex G149 = G147 + G148;
    magmaFloatComplex G150 = G146 + G149;
    magmaFloatComplex G151 = G56 * s_track[8];
    magmaFloatComplex G152 = G57 * s_track[8];
    magmaFloatComplex G153 = G151 + G152;
    magmaFloatComplex G154 = G59 * s_track[8];
    magmaFloatComplex G155 = G153 + G154;
    magmaFloatComplex G156 = G1 * s_startCoefs[7];
    magmaFloatComplex G157 = t * s_targetCoefs[7];
    magmaFloatComplex G158 = G156 + G157;
    magmaFloatComplex G159 = G155 + G158;
    magmaFloatComplex G160 = G53 * s_track[9];
    magmaFloatComplex G161 = G61 * s_track[8];
    magmaFloatComplex G162 = G160 + G161;
    magmaFloatComplex G163 = G62 * s_track[8];
    magmaFloatComplex G164 = G162 + G163;
    magmaFloatComplex G165 = G64 * s_track[8];
    magmaFloatComplex G166 = G164 + G165;
    magmaFloatComplex G167 = G1 * s_startCoefs[8];
    magmaFloatComplex G168 = t * s_targetCoefs[8];
    magmaFloatComplex G169 = G167 + G168;
    magmaFloatComplex G170 = G166 + G169;
    magmaFloatComplex G171 = G78 * s_track[10];
    magmaFloatComplex G172 = G85 * s_track[11];
    magmaFloatComplex G173 = G171 + G172;
    magmaFloatComplex G174 = G79 * s_track[10];
    magmaFloatComplex G175 = G173 + G174;
    magmaFloatComplex G176 = G86 * s_track[11];
    magmaFloatComplex G177 = G175 + G176;
    magmaFloatComplex G178 = G81 * s_track[10];
    magmaFloatComplex G179 = G177 + G178;
    magmaFloatComplex G180 = G88 * s_track[11];
    magmaFloatComplex G181 = G179 + G180;
    magmaFloatComplex G182 = G1 * s_startCoefs[9];
    magmaFloatComplex G183 = t * s_targetCoefs[9];
    magmaFloatComplex G184 = G182 + G183;
    magmaFloatComplex G185 = G181 + G184;

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
    r_cgesvB[0] = G95;

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
    r_cgesvB[1] = G101;

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
    r_cgesvB[2] = G107;

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
    r_cgesvB[3] = G113;

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
    r_cgesvB[4] = G119;

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
    r_cgesvB[5] = G125;

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
    r_cgesvB[6] = G133;

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
    r_cgesvB[7] = G142;

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
    r_cgesvB[8] = G150;

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
    r_cgesvB[9] = G159;

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
    r_cgesvB[10] = G170;

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
    r_cgesvB[11] = G185;
  }
}

#endif
