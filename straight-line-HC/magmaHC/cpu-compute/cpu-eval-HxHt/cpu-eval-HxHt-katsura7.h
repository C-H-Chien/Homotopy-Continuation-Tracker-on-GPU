#ifndef cpu_eval_HxHt_katsura7_h
#define cpu_eval_HxHt_katsura7_h
// ============================================================================
// partial derivative evaluations of the katsura7 problem for cpu HC computation
//
// Modifications
//    Chien  21-12-17:   Originally created
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
  void cpu_eval_HxHt_katsura7(
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
    magmaFloatComplex G20 = s_track[1] + s_track[1];
    magmaFloatComplex G21 = G13 * G20;
    magmaFloatComplex G22 = G13 * s_track[0];
    magmaFloatComplex G23 = G22 + G15;
    magmaFloatComplex G24 = G23 + G9;
    magmaFloatComplex G25 = G4 * G20;
    magmaFloatComplex G26 = G25 + G16;
    magmaFloatComplex G27 = G15 + G17;
    magmaFloatComplex G28 = G16 + G18;
    magmaFloatComplex G29 = G17 + G19;
    magmaFloatComplex G30 = s_track[7] * G13;
    magmaFloatComplex G31 = G18 + G30;
    magmaFloatComplex G32 = s_track[2] + s_track[2];
    magmaFloatComplex G33 = G13 * G32;
    magmaFloatComplex G34 = G13 * s_track[1];
    magmaFloatComplex G35 = G34 + G16;
    magmaFloatComplex G36 = G22 + G17;
    magmaFloatComplex G37 = G36 + G9;
    magmaFloatComplex G38 = G34 + G18;
    magmaFloatComplex G39 = G4 * G32;
    magmaFloatComplex G40 = G39 + G19;
    magmaFloatComplex G41 = G16 + G30;
    magmaFloatComplex G42 = s_track[3] + s_track[3];
    magmaFloatComplex G43 = G13 * G42;
    magmaFloatComplex G44 = G13 * s_track[2];
    magmaFloatComplex G45 = G44 + G17;
    magmaFloatComplex G46 = G22 + G19;
    magmaFloatComplex G47 = G46 + G9;
    magmaFloatComplex G48 = G34 + G30;
    magmaFloatComplex G49 = G4 * G42;
    magmaFloatComplex G50 = s_track[4] + s_track[4];
    magmaFloatComplex G51 = G13 * G50;
    magmaFloatComplex G52 = G13 * s_track[3];
    magmaFloatComplex G53 = G52 + G18;
    magmaFloatComplex G54 = G44 + G19;
    magmaFloatComplex G55 = G22 + G9;
    magmaFloatComplex G56 = s_track[5] + s_track[5];
    magmaFloatComplex G57 = G13 * G56;
    magmaFloatComplex G58 = G13 * s_track[4];
    magmaFloatComplex G59 = G58 + G19;
    magmaFloatComplex G60 = G52 + G30;
    magmaFloatComplex G61 = s_track[6] + s_track[6];
    magmaFloatComplex G62 = G13 * G61;
    magmaFloatComplex G63 = G13 * s_track[5];
    magmaFloatComplex G64 = G63 + G30;
    magmaFloatComplex G65 = s_track[7] + s_track[7];
    magmaFloatComplex G66 = G13 * G65;
    magmaFloatComplex G67 = G13 * s_track[6];
    magmaFloatComplex G68 = s_track[0] * s_track[0];
    magmaFloatComplex G70 = s_targetCoefs[0] - s_startCoefs[0];
    magmaFloatComplex G71 = G68 * G70;
    magmaFloatComplex G73 = s_targetCoefs[1] - s_startCoefs[1];
    magmaFloatComplex G74 = s_track[0] * G73;
    magmaFloatComplex G75 = G71 + G74;
    magmaFloatComplex G76 = s_track[1] * s_track[1];
    magmaFloatComplex G78 = s_targetCoefs[2] - s_startCoefs[2];
    magmaFloatComplex G79 = G76 * G78;
    magmaFloatComplex G80 = G75 + G79;
    magmaFloatComplex G81 = s_track[2] * s_track[2];
    magmaFloatComplex G82 = G81 * G78;
    magmaFloatComplex G83 = G80 + G82;
    magmaFloatComplex G84 = s_track[3] * s_track[3];
    magmaFloatComplex G85 = G84 * G78;
    magmaFloatComplex G86 = G83 + G85;
    magmaFloatComplex G87 = s_track[4] * s_track[4];
    magmaFloatComplex G88 = G87 * G78;
    magmaFloatComplex G89 = G86 + G88;
    magmaFloatComplex G90 = s_track[5] * s_track[5];
    magmaFloatComplex G91 = G90 * G78;
    magmaFloatComplex G92 = G89 + G91;
    magmaFloatComplex G93 = s_track[6] * s_track[6];
    magmaFloatComplex G94 = G93 * G78;
    magmaFloatComplex G95 = G92 + G94;
    magmaFloatComplex G96 = s_track[7] * s_track[7];
    magmaFloatComplex G97 = G96 * G78;
    magmaFloatComplex G98 = G95 + G97;
    magmaFloatComplex G99 = s_track[0] * G78;
    magmaFloatComplex G100 = s_track[1] * G99;
    magmaFloatComplex G101 = s_track[1] * G78;
    magmaFloatComplex G102 = s_track[2] * G101;
    magmaFloatComplex G103 = G100 + G102;
    magmaFloatComplex G104 = s_track[1] * G73;
    magmaFloatComplex G105 = G103 + G104;
    magmaFloatComplex G106 = s_track[2] * G78;
    magmaFloatComplex G107 = s_track[3] * G106;
    magmaFloatComplex G108 = G105 + G107;
    magmaFloatComplex G109 = s_track[3] * G78;
    magmaFloatComplex G110 = s_track[4] * G109;
    magmaFloatComplex G111 = G108 + G110;
    magmaFloatComplex G112 = s_track[4] * G78;
    magmaFloatComplex G113 = s_track[5] * G112;
    magmaFloatComplex G114 = G111 + G113;
    magmaFloatComplex G115 = s_track[5] * G78;
    magmaFloatComplex G116 = s_track[6] * G115;
    magmaFloatComplex G117 = G114 + G116;
    magmaFloatComplex G118 = s_track[6] * G78;
    magmaFloatComplex G119 = s_track[7] * G118;
    magmaFloatComplex G120 = G117 + G119;
    magmaFloatComplex G121 = s_track[2] * G99;
    magmaFloatComplex G122 = G76 * G70;
    magmaFloatComplex G123 = G121 + G122;
    magmaFloatComplex G124 = s_track[3] * G101;
    magmaFloatComplex G125 = G123 + G124;
    magmaFloatComplex G126 = s_track[4] * G106;
    magmaFloatComplex G127 = G125 + G126;
    magmaFloatComplex G128 = s_track[2] * G73;
    magmaFloatComplex G129 = G127 + G128;
    magmaFloatComplex G130 = s_track[5] * G109;
    magmaFloatComplex G131 = G129 + G130;
    magmaFloatComplex G132 = s_track[6] * G112;
    magmaFloatComplex G133 = G131 + G132;
    magmaFloatComplex G134 = s_track[7] * G115;
    magmaFloatComplex G135 = G133 + G134;
    magmaFloatComplex G136 = s_track[3] * G99;
    magmaFloatComplex G137 = G136 + G102;
    magmaFloatComplex G138 = s_track[4] * G101;
    magmaFloatComplex G139 = G137 + G138;
    magmaFloatComplex G140 = s_track[5] * G106;
    magmaFloatComplex G141 = G139 + G140;
    magmaFloatComplex G142 = s_track[6] * G109;
    magmaFloatComplex G143 = G141 + G142;
    magmaFloatComplex G144 = s_track[3] * G73;
    magmaFloatComplex G145 = G143 + G144;
    magmaFloatComplex G146 = s_track[7] * G112;
    magmaFloatComplex G147 = G145 + G146;
    magmaFloatComplex G148 = s_track[4] * G99;
    magmaFloatComplex G149 = G148 + G124;
    magmaFloatComplex G150 = s_track[5] * G101;
    magmaFloatComplex G151 = G149 + G150;
    magmaFloatComplex G152 = G81 * G70;
    magmaFloatComplex G153 = G151 + G152;
    magmaFloatComplex G154 = s_track[6] * G106;
    magmaFloatComplex G155 = G153 + G154;
    magmaFloatComplex G156 = s_track[7] * G109;
    magmaFloatComplex G157 = G155 + G156;
    magmaFloatComplex G158 = s_track[4] * G73;
    magmaFloatComplex G159 = G157 + G158;
    magmaFloatComplex G160 = s_track[5] * G99;
    magmaFloatComplex G161 = G160 + G138;
    magmaFloatComplex G162 = s_track[6] * G101;
    magmaFloatComplex G163 = G161 + G162;
    magmaFloatComplex G164 = G163 + G107;
    magmaFloatComplex G165 = s_track[7] * G106;
    magmaFloatComplex G166 = G164 + G165;
    magmaFloatComplex G167 = s_track[5] * G73;
    magmaFloatComplex G168 = G166 + G167;
    magmaFloatComplex G169 = s_track[6] * G99;
    magmaFloatComplex G170 = G169 + G150;
    magmaFloatComplex G171 = s_track[7] * G101;
    magmaFloatComplex G172 = G170 + G171;
    magmaFloatComplex G173 = G172 + G126;
    magmaFloatComplex G174 = G84 * G70;
    magmaFloatComplex G175 = G173 + G174;
    magmaFloatComplex G176 = s_track[6] * G73;
    magmaFloatComplex G177 = G175 + G176;
    magmaFloatComplex G178 = s_track[0] * G70;
    magmaFloatComplex G179 = G178 + G101;
    magmaFloatComplex G180 = G179 + G106;
    magmaFloatComplex G181 = G180 + G109;
    magmaFloatComplex G182 = G181 + G112;
    magmaFloatComplex G183 = G182 + G115;
    magmaFloatComplex G184 = G183 + G118;
    magmaFloatComplex G185 = s_track[7] * G78;
    magmaFloatComplex G186 = G184 + G185;
    magmaFloatComplex G187 = G186 + G73;

    r_cgesvA[0] = G10;
    r_cgesvA[1] = G14;
    r_cgesvA[2] = G15;
    r_cgesvA[3] = G16;
    r_cgesvA[4] = G17;
    r_cgesvA[5] = G18;
    r_cgesvA[6] = G19;
    r_cgesvA[7] = G4;
    r_cgesvB[0] = -G98;

    r_cgesvA[8] = G21;
    r_cgesvA[9] = G24;
    r_cgesvA[10] = G26;
    r_cgesvA[11] = G27;
    r_cgesvA[12] = G28;
    r_cgesvA[13] = G29;
    r_cgesvA[14] = G31;
    r_cgesvA[15] = G13;
    r_cgesvB[1] = -G120;

    r_cgesvA[16] = G33;
    r_cgesvA[17] = G35;
    r_cgesvA[18] = G37;
    r_cgesvA[19] = G38;
    r_cgesvA[20] = G40;
    r_cgesvA[21] = G41;
    r_cgesvA[22] = G17;
    r_cgesvA[23] = G13;
    r_cgesvB[2] = -G135;

    r_cgesvA[24] = G43;
    r_cgesvA[25] = G45;
    r_cgesvA[26] = G38;
    r_cgesvA[27] = G47;
    r_cgesvA[28] = G48;
    r_cgesvA[29] = G44;
    r_cgesvA[30] = G49;
    r_cgesvA[31] = G13;
    r_cgesvB[3] = -G147;

    r_cgesvA[32] = G51;
    r_cgesvA[33] = G53;
    r_cgesvA[34] = G54;
    r_cgesvA[35] = G48;
    r_cgesvA[36] = G55;
    r_cgesvA[37] = G34;
    r_cgesvA[38] = G44;
    r_cgesvA[39] = G13;
    r_cgesvB[4] = -G159;

    r_cgesvA[40] = G57;
    r_cgesvA[41] = G59;
    r_cgesvA[42] = G60;
    r_cgesvA[43] = G44;
    r_cgesvA[44] = G34;
    r_cgesvA[45] = G55;
    r_cgesvA[46] = G34;
    r_cgesvA[47] = G13;
    r_cgesvB[5] = -G168;

    r_cgesvA[48] = G62;
    r_cgesvA[49] = G64;
    r_cgesvA[50] = G58;
    r_cgesvA[51] = G52;
    r_cgesvA[52] = G44;
    r_cgesvA[53] = G34;
    r_cgesvA[54] = G55;
    r_cgesvA[55] = G13;
    r_cgesvB[6] = -G177;

    r_cgesvA[56] = G66;
    r_cgesvA[57] = G67;
    r_cgesvA[58] = G63;
    r_cgesvA[59] = G58;
    r_cgesvA[60] = G52;
    r_cgesvA[61] = G44;
    r_cgesvA[62] = G34;
    r_cgesvA[63] = G13;
    r_cgesvB[7] = -G187;

  }
}

#endif
