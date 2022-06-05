#ifndef cpu_eval_HxHt_katsura6_h
#define cpu_eval_HxHt_katsura6_h
// ============================================================================
// partial derivative evaluations of the katsura6 problem for cpu HC computation
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
  void cpu_eval_HxHt_katsura6(
      magma_int_t s, float t, int N, magmaFloatComplex* s_track,
      const magmaFloatComplex &C0, const magmaFloatComplex &C1, const magmaFloatComplex &C2,
      magmaFloatComplex* r_startCoefs, magmaFloatComplex* r_targetCoefs,
      magmaFloatComplex* r_cgesvA, magmaFloatComplex* r_cgesvB )
  {
    magmaFloatComplex G1 = C0 - t;
    magmaFloatComplex G2 = G1 * r_startCoefs[0];
    magmaFloatComplex G3 = t * r_targetCoefs[0];
    magmaFloatComplex G4 = G2 + G3;
    magmaFloatComplex G5 = G1 * r_startCoefs[1];
    magmaFloatComplex G6 = t * r_targetCoefs[1];
    magmaFloatComplex G7 = G5 + G6;
    magmaFloatComplex G8 = s_track[5] * G7;
    magmaFloatComplex G9 = s_track[4] * G7;
    magmaFloatComplex G10 = s_track[3] * G7;
    magmaFloatComplex G11 = s_track[2] * G7;
    magmaFloatComplex G12 = s_track[1] * G7;
    magmaFloatComplex G13 = s_track[0] + s_track[0];
    magmaFloatComplex G14 = G4 * G13;
    magmaFloatComplex G15 = G1 * r_startCoefs[2];
    magmaFloatComplex G16 = t * r_targetCoefs[2];
    magmaFloatComplex G17 = G15 + G16;
    magmaFloatComplex G18 = G14 + G17;
    magmaFloatComplex G19 = s_track[6] * G7;
    magmaFloatComplex G20 = G9 + G19;
    magmaFloatComplex G21 = G10 + G8;
    magmaFloatComplex G22 = G11 + G9;
    magmaFloatComplex G23 = s_track[1] + s_track[1];
    magmaFloatComplex G24 = G4 * G23;
    magmaFloatComplex G25 = G24 + G10;
    magmaFloatComplex G26 = G7 * s_track[0];
    magmaFloatComplex G27 = G26 + G11;
    magmaFloatComplex G28 = G27 + G17;
    magmaFloatComplex G29 = G7 * G23;
    magmaFloatComplex G30 = s_track[2] + s_track[2];
    magmaFloatComplex G31 = G4 * G30;
    magmaFloatComplex G32 = G31 + G19;
    magmaFloatComplex G33 = G7 * s_track[1];
    magmaFloatComplex G34 = G33 + G8;
    magmaFloatComplex G35 = G26 + G9;
    magmaFloatComplex G36 = G35 + G17;
    magmaFloatComplex G37 = G33 + G10;
    magmaFloatComplex G38 = G7 * G30;
    magmaFloatComplex G39 = G7 * s_track[2];
    magmaFloatComplex G40 = G26 + G19;
    magmaFloatComplex G41 = G40 + G17;
    magmaFloatComplex G42 = G39 + G9;
    magmaFloatComplex G43 = s_track[3] + s_track[3];
    magmaFloatComplex G44 = G7 * G43;
    magmaFloatComplex G45 = G26 + G17;
    magmaFloatComplex G46 = G39 + G19;
    magmaFloatComplex G47 = G7 * s_track[3];
    magmaFloatComplex G48 = G47 + G8;
    magmaFloatComplex G49 = s_track[4] + s_track[4];
    magmaFloatComplex G50 = G7 * G49;
    magmaFloatComplex G51 = G7 * s_track[4];
    magmaFloatComplex G52 = G51 + G19;
    magmaFloatComplex G53 = s_track[5] + s_track[5];
    magmaFloatComplex G54 = G7 * G53;
    magmaFloatComplex G55 = G7 * s_track[5];
    magmaFloatComplex G56 = s_track[6] + s_track[6];
    magmaFloatComplex G57 = G7 * G56;
    magmaFloatComplex G59 = r_targetCoefs[0] - r_startCoefs[0];
    magmaFloatComplex G60 = s_track[0] * G59;
    magmaFloatComplex G62 = r_targetCoefs[1] - r_startCoefs[1];
    magmaFloatComplex G63 = s_track[1] * G62;
    magmaFloatComplex G64 = G60 + G63;
    magmaFloatComplex G65 = s_track[2] * G62;
    magmaFloatComplex G66 = G64 + G65;
    magmaFloatComplex G67 = s_track[3] * G62;
    magmaFloatComplex G68 = G66 + G67;
    magmaFloatComplex G69 = s_track[4] * G62;
    magmaFloatComplex G70 = G68 + G69;
    magmaFloatComplex G71 = s_track[5] * G62;
    magmaFloatComplex G72 = G70 + G71;
    magmaFloatComplex G73 = s_track[6] * G62;
    magmaFloatComplex G74 = G72 + G73;
    magmaFloatComplex G76 = r_targetCoefs[2] - r_startCoefs[2];
    magmaFloatComplex G77 = G74 + G76;
    magmaFloatComplex G78 = s_track[0] * G62;
    magmaFloatComplex G79 = s_track[5] * G78;
    magmaFloatComplex G80 = s_track[4] * G63;
    magmaFloatComplex G81 = G79 + G80;
    magmaFloatComplex G82 = s_track[6] * G63;
    magmaFloatComplex G83 = G81 + G82;
    magmaFloatComplex G84 = s_track[3] * G65;
    magmaFloatComplex G85 = G83 + G84;
    magmaFloatComplex G86 = s_track[5] * G76;
    magmaFloatComplex G87 = G85 + G86;
    magmaFloatComplex G88 = s_track[4] * G78;
    magmaFloatComplex G89 = s_track[3] * G63;
    magmaFloatComplex G90 = G88 + G89;
    magmaFloatComplex G91 = s_track[5] * G63;
    magmaFloatComplex G92 = G90 + G91;
    magmaFloatComplex G93 = s_track[2] * s_track[2];
    magmaFloatComplex G94 = G93 * G59;
    magmaFloatComplex G95 = G92 + G94;
    magmaFloatComplex G96 = s_track[6] * G65;
    magmaFloatComplex G97 = G95 + G96;
    magmaFloatComplex G98 = s_track[4] * G76;
    magmaFloatComplex G99 = G97 + G98;
    magmaFloatComplex G100 = s_track[3] * G78;
    magmaFloatComplex G101 = s_track[2] * G63;
    magmaFloatComplex G102 = G100 + G101;
    magmaFloatComplex G103 = G102 + G80;
    magmaFloatComplex G104 = s_track[5] * G65;
    magmaFloatComplex G105 = G103 + G104;
    magmaFloatComplex G106 = s_track[6] * G67;
    magmaFloatComplex G107 = G105 + G106;
    magmaFloatComplex G108 = s_track[3] * G76;
    magmaFloatComplex G109 = G107 + G108;
    magmaFloatComplex G110 = s_track[2] * G78;
    magmaFloatComplex G111 = s_track[1] * s_track[1];
    magmaFloatComplex G112 = G111 * G59;
    magmaFloatComplex G113 = G110 + G112;
    magmaFloatComplex G114 = G113 + G89;
    magmaFloatComplex G115 = s_track[4] * G65;
    magmaFloatComplex G116 = G114 + G115;
    magmaFloatComplex G117 = s_track[2] * G76;
    magmaFloatComplex G118 = G116 + G117;
    magmaFloatComplex G119 = s_track[5] * G67;
    magmaFloatComplex G120 = G118 + G119;
    magmaFloatComplex G121 = s_track[6] * G69;
    magmaFloatComplex G122 = G120 + G121;
    magmaFloatComplex G123 = s_track[1] * G78;
    magmaFloatComplex G124 = G123 + G101;
    magmaFloatComplex G125 = s_track[1] * G76;
    magmaFloatComplex G126 = G124 + G125;
    magmaFloatComplex G127 = G126 + G84;
    magmaFloatComplex G128 = s_track[4] * G67;
    magmaFloatComplex G129 = G127 + G128;
    magmaFloatComplex G130 = s_track[5] * G69;
    magmaFloatComplex G131 = G129 + G130;
    magmaFloatComplex G132 = s_track[6] * G71;
    magmaFloatComplex G133 = G131 + G132;
    magmaFloatComplex G134 = s_track[0] * s_track[0];
    magmaFloatComplex G135 = G134 * G59;
    magmaFloatComplex G136 = s_track[0] * G76;
    magmaFloatComplex G137 = G135 + G136;
    magmaFloatComplex G138 = G111 * G62;
    magmaFloatComplex G139 = G137 + G138;
    magmaFloatComplex G140 = G93 * G62;
    magmaFloatComplex G141 = G139 + G140;
    magmaFloatComplex G142 = s_track[3] * s_track[3];
    magmaFloatComplex G143 = G142 * G62;
    magmaFloatComplex G144 = G141 + G143;
    magmaFloatComplex G145 = s_track[4] * s_track[4];
    magmaFloatComplex G146 = G145 * G62;
    magmaFloatComplex G147 = G144 + G146;
    magmaFloatComplex G148 = s_track[5] * s_track[5];
    magmaFloatComplex G149 = G148 * G62;
    magmaFloatComplex G150 = G147 + G149;
    magmaFloatComplex G151 = s_track[6] * s_track[6];
    magmaFloatComplex G152 = G151 * G62;
    magmaFloatComplex G153 = G150 + G152;

    r_cgesvA[0] = G4;
    r_cgesvA[1] = G8;
    r_cgesvA[2] = G9;
    r_cgesvA[3] = G10;
    r_cgesvA[4] = G11;
    r_cgesvA[5] = G12;
    r_cgesvA[6] = G18;
    r_cgesvB[0] = -G77;

    r_cgesvA[7] = G7;
    r_cgesvA[8] = G20;
    r_cgesvA[9] = G21;
    r_cgesvA[10] = G22;
    r_cgesvA[11] = G25;
    r_cgesvA[12] = G28;
    r_cgesvA[13] = G29;
    r_cgesvB[1] = -G87;

    r_cgesvA[14] = G7;
    r_cgesvA[15] = G10;
    r_cgesvA[16] = G32;
    r_cgesvA[17] = G34;
    r_cgesvA[18] = G36;
    r_cgesvA[19] = G37;
    r_cgesvA[20] = G38;
    r_cgesvB[2] = -G99;

    r_cgesvA[21] = G7;
    r_cgesvA[22] = G39;
    r_cgesvA[23] = G33;
    r_cgesvA[24] = G41;
    r_cgesvA[25] = G34;
    r_cgesvA[26] = G42;
    r_cgesvA[27] = G44;
    r_cgesvB[3] = -G109;

    r_cgesvA[28] = G7;
    r_cgesvA[29] = G33;
    r_cgesvA[30] = G45;
    r_cgesvA[31] = G33;
    r_cgesvA[32] = G46;
    r_cgesvA[33] = G48;
    r_cgesvA[34] = G50;
    r_cgesvB[4] = -G122;

    r_cgesvA[35] = G7;
    r_cgesvA[36] = G45;
    r_cgesvA[37] = G33;
    r_cgesvA[38] = G39;
    r_cgesvA[39] = G47;
    r_cgesvA[40] = G52;
    r_cgesvA[41] = G54;
    r_cgesvB[5] = -G133;

    r_cgesvA[42] = G7;
    r_cgesvA[43] = G33;
    r_cgesvA[44] = G39;
    r_cgesvA[45] = G47;
    r_cgesvA[46] = G51;
    r_cgesvA[47] = G55;
    r_cgesvA[48] = G57;
    r_cgesvB[6] = -G153;

  }
}

#endif
