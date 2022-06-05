#ifndef cpu_eval_HxH_katsura6_h
#define cpu_eval_HxH_katsura6_h
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
  void cpu_eval_HxH_katsura6(
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
    magmaFloatComplex G58 = G4 * s_track[0];
    magmaFloatComplex G59 = G58 + G33;
    magmaFloatComplex G60 = G59 + G39;
    magmaFloatComplex G61 = G60 + G47;
    magmaFloatComplex G62 = G61 + G51;
    magmaFloatComplex G63 = G62 + G55;
    magmaFloatComplex G64 = G7 * s_track[6];
    magmaFloatComplex G65 = G63 + G64;
    magmaFloatComplex G66 = G65 + G17;
    magmaFloatComplex G67 = G26 * s_track[5];
    magmaFloatComplex G68 = G33 * s_track[4];
    magmaFloatComplex G69 = G67 + G68;
    magmaFloatComplex G70 = G33 * s_track[6];
    magmaFloatComplex G71 = G69 + G70;
    magmaFloatComplex G72 = G39 * s_track[3];
    magmaFloatComplex G73 = G71 + G72;
    magmaFloatComplex G74 = G17 * s_track[5];
    magmaFloatComplex G75 = G73 + G74;
    magmaFloatComplex G76 = G26 * s_track[4];
    magmaFloatComplex G77 = G33 * s_track[3];
    magmaFloatComplex G78 = G76 + G77;
    magmaFloatComplex G79 = G33 * s_track[5];
    magmaFloatComplex G80 = G78 + G79;
    magmaFloatComplex G81 = s_track[2] * s_track[2];
    magmaFloatComplex G82 = G4 * G81;
    magmaFloatComplex G83 = G80 + G82;
    magmaFloatComplex G84 = G39 * s_track[6];
    magmaFloatComplex G85 = G83 + G84;
    magmaFloatComplex G86 = G17 * s_track[4];
    magmaFloatComplex G87 = G85 + G86;
    magmaFloatComplex G88 = G26 * s_track[3];
    magmaFloatComplex G89 = G33 * s_track[2];
    magmaFloatComplex G90 = G88 + G89;
    magmaFloatComplex G91 = G90 + G68;
    magmaFloatComplex G92 = G39 * s_track[5];
    magmaFloatComplex G93 = G91 + G92;
    magmaFloatComplex G94 = G47 * s_track[6];
    magmaFloatComplex G95 = G93 + G94;
    magmaFloatComplex G96 = G17 * s_track[3];
    magmaFloatComplex G97 = G95 + G96;
    magmaFloatComplex G98 = G26 * s_track[2];
    magmaFloatComplex G99 = s_track[1] * s_track[1];
    magmaFloatComplex G100 = G4 * G99;
    magmaFloatComplex G101 = G98 + G100;
    magmaFloatComplex G102 = G101 + G77;
    magmaFloatComplex G103 = G39 * s_track[4];
    magmaFloatComplex G104 = G102 + G103;
    magmaFloatComplex G105 = G17 * s_track[2];
    magmaFloatComplex G106 = G104 + G105;
    magmaFloatComplex G107 = G47 * s_track[5];
    magmaFloatComplex G108 = G106 + G107;
    magmaFloatComplex G109 = G51 * s_track[6];
    magmaFloatComplex G110 = G108 + G109;
    magmaFloatComplex G111 = G26 * s_track[1];
    magmaFloatComplex G112 = G111 + G89;
    magmaFloatComplex G113 = G17 * s_track[1];
    magmaFloatComplex G114 = G112 + G113;
    magmaFloatComplex G115 = G114 + G72;
    magmaFloatComplex G116 = G47 * s_track[4];
    magmaFloatComplex G117 = G115 + G116;
    magmaFloatComplex G118 = G51 * s_track[5];
    magmaFloatComplex G119 = G117 + G118;
    magmaFloatComplex G120 = G55 * s_track[6];
    magmaFloatComplex G121 = G119 + G120;
    magmaFloatComplex G122 = s_track[0] * s_track[0];
    magmaFloatComplex G123 = G4 * G122;
    magmaFloatComplex G124 = G17 * s_track[0];
    magmaFloatComplex G125 = G123 + G124;
    magmaFloatComplex G126 = G7 * G99;
    magmaFloatComplex G127 = G125 + G126;
    magmaFloatComplex G128 = G7 * G81;
    magmaFloatComplex G129 = G127 + G128;
    magmaFloatComplex G130 = s_track[3] * s_track[3];
    magmaFloatComplex G131 = G7 * G130;
    magmaFloatComplex G132 = G129 + G131;
    magmaFloatComplex G133 = s_track[4] * s_track[4];
    magmaFloatComplex G134 = G7 * G133;
    magmaFloatComplex G135 = G132 + G134;
    magmaFloatComplex G136 = s_track[5] * s_track[5];
    magmaFloatComplex G137 = G7 * G136;
    magmaFloatComplex G138 = G135 + G137;
    magmaFloatComplex G139 = s_track[6] * s_track[6];
    magmaFloatComplex G140 = G7 * G139;
    magmaFloatComplex G141 = G138 + G140;

    r_cgesvA[0] = G4;
    r_cgesvA[1] = G8;
    r_cgesvA[2] = G9;
    r_cgesvA[3] = G10;
    r_cgesvA[4] = G11;
    r_cgesvA[5] = G12;
    r_cgesvA[6] = G18;
    r_cgesvB[0] = G66;

    r_cgesvA[7] = G7;
    r_cgesvA[8] = G20;
    r_cgesvA[9] = G21;
    r_cgesvA[10] = G22;
    r_cgesvA[11] = G25;
    r_cgesvA[12] = G28;
    r_cgesvA[13] = G29;
    r_cgesvB[1] = G75;

    r_cgesvA[14] = G7;
    r_cgesvA[15] = G10;
    r_cgesvA[16] = G32;
    r_cgesvA[17] = G34;
    r_cgesvA[18] = G36;
    r_cgesvA[19] = G37;
    r_cgesvA[20] = G38;
    r_cgesvB[2] = G87;

    r_cgesvA[21] = G7;
    r_cgesvA[22] = G39;
    r_cgesvA[23] = G33;
    r_cgesvA[24] = G41;
    r_cgesvA[25] = G34;
    r_cgesvA[26] = G42;
    r_cgesvA[27] = G44;
    r_cgesvB[3] = G97;

    r_cgesvA[28] = G7;
    r_cgesvA[29] = G33;
    r_cgesvA[30] = G45;
    r_cgesvA[31] = G33;
    r_cgesvA[32] = G46;
    r_cgesvA[33] = G48;
    r_cgesvA[34] = G50;
    r_cgesvB[4] = G110;

    r_cgesvA[35] = G7;
    r_cgesvA[36] = G45;
    r_cgesvA[37] = G33;
    r_cgesvA[38] = G39;
    r_cgesvA[39] = G47;
    r_cgesvA[40] = G52;
    r_cgesvA[41] = G54;
    r_cgesvB[5] = G121;

    r_cgesvA[42] = G7;
    r_cgesvA[43] = G33;
    r_cgesvA[44] = G39;
    r_cgesvA[45] = G47;
    r_cgesvA[46] = G51;
    r_cgesvA[47] = G55;
    r_cgesvA[48] = G57;
    r_cgesvB[6] = G141;
 
  }
}

#endif
