#ifndef DEV_TRIFOCAL_2OP1P_EVAL_CUH_
#define DEV_TRIFOCAL_2OP1P_EVAL_CUH_
// ==============================================================================================================
// Device function for evaluating the covnerged, real solutions of trifocal pose estimation problem
//
// Modifications
//    Chien   24-08-17:   Initially created for the performance improvement of GPU-HC on trifocal pose estimation
//
//> (c) LEMS, Brown University
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// ==============================================================================================================
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>

#include <cuda_runtime.h>

#include "magma_v2.h"
#include "magma_lapack.h"
#include "magma_internal.h"

#include "definitions.hpp"

//> Macros
#define d_Triplet_Edges(i,j)    d_Triplet_Edges[(i) * 6 + (j)]

template< int Num_Of_Vars >
__device__ __inline__ void
evaluate_trifocal_2op1p_30x30_sol(
  const int tx, const int batchid, \
  const float* __restrict__ d_Triplet_Edge_Locations, \
  const int Num_Of_Edgels, const int Num_Of_Edgels_Process_Rounds, \
  const int Num_Of_Remaining_Edgels, const int Edgels_Index_Offset, \
  int *s_RK_Coeffs, magmaFloatComplex *s_track, float *dsx, \
  float *s_hypothesis_Rs, float *s_intrinsic_matrix, \
  int *s_local_num_of_inlier_supports_21, int *s_local_num_of_inlier_supports_31, \
  float one_half_delta_t, float r_sqrt_sols, float r_sqrt_corr, \
  bool* &d_Found_Trifocal_Sols, int* &d_Trifocal_Sols_Batch_Index )
{
    if (tx == 0) {
      //> reset the value of s_RK_Coeffs[0]
      s_RK_Coeffs[0] = 0;

      //> Check the imaginary part of the two translations and two relative rotations
      if (fabs(MAGMA_C_IMAG(s_track[18])) < IMAG_PART_TOL && fabs(MAGMA_C_IMAG(s_track[19])) < IMAG_PART_TOL && \
          fabs(MAGMA_C_IMAG(s_track[20])) < IMAG_PART_TOL && fabs(MAGMA_C_IMAG(s_track[21])) < IMAG_PART_TOL && \
          fabs(MAGMA_C_IMAG(s_track[22])) < IMAG_PART_TOL && fabs(MAGMA_C_IMAG(s_track[23])) < IMAG_PART_TOL && \
          fabs(MAGMA_C_IMAG(s_track[24])) < IMAG_PART_TOL && fabs(MAGMA_C_IMAG(s_track[25])) < IMAG_PART_TOL && \
          fabs(MAGMA_C_IMAG(s_track[26])) < IMAG_PART_TOL && fabs(MAGMA_C_IMAG(s_track[27])) < IMAG_PART_TOL && \
          fabs(MAGMA_C_IMAG(s_track[28])) < IMAG_PART_TOL && fabs(MAGMA_C_IMAG(s_track[29])) < IMAG_PART_TOL) {
        s_RK_Coeffs[0] = 1;
      }
    }
    magmablas_syncwarp();

    //> Continue if the previous step passed
    if (s_RK_Coeffs[0] == 1) {

      //> Reusing dsx for real parts of the solution
      dsx[ tx ] = MAGMA_C_REAL(s_track[ tx ]);

      //> Use 2 threads to get normalized rotation matrices
      if (tx < 2) {
        
        //> Translation vector can be up to a scale, but rotation has to be normalized. The first and second 9 entries are R1 and R2, respectively.
        s_hypothesis_Rs[ tx*9 + 0 ] = 1 + dsx[tx*3 + 24]*dsx[tx*3 + 24] - (dsx[tx*3 + 25]*dsx[tx*3 + 25] + dsx[tx*3 + 26]*dsx[tx*3 + 26]);
        s_hypothesis_Rs[ tx*9 + 1 ] = 2*(dsx[tx*3 + 24]*dsx[tx*3 + 25] - dsx[tx*3 + 26]);
        s_hypothesis_Rs[ tx*9 + 2 ] = 2*(dsx[tx*3 + 24]*dsx[tx*3 + 26] + dsx[tx*3 + 25]);
        s_hypothesis_Rs[ tx*9 + 3 ] = 2*(dsx[tx*3 + 24]*dsx[tx*3 + 25] + dsx[tx*3 + 26]);
        s_hypothesis_Rs[ tx*9 + 4 ] = 1 + dsx[tx*3 + 25]*dsx[tx*3 + 25] - (dsx[tx*3 + 24]*dsx[tx*3 + 24] + dsx[tx*3 + 26]*dsx[tx*3 + 26]);
        s_hypothesis_Rs[ tx*9 + 5 ] = 2*(dsx[tx*3 + 25]*dsx[tx*3 + 26] - dsx[tx*3 + 24]);
        s_hypothesis_Rs[ tx*9 + 6 ] = 2*(dsx[tx*3 + 24]*dsx[tx*3 + 26] - dsx[tx*3 + 25]);
        s_hypothesis_Rs[ tx*9 + 7 ] = 2*(dsx[tx*3 + 25]*dsx[tx*3 + 26] + dsx[tx*3 + 24]);
        s_hypothesis_Rs[ tx*9 + 8 ] = 1 + dsx[tx*3 + 26]*dsx[tx*3 + 26] - (dsx[tx*3 + 24]*dsx[tx*3 + 24] + dsx[tx*3 + 25]*dsx[tx*3 + 25]);

        //> Reusing one_half_delta_t, r_sqrt_sols, and r_sqrt_corr for norm of the three column vectors of the rotation matrix
        one_half_delta_t = rnorm3df(s_hypothesis_Rs[ tx*9 + 0 ], s_hypothesis_Rs[ tx*9 + 3 ], s_hypothesis_Rs[ tx*9 + 6 ]);
        r_sqrt_sols      = rnorm3df(s_hypothesis_Rs[ tx*9 + 1 ], s_hypothesis_Rs[ tx*9 + 4 ], s_hypothesis_Rs[ tx*9 + 7 ]);
        r_sqrt_corr      = rnorm3df(s_hypothesis_Rs[ tx*9 + 2 ], s_hypothesis_Rs[ tx*9 + 5 ], s_hypothesis_Rs[ tx*9 + 8 ]);

        //> Normalizing the rotation matrix
        s_hypothesis_Rs[ tx*9 + 0 ] *= one_half_delta_t;
        s_hypothesis_Rs[ tx*9 + 1 ] *= one_half_delta_t;
        s_hypothesis_Rs[ tx*9 + 2 ] *= one_half_delta_t;
        s_hypothesis_Rs[ tx*9 + 3 ] *= r_sqrt_sols;
        s_hypothesis_Rs[ tx*9 + 4 ] *= r_sqrt_sols;
        s_hypothesis_Rs[ tx*9 + 5 ] *= r_sqrt_sols;
        s_hypothesis_Rs[ tx*9 + 6 ] *= r_sqrt_corr;
        s_hypothesis_Rs[ tx*9 + 7 ] *= r_sqrt_corr;
        s_hypothesis_Rs[ tx*9 + 8 ] *= r_sqrt_corr;

        //> Normalizing translations. one_falg_delta_t is the norm of the translation vector.
        // one_half_delta_t = rnorm3df(dsx[ tx*3 + 0 ], dsx[ tx*3 + 1 ], dsx[ tx*3 + 2 ]);
        // dsx[tx*3 + 18] = fdividef( dsx[tx*3 + 18], one_half_delta_t );
        // dsx[tx*3 + 19] = fdividef( dsx[tx*3 + 19], one_half_delta_t );
        // dsx[tx*3 + 20] = fdividef( dsx[tx*3 + 20], one_half_delta_t );
      }
      magmablas_syncwarp();

      const float* __restrict__ d_Triplet_Edges = d_Triplet_Edge_Locations;
      float vec_3d[3] = {0.0};

      //> Read each edgel and validate the solution in parallel
      for (int i = 0; i < Num_Of_Edgels_Process_Rounds; i++) {
        
        // R21 = [s_hypothesis_Rs[0], s_hypothesis_Rs[1], s_hypothesis_Rs[2]; 
        //        s_hypothesis_Rs[3], s_hypothesis_Rs[4], s_hypothesis_Rs[5];
        //        s_hypothesis_Rs[6], s_hypothesis_Rs[7], s_hypothesis_Rs[8] ]
        // R31 = [s_hypothesis_Rs[9],  s_hypothesis_Rs[10], s_hypothesis_Rs[11]; 
        //        s_hypothesis_Rs[12], s_hypothesis_Rs[13], s_hypothesis_Rs[14];
        //        s_hypothesis_Rs[15], s_hypothesis_Rs[16], s_hypothesis_Rs[17] ]
        // T21 = (dsx[18], dsx[19], dsx[20]), T31 = (dsx[21], dsx[22], dsx[23])
        // g1 = [d_Triplet_Edges(tx + Num_Of_Vars*(i), 0); d_Triplet_Edges(tx + Num_Of_Vars*(i), 1); 1]
        // g2 = [d_Triplet_Edges(tx + Num_Of_Vars*(i), 2); d_Triplet_Edges(tx + Num_Of_Vars*(i), 3); 1]
        // g3 = [d_Triplet_Edges(tx + Num_Of_Vars*(i), 4); d_Triplet_Edges(tx + Num_Of_Vars*(i), 5); 1]
        // fx = s_intrinsic_matrix[0]
        // fy = s_intrinsic_matrix[1]
        // cx = s_intrinsic_matrix[2]
        // cy = s_intrinsic_matrix[3]

        //> ================================================================
        //> From image 1 to 2
        //> ================================================================
        //> Reusing the register file again
        //> r_sqrt_sols is the numerator of rho equation (rho_numer)
        //> r_sqrt_corr is the denominator of the rho equation (rho_denom)
        //> one_half_delta_t is the rho
        r_sqrt_sols = dsx[20]*(s_hypothesis_Rs[2]*d_Triplet_Edges(tx + Num_Of_Vars*(i), 2) + s_hypothesis_Rs[5]*d_Triplet_Edges(tx + Num_Of_Vars*(i), 3) + s_hypothesis_Rs[8]) \
                      - (s_hypothesis_Rs[2]*dsx[18] + s_hypothesis_Rs[5]*dsx[19] + s_hypothesis_Rs[8]*dsx[20]);
        r_sqrt_corr = 1 - (s_hypothesis_Rs[6]*d_Triplet_Edges(tx + Num_Of_Vars*(i), 0) + s_hypothesis_Rs[7]*d_Triplet_Edges(tx + Num_Of_Vars*(i), 1) + s_hypothesis_Rs[8]) \
                          *(s_hypothesis_Rs[2]*d_Triplet_Edges(tx + Num_Of_Vars*(i), 2) + s_hypothesis_Rs[5]*d_Triplet_Edges(tx + Num_Of_Vars*(i), 3) + s_hypothesis_Rs[8]);

        //> Now vec_3d is the projected point in meters from rho_numer*R21*gamma1 + rho_denom*T21
        vec_3d[2] = r_sqrt_sols*(s_hypothesis_Rs[6]*d_Triplet_Edges(tx + Num_Of_Vars*(i), 0) + s_hypothesis_Rs[7]*d_Triplet_Edges(tx + Num_Of_Vars*(i), 1) + s_hypothesis_Rs[8]) + r_sqrt_corr*dsx[20];
        vec_3d[0] = (r_sqrt_sols*(s_hypothesis_Rs[0]*d_Triplet_Edges(tx + Num_Of_Vars*(i), 0) + s_hypothesis_Rs[1]*d_Triplet_Edges(tx + Num_Of_Vars*(i), 1) + s_hypothesis_Rs[2]) + r_sqrt_corr*dsx[18]) / vec_3d[2];
        vec_3d[1] = (r_sqrt_sols*(s_hypothesis_Rs[3]*d_Triplet_Edges(tx + Num_Of_Vars*(i), 0) + s_hypothesis_Rs[4]*d_Triplet_Edges(tx + Num_Of_Vars*(i), 1) + s_hypothesis_Rs[5]) + r_sqrt_corr*dsx[19]) / vec_3d[2];
        //> Now r_sqrt_sols and r_sqrt_corr are reprojection error in x and y
        r_sqrt_sols = (vec_3d[0] * s_intrinsic_matrix[0] + s_intrinsic_matrix[2]) - (d_Triplet_Edges(tx + Num_Of_Vars*(i), 2) * s_intrinsic_matrix[0] + s_intrinsic_matrix[2]);
        r_sqrt_corr = (vec_3d[1] * s_intrinsic_matrix[1] + s_intrinsic_matrix[3]) - (d_Triplet_Edges(tx + Num_Of_Vars*(i), 3) * s_intrinsic_matrix[1] + s_intrinsic_matrix[3]);
        one_half_delta_t = hypotf(r_sqrt_sols, r_sqrt_corr);
        
        //> Add up the local value of s_local_num_of_inlier_supports_21 if it passes the reprojection error threshold
        if (one_half_delta_t < REPROJ_ERROR_INLIER_THRESH) 
          s_local_num_of_inlier_supports_21[tx]++;

        //> ================================================================
        //> From image 1 to 3
        //> ================================================================
        //> Reusing the register file again
        // rho = ((e3' * T31)*(e3' * R31' * g3) - e3'*R31'*T31) / (1 - (e3'*R31*g1)*(e3'*R31'*g3))
        // one_half_delta_t;
        //> r_sqrt_sols is the numerator of rho equation (rho_numer)
        //> r_sqrt_corr is the denominator of the rho equation (rho_denom)
        //> one_half_delta_t is the rho
        r_sqrt_sols = dsx[23]*(s_hypothesis_Rs[11]*d_Triplet_Edges(tx + Num_Of_Vars*(i), 4) + s_hypothesis_Rs[14]*d_Triplet_Edges(tx + Num_Of_Vars*(i), 5) + s_hypothesis_Rs[17]) \
                      - (s_hypothesis_Rs[11]*dsx[21] + s_hypothesis_Rs[14]*dsx[22] + s_hypothesis_Rs[17]*dsx[23]);
        r_sqrt_corr = 1 - (s_hypothesis_Rs[15]*d_Triplet_Edges(tx + Num_Of_Vars*(i), 0) + s_hypothesis_Rs[16]*d_Triplet_Edges(tx + Num_Of_Vars*(i), 1) + s_hypothesis_Rs[17]) \
                          *(s_hypothesis_Rs[11]*d_Triplet_Edges(tx + Num_Of_Vars*(i), 4) + s_hypothesis_Rs[14]*d_Triplet_Edges(tx + Num_Of_Vars*(i), 5) + s_hypothesis_Rs[17]);

        //> Now vec_3d is the projected point in meters from rho_numer*R31*gamma1 + rho_denom*T31
        vec_3d[2] = r_sqrt_sols*(s_hypothesis_Rs[15]*d_Triplet_Edges(tx + Num_Of_Vars*(i), 0) + s_hypothesis_Rs[16]*d_Triplet_Edges(tx + Num_Of_Vars*(i), 1) + s_hypothesis_Rs[17]) + r_sqrt_corr*dsx[23];
        vec_3d[0] = (r_sqrt_sols*(s_hypothesis_Rs[9]*d_Triplet_Edges(tx + Num_Of_Vars*(i), 0) + s_hypothesis_Rs[10]*d_Triplet_Edges(tx + Num_Of_Vars*(i), 1) + s_hypothesis_Rs[11]) + r_sqrt_corr*dsx[21]) / vec_3d[2];
        vec_3d[1] = (r_sqrt_sols*(s_hypothesis_Rs[12]*d_Triplet_Edges(tx + Num_Of_Vars*(i), 0) + s_hypothesis_Rs[13]*d_Triplet_Edges(tx + Num_Of_Vars*(i), 1) + s_hypothesis_Rs[14]) + r_sqrt_corr*dsx[22]) / vec_3d[2];
        //> Now r_sqrt_sols and r_sqrt_corr are reprojection error in x and y
        r_sqrt_sols = (vec_3d[0] * s_intrinsic_matrix[0] + s_intrinsic_matrix[2]) - (d_Triplet_Edges(tx + Num_Of_Vars*(i), 4) * s_intrinsic_matrix[0] + s_intrinsic_matrix[2]);
        r_sqrt_corr = (vec_3d[1] * s_intrinsic_matrix[1] + s_intrinsic_matrix[3]) - (d_Triplet_Edges(tx + Num_Of_Vars*(i), 5) * s_intrinsic_matrix[1] + s_intrinsic_matrix[3]);
        one_half_delta_t = hypotf(r_sqrt_sols, r_sqrt_corr);
        
        //> Add up the local value of s_local_num_of_inlier_supports_31 if it passes the reprojection error threshold
        if (one_half_delta_t < REPROJ_ERROR_INLIER_THRESH) 
          s_local_num_of_inlier_supports_31[tx]++;
      }
      magmablas_syncwarp();

      //> For the remaining edgels whose number is less than the number of variables, only partial threads are necessary to process
      if (tx < Num_Of_Remaining_Edgels) {

        //> ================================================================
        //> From image 1 to 2
        //> ================================================================
        //> Reusing the register file again
        //> r_sqrt_sols is the numerator of rho equation (rho_numer)
        //> r_sqrt_corr is the denominator of the rho equation (rho_denom)
        //> one_half_delta_t is the rho
        r_sqrt_sols = dsx[20]*(s_hypothesis_Rs[2]*d_Triplet_Edges(tx + Edgels_Index_Offset, 2) + s_hypothesis_Rs[5]*d_Triplet_Edges(tx + Edgels_Index_Offset, 3) + s_hypothesis_Rs[8]) \
                      - (s_hypothesis_Rs[2]*dsx[18] + s_hypothesis_Rs[5]*dsx[19] + s_hypothesis_Rs[8]*dsx[20]);
        r_sqrt_corr = 1 - (s_hypothesis_Rs[6]*d_Triplet_Edges(tx + Edgels_Index_Offset, 0) + s_hypothesis_Rs[7]*d_Triplet_Edges(tx + Edgels_Index_Offset, 1) + s_hypothesis_Rs[8]) \
                          *(s_hypothesis_Rs[2]*d_Triplet_Edges(tx + Edgels_Index_Offset, 2) + s_hypothesis_Rs[5]*d_Triplet_Edges(tx + Edgels_Index_Offset, 3) + s_hypothesis_Rs[8]);
        // one_half_delta_t = r_sqrt_sols / r_sqrt_corr; //> this is rho

        //> Now vec_3d is the projected point in meters from rho_numer*R21*gamma1 + rho_denom*T21
        vec_3d[2] = r_sqrt_sols*(s_hypothesis_Rs[6]*d_Triplet_Edges(tx + Edgels_Index_Offset, 0) + s_hypothesis_Rs[7]*d_Triplet_Edges(tx + Edgels_Index_Offset, 1) + s_hypothesis_Rs[8]) + r_sqrt_corr*dsx[20];
        vec_3d[0] = (r_sqrt_sols*(s_hypothesis_Rs[0]*d_Triplet_Edges(tx + Edgels_Index_Offset, 0) + s_hypothesis_Rs[1]*d_Triplet_Edges(tx + Edgels_Index_Offset, 1) + s_hypothesis_Rs[2]) + r_sqrt_corr*dsx[18]) / vec_3d[2];
        vec_3d[1] = (r_sqrt_sols*(s_hypothesis_Rs[3]*d_Triplet_Edges(tx + Edgels_Index_Offset, 0) + s_hypothesis_Rs[4]*d_Triplet_Edges(tx + Edgels_Index_Offset, 1) + s_hypothesis_Rs[5]) + r_sqrt_corr*dsx[19]) / vec_3d[2];
        //> Now r_sqrt_sols and r_sqrt_corr are reprojection error in x and y
        r_sqrt_sols = (vec_3d[0] * s_intrinsic_matrix[0] + s_intrinsic_matrix[2]) - (d_Triplet_Edges(tx + Edgels_Index_Offset, 2) * s_intrinsic_matrix[0] + s_intrinsic_matrix[2]);
        r_sqrt_corr = (vec_3d[1] * s_intrinsic_matrix[1] + s_intrinsic_matrix[3]) - (d_Triplet_Edges(tx + Edgels_Index_Offset, 3) * s_intrinsic_matrix[1] + s_intrinsic_matrix[3]);
        one_half_delta_t = hypotf(r_sqrt_sols, r_sqrt_corr);
        
        //> Add up the local value of s_local_num_of_inlier_supports_21 if it passes the reprojection error threshold
        if (one_half_delta_t < REPROJ_ERROR_INLIER_THRESH) 
          s_local_num_of_inlier_supports_21[tx]++;

        //> ================================================================
        //> From image 1 to 3
        //> ================================================================
        //> Reusing the register file again
        // rho = ((e3' * T31)*(e3' * R31' * g3) - e3'*R31'*T31) / (1 - (e3'*R31*g1)*(e3'*R31'*g3))
        // one_half_delta_t;
        //> r_sqrt_sols is the numerator of rho equation (rho_numer)
        //> r_sqrt_corr is the denominator of the rho equation (rho_denom)
        //> one_half_delta_t is the rho
        r_sqrt_sols = dsx[23]*(s_hypothesis_Rs[11]*d_Triplet_Edges(tx + Edgels_Index_Offset, 4) + s_hypothesis_Rs[14]*d_Triplet_Edges(tx + Edgels_Index_Offset, 5) + s_hypothesis_Rs[17]) \
                      - (s_hypothesis_Rs[11]*dsx[21] + s_hypothesis_Rs[14]*dsx[22] + s_hypothesis_Rs[17]*dsx[23]);
        r_sqrt_corr = 1 - (s_hypothesis_Rs[15]*d_Triplet_Edges(tx + Edgels_Index_Offset, 0) + s_hypothesis_Rs[16]*d_Triplet_Edges(tx + Edgels_Index_Offset, 1) + s_hypothesis_Rs[17]) \
                          *(s_hypothesis_Rs[11]*d_Triplet_Edges(tx + Edgels_Index_Offset, 4) + s_hypothesis_Rs[14]*d_Triplet_Edges(tx + Edgels_Index_Offset, 5) + s_hypothesis_Rs[17]);

        //> Now vec_3d is the projected point in meters from rho_numer*R31*gamma1 + rho_denom*T31
        vec_3d[2] = r_sqrt_sols*(s_hypothesis_Rs[15]*d_Triplet_Edges(tx + Edgels_Index_Offset, 0) + s_hypothesis_Rs[16]*d_Triplet_Edges(tx + Edgels_Index_Offset, 1) + s_hypothesis_Rs[17]) + r_sqrt_corr*dsx[23];
        vec_3d[0] = (r_sqrt_sols*(s_hypothesis_Rs[9]*d_Triplet_Edges(tx + Edgels_Index_Offset, 0) + s_hypothesis_Rs[10]*d_Triplet_Edges(tx + Edgels_Index_Offset, 1) + s_hypothesis_Rs[11]) + r_sqrt_corr*dsx[21]) / vec_3d[2];
        vec_3d[1] = (r_sqrt_sols*(s_hypothesis_Rs[12]*d_Triplet_Edges(tx + Edgels_Index_Offset, 0) + s_hypothesis_Rs[13]*d_Triplet_Edges(tx + Edgels_Index_Offset, 1) + s_hypothesis_Rs[14]) + r_sqrt_corr*dsx[22]) / vec_3d[2];
        //> Now r_sqrt_sols and r_sqrt_corr are reprojection error in x and y
        r_sqrt_sols = (vec_3d[0] * s_intrinsic_matrix[0] + s_intrinsic_matrix[2]) - (d_Triplet_Edges(tx + Edgels_Index_Offset, 4) * s_intrinsic_matrix[0] + s_intrinsic_matrix[2]);
        r_sqrt_corr = (vec_3d[1] * s_intrinsic_matrix[1] + s_intrinsic_matrix[3]) - (d_Triplet_Edges(tx + Edgels_Index_Offset, 5) * s_intrinsic_matrix[1] + s_intrinsic_matrix[3]);
        one_half_delta_t = hypotf(r_sqrt_sols, r_sqrt_corr);
        
        //> Add up the local value of s_local_num_of_inlier_supports_31 if it passes the reprojection error threshold
        if (one_half_delta_t < REPROJ_ERROR_INLIER_THRESH) 
          s_local_num_of_inlier_supports_31[tx]++;
      }
      magmablas_syncwarp();

      //> Collect all local number of inlier supports and return whether the converged GPU-HC solution passes the RANSAC inlier support test
      if (tx == 0) {
        for (int i = 1; i < Num_Of_Vars; i++) {
          s_local_num_of_inlier_supports_21[0] += s_local_num_of_inlier_supports_21[i];
          s_local_num_of_inlier_supports_31[0] += s_local_num_of_inlier_supports_31[i];
        }

        if ( fdividef( (float)s_local_num_of_inlier_supports_21[0], (float)Num_Of_Edgels ) >= PASS_RANSAC_INLIER_SUPPORT_RATIO && \
             fdividef( (float)s_local_num_of_inlier_supports_31[0], (float)Num_Of_Edgels ) >= PASS_RANSAC_INLIER_SUPPORT_RATIO ) {

          d_Found_Trifocal_Sols[0] = 1; 
          d_Trifocal_Sols_Batch_Index[batchid] = batchid;
        }
      }
      
    }
}

#endif
