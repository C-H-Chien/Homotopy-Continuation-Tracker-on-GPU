#ifndef UTIL_HPP
#define UTIL_HPP
// =============================================================================
//
// Modifications
//    Chiang-Heng Chien  23-07-14:   Intiailly Created for Multiview Geometry
//
//> (c) LEMS, Brown University
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// ==============================================================================
#include <cmath>
#include <math.h>
#include <algorithm>
#include <string.h>
#include <vector>

#include "definitions.hpp"

class util {

public:
    util() {
        T_x     = new float[9];
        E       = new float[9];
        F       = new float[9];
        inv_K   = new float[9];
    }

    void Cayley_To_Rotation_Matrix( float *r, float* &R_ ) {
        #define R_(i,j)  R_[(i)*3 + (j)]
        R_(0,0) = 1 + r[0]*r[0] - (r[1]*r[1] + r[2]*r[2]);
        R_(0,1) = 2*(r[0]*r[1] - r[2]);
        R_(0,2) = 2*(r[0]*r[2] + r[1]);
        R_(1,0) = 2*(r[0]*r[1] + r[2]);
        R_(1,1) = 1 + r[1]*r[1] - (r[0]*r[0] + r[2]*r[2]);
        R_(1,2) = 2*(r[1]*r[2] - r[0]);
        R_(2,0) = 2*(r[0]*r[2] - r[1]);
        R_(2,1) = 2*(r[1]*r[2] + r[0]);
        R_(2,2) = 1 + r[2]*r[2] - (r[0]*r[0] + r[1]*r[1]);

        //> The rotation matrix R now is up to some scale which needs to be normalized.
        R_ = Normalize_Rotation_Matrix( R_ );
    }

    float *Normalize_Rotation_Matrix( float *R ) {
        #define R(i,j)  R[(i)*3 + (j)]
        //> Compute the column vector norms
        float norm_col1 = sqrt(R(0,0)*R(0,0) + R(1,0)*R(1,0) + R(2,0)*R(2,0));
        float norm_col2 = sqrt(R(0,1)*R(0,1) + R(1,1)*R(1,1) + R(2,1)*R(2,1));
        float norm_col3 = sqrt(R(0,2)*R(0,2) + R(1,2)*R(1,2) + R(2,2)*R(2,2));
        R(0,0) /= norm_col1;
        R(1,0) /= norm_col1;
        R(2,0) /= norm_col1;
        R(0,1) /= norm_col2;
        R(1,1) /= norm_col2;
        R(2,1) /= norm_col2;
        R(0,2) /= norm_col3;
        R(1,2) /= norm_col3;
        R(2,2) /= norm_col3;

        //> Check whether R is normalized such that det(R)=1
        assert( fabs(get_Matrix_determinant(R) - 1.0) <= IS_SO3_DET_R_TOL );

        return R;
    }

    void Normalize_Translation_Vector( float* &transl ) {
        float norm_Transl = get_Vector_Norm(transl);
        transl[0] /= norm_Transl;
        transl[1] /= norm_Transl;
        transl[2] /= norm_Transl;
    }

    float get_Vector_Norm( float *vec ) {
        return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    }

    float get_Matrix_3x3_determinant(float *R) {
        return R(0,0)*R(1,1)*R(2,2) + R(0,1)*R(1,2)*R(2,0) + R(0,2)*R(1,0)*R(2,1) - \
                R(0,2)*R(1,1)*R(2,0) - R(0,1)*R(1,0)*R(2,2) - R(0,0)*R(1,2)*R(2,1);
    }

    template< int n >
    void get_Matrix_Matrix_Product(float* M1, float* M2, float* &M1M2) {
        #define M1(i,j)         M1[(i)*n + (j)]
        #define M2(i,j)         M2[(i)*n + (j)]
        #define M1M2(i,j)       M1M2[(i)*n + (j)]

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                M1M2(i,j) = 0;
                for (int k = 0; k < n; k++) {
                    M1M2(i,j) += M1(i,k) * M2(k,j);
                }
            }
        }
    }

    template< int n >
    void get_Matrix_Transpose(float* &R) {

        float transpose_R[ n*n ];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                transpose_R[(i)*n + (j)] = R[(j)*n + (i)];
            }
        }

        for (int i = 0; i < n*n; i++) R[i] = transpose_R[i];
    }

    template< int n >
    float get_Vector_Dot_Product(float* vec1, float* vec2) {
        float dot_prod = 0.0;
        for (int i = 0; i < n; i++) dot_prod += vec1[i] * vec2[i];
        return dot_prod;
    }

    template< int n >
    float get_Matrix_Trace(float* M) {
        float trace_ = 0.0;
        for (int i = 0; i < n; i++) trace_ += M[(i)*n + (i)];
        return trace_;
    }

    void get_Matrix_3x3_Inverse(float* M, float* &inv_M) {
        #define M(i,j)  M[(i)*3 + (j)]
        #define inv_M(i,j)  inv_M[(i)*3 + (j)]
        float det_of_M = get_Matrix_3x3_determinant(M);
        float inv_det = 1.0 / det_of_M;
        
        inv_M(0, 0) = (M(1, 1) * M(2, 2) - M(2, 1) * M(1, 2)) * inv_det;
        inv_M(0, 1) = (M(0, 2) * M(2, 1) - M(0, 1) * M(2, 2)) * inv_det;
        inv_M(0, 2) = (M(0, 1) * M(1, 2) - M(0, 2) * M(1, 1)) * inv_det;
        inv_M(1, 0) = (M(1, 2) * M(2, 0) - M(1, 0) * M(2, 2)) * inv_det;
        inv_M(1, 1) = (M(0, 0) * M(2, 2) - M(0, 2) * M(2, 0)) * inv_det;
        inv_M(1, 2) = (M(1, 0) * M(0, 2) - M(0, 0) * M(1, 2)) * inv_det;
        inv_M(2, 0) = (M(1, 0) * M(2, 1) - M(2, 0) * M(1, 1)) * inv_det;
        inv_M(2, 1) = (M(2, 0) * M(0, 1) - M(0, 0) * M(2, 1)) * inv_det;
        inv_M(2, 2) = (M(0, 0) * M(1, 1) - M(1, 0) * M(0, 1)) * inv_det;
    }

    void get_Skew_Symmetric_Matrix(float* T) {
        #define T_x(i,j)  T_x[(i)*3 + (j)]
        T_x(0,0) = 0.0;
        T_x(0,1) = -T[2];
        T_x(0,2) = T[1];
        T_x(1,0) = T[2];
        T_x(1,1) = 0.0;
        T_x(1,2) = -T[0];
        T_x(2,0) = -T[1];
        T_x(2,1) = T[0];
        T_x(2,2) = 0.0;
    }

    void get_Essential_Matrix( float* R21, float* T21 ) {
        //> E21 = (skew_T(T21)*R21);
        get_Skew_Symmetric_Matrix(T21);
        get_Matrix_Matrix_Product<3>(T_x, R21, E);
    }

    void get_Fundamental_Matrix( float* K, float* R21, float* T21 ) {
        //> First get the essential matrix
        get_Essential_Matrix(R21, T21);
        get_Matrix_3x3_Inverse(K, inv_K);

        get_Matrix_Transpose< 3 >(inv_K);

        get_Matrix_Matrix_Product< 3 >( inv_K, E, F);

        get_Matrix_Transpose< 3 >(inv_K);
        get_Matrix_Matrix_Product< 3 >( F, inv_K, F );
    }

    ~util() {
        delete [] T_x;
        delete [] E;
        delete [] F;
        delete [] inv_K;
    }

private:
    float *T_x;     //> skew symmetric matrix
    float *E;
    float *inv_K;
    float *F;
};



#endif
