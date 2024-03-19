#ifndef magmaHC_dev_sm_cuh_
#define magmaHC_dev_sm_cuh_
// ============================================================================
// Device function for magma cgesv batched small
//
// Modifications
//    Chiang-Heng Chien  21-05-03:   Include device function from cgesv_batched_small.cu
//
//> (c) LEMS, Brown University
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// ============================================================================
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>

// -- cuda included --
#include <cuda_runtime.h>

// -- magma included --
#include "flops.h"
#include "magma_v2.h"
#include "magma_lapack.h"
#include "magma_internal.h"
#undef max
#undef min
#include "magma_templates.h"
#include "sync.cuh"
#undef max
#undef min
#include "shuffle.cuh"
#undef max
#undef min
#include "batched_kernel_param.h"

#define sA(i,j)  sA[(j)*slda + (i)]
#define sB(i,j)  sB[(j)*sldb + (i)]

__device__ __inline__ void
cgesv_batched_small_sm_device(
    const int tx,
    int n, int nrhs,
    magmaFloatComplex *sA, int slda, int *sipiv,
    magmaFloatComplex *sB, int sldb,
    magmaFloatComplex *sx, float *dsx,
    int &linfo )
{
    int max_id;
    float rx_abs_max = MAGMA_D_ZERO;
    magmaFloatComplex reg    = MAGMA_C_ZERO;
    magmaFloatComplex update = MAGMA_C_ZERO;

    linfo = 0;

    #pragma unroll
    for(int i = 0; i < n; i++) {
        // icamax and find pivot
        dsx[ tx ] = fabs(MAGMA_C_REAL( sA(tx,i) )) + fabs(MAGMA_C_IMAG( sA(tx,i) ));
        __syncthreads();
        rx_abs_max = dsx[i];
        max_id = i;
        for(int j = i+1; j < n; j++){
            if( dsx[j] > rx_abs_max){
                max_id = j;
                rx_abs_max = dsx[j];
            }
        }
        linfo  = ( rx_abs_max == MAGMA_D_ZERO && linfo == 0) ? (i+1) : linfo;
        update = ( rx_abs_max == MAGMA_D_ZERO ) ? MAGMA_C_ZERO : MAGMA_C_ONE;

        // write pivot index
        if(tx == 0){
            sipiv[i] = max_id;
        }

        // swap
        if( max_id != i) {
            reg            = sA(i, tx);
            sA(i, tx)      = sA(max_id, tx);
            sA(max_id, tx) = reg;

            for(int itx = tx; itx < nrhs; itx+=blockDim.x) {
                reg             = sB(i, itx);
                sB(i, itx)      = sB(max_id, itx);
                sB(max_id, itx) = reg;
            }
        }
        __syncthreads();

        reg = ( rx_abs_max == MAGMA_D_ZERO ) ? MAGMA_C_ONE : MAGMA_C_DIV(MAGMA_C_ONE, sA(i,i) );
        // scal and ger
        if( tx > i ){
            sA(tx,i) *= reg;
            for(int j = i+1; j < n; j++) {
                sA(tx, j) -= sA(tx, i) * ( update * sA(i, j) );
            }

            for(int j = 0; j < nrhs; j++) {
                sB(tx, j) -= sA(tx, i) * ( update * sB(i, j) );
            }
        }
        __syncthreads();
    }

    for(int i = n-1; i >= 0; i--) {
        for(int j = 0; j < nrhs; j++) {
            reg       = MAGMA_C_DIV(sB(i, j), sA(i,i));
            __syncthreads();
            sB(tx, j) = (tx <  i) ? sB(tx, j) - reg * sA(tx,i): sB(tx, j);
            sB(tx, j) = (tx == i) ? reg : sB(tx, j);
            __syncthreads();
        }
    }
}

#undef sA
#undef sB

#endif //magmaHC_dev_sm_cuh_
