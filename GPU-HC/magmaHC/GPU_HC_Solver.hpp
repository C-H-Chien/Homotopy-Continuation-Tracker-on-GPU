#ifndef magmaHC_kernels_h
#define magmaHC_kernels_h
// ============================================================================
// Header file declaring all kernels
//
// Modifications
//    Chiang-Heng Chien  22-10-31:   Initially Created (Copied from other repos)
//
//> (c) LEMS, Brown University
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// ============================================================================
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <chrono>

//> MAGMA
#include "magma_v2.h"
#include "Evaluations.hpp"
#include "gpu-kernels/magmaHC-kernels.hpp"

class GPU_HC_Solver {
    
    magma_device_t cdev;       // variable to indicate current gpu id
    magma_queue_t my_queue;    // magma queue variable, internally holds a cuda stream and a cublas handle

    //> Varaibles as sizes of arrays
    magma_int_t             ldda, lddb, ldd_params;
    magma_int_t             size_Hx;
    magma_int_t             size_Ht;
    magma_int_t             phc_coeffs_Hx_size;
    magma_int_t             phc_coeffs_Ht_size;
    magma_int_t             ldd_phc_Params_Hx;
    magma_int_t             ldd_phc_Params_Ht;

    //> Variables and arrays on the CPU side
    magmaFloatComplex       *h_GPU_HC_Track_Sols;
    magmaFloatComplex       *h_params_diff;
    bool                    *h_is_GPU_HC_Sol_Converge;
    bool                    *h_is_GPU_HC_Sol_Infinity;
    
    magmaFloatComplex       *h_startSols;
    magmaFloatComplex       *h_Track;
    magmaFloatComplex       *h_startParams;
    magmaFloatComplex       *h_targetParams;
    magma_int_t             *h_Hx_idx;
    magma_int_t             *h_Ht_idx;
    magmaFloatComplex       *h_phc_coeffs_Hx;
    magmaFloatComplex       *h_phc_coeffs_Ht;

    //> Variables and arrays on the GPU side
    bool                    *d_is_GPU_HC_Sol_Converge;
    bool                    *d_is_GPU_HC_Sol_Infinity;
    magmaFloatComplex_ptr   d_startSols, d_Track;
    magmaFloatComplex_ptr   d_startParams, d_targetParams;
    magmaFloatComplex_ptr   d_phc_coeffs_Hx;
    magmaFloatComplex_ptr   d_phc_coeffs_Ht;
    magma_int_t             *d_Hx_idx;
    magma_int_t             *d_Ht_idx;
    magmaFloatComplex       **d_startSols_array;
    magmaFloatComplex       **d_Track_array;
    magmaFloatComplex       *d_diffParams;

    magmaFloatComplex       *h_Debug_Purpose;
    magmaFloatComplex       *d_Debug_Purpose;

public:
    //> The timers
    real_Double_t           gpu_time;
    real_Double_t           data_h2d_time, data_d2h_time;
    
    GPU_HC_Solver( magmaFloatComplex*, magmaFloatComplex*, \
                   magmaFloatComplex*, magmaFloatComplex*, \
                   magma_int_t*, magma_int_t*, \
                   magmaFloatComplex*, magmaFloatComplex* );
    
    //> Member functions
    //void Prepare_Files_for_Write();
    void Allocate_Arrays();
    void Data_Transfer_From_Host_To_Device();
    void Solve_by_GPU_HC();

    // friend class Evaluations;

    ~GPU_HC_Solver();
};

#endif
