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
#include "gpu-kernels/magmaHC-kernels.hpp"

class GPU_HC_Solver {
    
    magma_device_t cdev;       // variable to indicate current gpu id
    magma_queue_t my_queue;    // magma queue variable, internally holds a cuda stream and a cublas handle

    //> Varaibles as sizes of arrays
    magma_int_t             ldda, lddb, ldd_params;
    magma_int_t             dHdx_Index_Size;
    magma_int_t             dHdt_Index_Size;
    magma_int_t             dHdx_PHC_Coeffs_Size;
    magma_int_t             dHdt_PHC_Coeffs_Size;
    magma_int_t             ldd_phc_Params_Hx;
    magma_int_t             ldd_phc_Params_Ht;

    //> Variables and arrays on the CPU side
    magmaFloatComplex       *h_GPU_HC_Track_Sols;
    magmaFloatComplex       *h_diffParams;
    magmaFloatComplex       *h_Start_Sols;
    magmaFloatComplex       *h_Homotopy_Sols;
    magmaFloatComplex       *h_Start_Params;
    magmaFloatComplex       *h_Target_Params;
    magmaFloatComplex       *h_dHdx_PHC_Coeffs;
    magmaFloatComplex       *h_dHdt_PHC_Coeffs;
    magma_int_t             *h_dHdx_Index;
    magma_int_t             *h_dHdt_Index;
    magmaFloatComplex       *h_Debug_Purpose;
    bool                    *h_is_GPU_HC_Sol_Converge;
    bool                    *h_is_GPU_HC_Sol_Infinity;

    //> Variables and arrays on the GPU side
    magmaFloatComplex_ptr   d_Start_Sols, d_Homotopy_Sols;
    magmaFloatComplex_ptr   d_Start_Params, d_Target_Params;
    magmaFloatComplex_ptr   d_dHdx_PHC_Coeffs;
    magmaFloatComplex_ptr   d_dHdt_PHC_Coeffs;
    magma_int_t             *d_dHdx_Index;
    magma_int_t             *d_dHdt_Index;
    magmaFloatComplex       **d_Start_Sols_array;
    magmaFloatComplex       **d_Homotopy_Sols_array;
    magmaFloatComplex       *d_diffParams;
    magmaFloatComplex       *d_Debug_Purpose;
    bool                    *d_is_GPU_HC_Sol_Converge;
    bool                    *d_is_GPU_HC_Sol_Infinity;

public:
    //> The timers
    real_Double_t           gpu_time;
    real_Double_t           transfer_h2d_time;
    real_Double_t           transfer_d2h_time;

    //> Constructor
    GPU_HC_Solver();
    
    //> Member functions
    bool Read_Problem_Data();
    void Allocate_Arrays();
    void Construct_Coeffs_From_Params();
    void Data_Transfer_From_Host_To_Device();
    void Solve_by_GPU_HC();

    //> Destructor
    ~GPU_HC_Solver();

private:
    std::string Problem_File_Path;

};

#endif
