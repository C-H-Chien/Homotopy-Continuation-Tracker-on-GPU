#ifndef CPUHC_GENERIC_SOLVER_EVAL_BY_INDX
#define CPUHC_GENERIC_SOLVER_EVAL_BY_INDX
// =======================================================================
// Solve homotopy continuation on cpu for the trifocal 2op1p 30x30 problem
//
// Major Modifications
//    Chiang-Heng Chien  22-09-29:   Initially Created
//
// =======================================================================
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <vector>
#include <omp.h>

#include "magma_v2.h"
#include "magma_lapack.h"
#include "../definitions.hpp"
#include "../CPU_HC_Solver.hpp"

real_Double_t CPU_HC_Solver::CPUHC_Generic_Solver_Eval_by_Indx(
    const int num_of_cores, 
    const Eval_dHdX& Eval_dHdX_func, 
    const Eval_dHdt& Eval_dHdt_func, 
    const Eval_H& Eval_H_func ) 
{
    CPU_HC_time = magma_wtime();

    // magma_int_t nthreads = num_of_cores;
    // omp_set_num_threads(nthreads);    
    magma_set_lapack_numthreads(1);
    magma_set_omp_numthreads(num_of_cores);

    // #pragma omp parallel for schedule(dynamic) collapse(2)
    #pragma omp parallel for schedule(dynamic) shared(h_Start_Params, h_diff_params, h_dHdx_Index, h_dHdt_Index)
    for (magma_int_t ix = 0; ix < NUM_OF_RANSAC_ITERATIONS * Num_Of_Tracks; ix++) {
    // for (magma_int_t ri = 0; ri < NUM_OF_RANSAC_ITERATIONS; ri++) {

        magma_int_t ri = ix / Num_Of_Tracks;
        magma_int_t s  = ix % Num_Of_Tracks;

        //> Loop over all homotopy paths
        // #pragma omp parallel for schedule(dynamic) shared(h_Start_Params, h_diff_params, h_dHdx_Index, h_dHdt_Index)
        // for (magma_int_t s = 0; s < Num_Of_Tracks; s++) {
            //> Local declarations for openmp parallelization to avoid race condition
            int pred_success_count = 0;
            float t0 = 0.0;
            float t_step = 0.0;
            float delta_t = 0.01;
            bool is_successful = false;
            bool is_inf_failed = false;
            float one_half_delta_t;   //> 1/2 \Delta t
            bool end_zone = false;
            magma_int_t locinfo;
            magma_int_t nrhs = 1;
            magmaFloatComplex h_param_homotopy[Num_Of_Params+1] = {MAGMA_C_ZERO};

            //> Pre-compute some offsets
            int hc_track_sol_offset = s*(Num_Of_Vars+1) + ri*RANSAC_Sol_Offset_with_Dummy_Variable;
            int cgsevA_offset       = s*(Num_Of_Vars*Num_Of_Vars) + ri*RANSAC_Sol_Offset*Num_Of_Vars;
            int cgesvB_offset       = s*(Num_Of_Vars) + ri*RANSAC_Sol_Offset;
            int ipiv_offset         = cgesvB_offset;
            int param_offset        = ri*(Num_Of_Params+1);
            
            //> Loop over all HC steps
            for(int step = 0; step <= CPUHC_Max_Steps; step++) {
                if (t0 < 1.0 && (1.0-t0 > 0.0000001)) {

                    // ===================================================================
                    // Decide delta t at end zone
                    // ===================================================================
                    if (!end_zone && fabs(1 - t0) <= (0.0500001)) end_zone = true;
                    if (end_zone) {
                        if (delta_t > fabs(1 - t0))
                            delta_t = fabs(1 - t0);
                    }
                    else if (delta_t > fabs(0.95 - t0)) {
                        delta_t = fabs(0.95 - t0);
                    }
                    //> ==================================================================

                    t_step = t0;
                    one_half_delta_t = 0.5 * delta_t;

                    // ===================================================================
                    // Prediction: 4-th order Runge-Kutta
                    // ===================================================================
                    //> (i) a pair of Jacobian evaluation + linear system solver
                    cpu_eval_compute_param_homotopy( t0, h_param_homotopy, h_Start_Params, h_Target_Params + param_offset );
                    Eval_dHdX_func( Num_Of_Vars, dHdx_Max_Terms, dHdx_Entry_Offset, dHdx_Max_Parts, h_dHdx_Index, h_CPU_HC_Track_Sols + hc_track_sol_offset, h_param_homotopy, h_cgesvA + cgsevA_offset );
                    Eval_dHdt_func( Num_Of_Vars, dHdt_Max_Terms, dHdt_Max_Parts, h_dHdt_Index, h_CPU_HC_Track_Sols + hc_track_sol_offset, h_param_homotopy, h_cgesvB + cgesvB_offset, h_diff_params + param_offset );
                    lapackf77_cgesv( &Num_Of_Vars, &nrhs, h_cgesvA + cgsevA_offset, &Num_Of_Vars, ipiv + ipiv_offset, h_cgesvB + cgesvB_offset, &Num_Of_Vars, &locinfo );
                    CPUHC_get_Runge_Kutta_x_for_k2( h_Intermediate_Sols + hc_track_sol_offset, h_CPU_HC_Track_Sols + hc_track_sol_offset, h_cgesvB + cgesvB_offset, one_half_delta_t, delta_t, t0 );

                    //> (ii) a pair of Jacobian evaluation + linear system solver
                    cpu_eval_compute_param_homotopy( t0, h_param_homotopy, h_Start_Params, h_Target_Params + param_offset );
                    Eval_dHdX_func( Num_Of_Vars, dHdx_Max_Terms, dHdx_Entry_Offset, dHdx_Max_Parts, h_dHdx_Index, h_CPU_HC_Track_Sols + hc_track_sol_offset, h_param_homotopy, h_cgesvA + cgsevA_offset );                    
                    Eval_dHdt_func( Num_Of_Vars, dHdt_Max_Terms, dHdt_Max_Parts, h_dHdt_Index, h_CPU_HC_Track_Sols + hc_track_sol_offset, h_param_homotopy, h_cgesvB + cgesvB_offset, h_diff_params + param_offset );
                    lapackf77_cgesv( &Num_Of_Vars, &nrhs, h_cgesvA + cgsevA_offset, &Num_Of_Vars, ipiv + ipiv_offset, h_cgesvB + cgesvB_offset, &Num_Of_Vars, &locinfo );
                    CPUHC_get_Runge_Kutta_x_for_k3( h_Intermediate_Sols + hc_track_sol_offset, h_CPU_HC_Track_Sols + hc_track_sol_offset, h_cgesvB + cgesvB_offset, h_Track_Last_Success + hc_track_sol_offset, one_half_delta_t, delta_t );

                    //> (iii) linear system solver (Jacobian evaluation is unnecessary)
                    cpu_eval_compute_param_homotopy( t0, h_param_homotopy, h_Start_Params, h_Target_Params + param_offset );
                    Eval_dHdX_func( Num_Of_Vars, dHdx_Max_Terms, dHdx_Entry_Offset, dHdx_Max_Parts, h_dHdx_Index, h_CPU_HC_Track_Sols + hc_track_sol_offset, h_param_homotopy, h_cgesvA + cgsevA_offset );                    
                    Eval_dHdt_func( Num_Of_Vars, dHdt_Max_Terms, dHdt_Max_Parts, h_dHdt_Index, h_CPU_HC_Track_Sols + hc_track_sol_offset, h_param_homotopy, h_cgesvB + cgesvB_offset, h_diff_params + param_offset );
                    lapackf77_cgesv( &Num_Of_Vars, &nrhs, h_cgesvA + cgsevA_offset, &Num_Of_Vars, ipiv + ipiv_offset, h_cgesvB + cgesvB_offset, &Num_Of_Vars, &locinfo );
                    CPUHC_get_Runge_Kutta_x_for_k4( h_Intermediate_Sols + hc_track_sol_offset, h_CPU_HC_Track_Sols + hc_track_sol_offset, h_cgesvB + cgesvB_offset, h_Track_Last_Success + hc_track_sol_offset, one_half_delta_t, t0, delta_t );

                    //> (iv) a pair of Jacobian evaluation + linear system solver
                    cpu_eval_compute_param_homotopy( t0, h_param_homotopy, h_Start_Params, h_Target_Params + param_offset );
                    Eval_dHdX_func( Num_Of_Vars, dHdx_Max_Terms, dHdx_Entry_Offset, dHdx_Max_Parts, h_dHdx_Index, h_CPU_HC_Track_Sols + hc_track_sol_offset, h_param_homotopy, h_cgesvA + cgsevA_offset );                    
                    Eval_dHdt_func( Num_Of_Vars, dHdt_Max_Terms, dHdt_Max_Parts, h_dHdt_Index, h_CPU_HC_Track_Sols + hc_track_sol_offset, h_param_homotopy, h_cgesvB + cgesvB_offset, h_diff_params + param_offset );
                    lapackf77_cgesv( &Num_Of_Vars, &nrhs, h_cgesvA + cgsevA_offset, &Num_Of_Vars, ipiv + ipiv_offset, h_cgesvB + cgesvB_offset, &Num_Of_Vars, &locinfo );

                    //> make prediction
                    CPUHC_make_Runge_Kutta_prediction( h_Intermediate_Sols + hc_track_sol_offset, h_CPU_HC_Track_Sols + hc_track_sol_offset, h_cgesvB + cgesvB_offset, delta_t );

                    // ===================================================================
                    // Correction: Newton's method
                    // ===================================================================
                    for(int c = 0; c < CPUHC_Max_Correction_Steps; c++) {
                        
                        //> evaluate dH/dX and H for solving the corrected step in the linear system
                        Eval_dHdX_func( Num_Of_Vars, dHdx_Max_Terms, dHdx_Entry_Offset, dHdx_Max_Parts, h_dHdx_Index, h_CPU_HC_Track_Sols + hc_track_sol_offset, h_param_homotopy, h_cgesvA + cgsevA_offset );
                        Eval_H_func( Num_Of_Vars, dHdt_Max_Terms, dHdt_Max_Parts, h_dHdt_Index, h_CPU_HC_Track_Sols + hc_track_sol_offset, h_param_homotopy, h_cgesvB + cgesvB_offset );
                        lapackf77_cgesv( &Num_Of_Vars, &nrhs, h_cgesvA + cgsevA_offset, &Num_Of_Vars, ipiv + ipiv_offset, h_cgesvB + cgesvB_offset, &Num_Of_Vars, &locinfo );
                        
                        //> make correction; see if the corrected solution is successful
                        CPUHC_make_correction( h_CPU_HC_Track_Sols + hc_track_sol_offset, h_cgesvB + cgesvB_offset, is_successful, is_inf_failed );

                        //> a successful prediction or an infinity failed solution stops the correction
                        if (is_successful) break;
                        if (is_inf_failed) break;
                    }

                    //> break the entire loop if the solution is already infinity failed
                    if (is_inf_failed) {
                        h_is_Track_Inf_Failed[s + ri*Num_Of_Tracks] = true;
                        break;
                    }

                    // ===================================================================
                    // Decide Homotopy Path Changes
                    // ===================================================================
                    if (!is_successful) {
                        pred_success_count = 0;
                        delta_t *= 0.5;
                        t0 = t_step;
                        //> should be the last successful tracked sols
                        for (int i = 0; i < Num_Of_Vars; i++) {
                            (h_CPU_HC_Track_Sols + hc_track_sol_offset)[i] = (h_Track_Last_Success + hc_track_sol_offset)[i];
                            (h_Intermediate_Sols + hc_track_sol_offset)[i] = (h_Track_Last_Success + hc_track_sol_offset)[i];
                        }
                    }
                    else {
                        for (int i = 0; i < Num_Of_Vars; i++) {
                            (h_Track_Last_Success + hc_track_sol_offset)[i] = (h_CPU_HC_Track_Sols + hc_track_sol_offset)[i];
                            (h_Intermediate_Sols + hc_track_sol_offset)[i]  = (h_CPU_HC_Track_Sols + hc_track_sol_offset)[i];
                        }
                        pred_success_count++;
                        if (pred_success_count >= CPUHC_delta_t_incremental_steps) {
                            pred_success_count = 0;
                            delta_t *= 2;
                        }
                    }
                }
                else {
                    break;
                }
            }   //> loop over HC steps
            h_is_Track_Converged[s + ri*Num_Of_Tracks] = (t0 >= 1.0 || (1.0-t0 <= 0.0000001)) ? true : false;
        // }   //> loop over all homotopy paths
    }   //> loop over RANSAC iterations
    
    //> return CPU-HC timing
    return magma_wtime() - CPU_HC_time;
}

void CPU_HC_Solver::CPUHC_get_Runge_Kutta_x_for_k2(
    magmaFloatComplex *h_Intermediate_Sols, magmaFloatComplex *h_CPU_HC_Track_Sols, magmaFloatComplex *h_cgesvB,
    float one_half_delta_t, float delta_t, float &t0 ) 
{
    //> prepare data for k2
    for (int i = 0; i < Num_Of_Vars; i++) {
        (h_Intermediate_Sols)[i] += (h_cgesvB)[i] * delta_t * 1.0/6.0;  //> x = x + (\Delta t) (k1/6)
        (h_cgesvB)[i]            *= one_half_delta_t;                   //> k1 * (\Delta t)/2
        (h_CPU_HC_Track_Sols)[i] += (h_cgesvB)[i];                      //> x = x + k1 * (\Delta t)/2
    }
    t0 += one_half_delta_t;                                             //> t = t + (\Delta t)/2 
}

void CPU_HC_Solver::CPUHC_get_Runge_Kutta_x_for_k3(
    magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_CPU_HC_Track_Sols, 
    magmaFloatComplex *h_cgesvB, magmaFloatComplex *h_Track_Last_Success,
    float one_half_delta_t, float delta_t ) 
{
    //> prepare data for k3
    for (int i = 0; i < Num_Of_Vars; i++) {
        (h_Sols_cpu)[i]             += (h_cgesvB)[i] * delta_t * 1.0/3.0;  //> x = x + (\Delta t) (k1/6 + k2/3)
        (h_CPU_HC_Track_Sols)[i]    =  (h_Track_Last_Success)[i];          //> copy the initial prior prediction solution
        (h_cgesvB)[i]               *= one_half_delta_t;                   //> k2 * (\Delta t)/2
        (h_CPU_HC_Track_Sols)[i]    += (h_cgesvB)[i];                      //> x = x + k2 * (\Delta t)/2
    }
}

void CPU_HC_Solver::CPUHC_get_Runge_Kutta_x_for_k4(
    magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_CPU_HC_Track_Sols, 
    magmaFloatComplex *h_cgesvB, magmaFloatComplex *h_Track_Last_Success,
    float one_half_delta_t, float &t0, float delta_t ) 
{
    //> prepare data for k4
    for (int i = 0; i < Num_Of_Vars; i++) {
        (h_Sols_cpu)[i]             += (h_cgesvB)[i] * delta_t * 1.0/3.0;    //> x = x + (\Delta t) (k1/6 + k2/3 + k3/3)
        (h_CPU_HC_Track_Sols)[i]    =  (h_Track_Last_Success)[i];            //> copy the initial prior prediction solution
        (h_cgesvB)[i]               *= delta_t;                              //> k3 * (\Delta t)
        (h_CPU_HC_Track_Sols)[i]    += (h_cgesvB)[i];                        //> x = x + k3 * (\Delta t)
    }
    t0 += one_half_delta_t;                                                                                                                 //> now t becomes t = t + (\Delta t)
}

void CPU_HC_Solver::CPUHC_make_Runge_Kutta_prediction(
    magmaFloatComplex *h_Sols_cpu, magmaFloatComplex *h_CPU_HC_Track_Sols, 
    magmaFloatComplex *h_cgesvB, float delta_t ) 
{
    for (int i = 0; i < Num_Of_Vars; i++) {
        (h_Sols_cpu)[i] += (h_cgesvB)[i] * delta_t * 1.0/6.0;
        (h_CPU_HC_Track_Sols)[i] =  (h_Sols_cpu)[i];
    }
}

void CPU_HC_Solver::CPUHC_make_correction( 
    magmaFloatComplex *h_CPU_HC_Track_Sols, magmaFloatComplex *h_cgesvB,
    bool &is_successful, bool &is_inf_failed ) 
{
    float sqrt_sols = 0.0;
    float sqrt_corr = 0.0;

    //> make correction and compute norm
    for (int i = 0; i < Num_Of_Vars; i++) {
        (h_CPU_HC_Track_Sols)[i] -= (h_cgesvB)[i];
        sqrt_sols += MAGMA_C_REAL((h_cgesvB)[i])*MAGMA_C_REAL((h_cgesvB)[i]) + MAGMA_C_IMAG((h_cgesvB)[i])*MAGMA_C_IMAG((h_cgesvB)[i]);
        sqrt_corr += MAGMA_C_REAL((h_CPU_HC_Track_Sols)[i])*MAGMA_C_REAL((h_CPU_HC_Track_Sols)[i]) + MAGMA_C_IMAG((h_CPU_HC_Track_Sols)[i])*MAGMA_C_IMAG((h_CPU_HC_Track_Sols)[i]);
    }

    //> check if the correction is successful
    is_successful = (sqrt_sols < 0.000001 * sqrt_corr);

    //> check if the correction is infinity-failed
    is_inf_failed = (sqrt_corr > 1e14);
}

#endif
