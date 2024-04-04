#ifndef DATA_READER_CPP
#define DATA_READER_CPP
// ============================================================================
// Data_Reader class CPP: read data from the problem files
//
// Changelogs
//    Chien  24-01-21:   Initially Created.
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
#include <cuComplex.h>

#include "magma_v2.h"

#include "definitions.hpp"
#include "Data_Reader.hpp"

Data_Reader::Data_Reader(std::string Problem_Filename, const int Num_Of_Tracks, const int Num_Of_Vars) \
    : num_of_tracks(Num_Of_Tracks), num_of_variables(Num_Of_Vars) {

    //> Define problem file names
    File_Name_Start_Params = Problem_Filename + std::string("/start_params.txt");
    File_Name_Target_Params = Problem_Filename + std::string("/target_params.txt");
    File_Name_Start_Sols = Problem_Filename + std::string("/start_sols.txt");

    File_Name_dHdx_Indx = Problem_Filename + std::string("/dHdx_indx.txt");
    File_Name_dHdt_Indx = Problem_Filename + std::string("/dHdt_indx.txt");
}

bool Data_Reader::Read_Start_Sols(magmaFloatComplex* &h_Start_Sols, magmaFloatComplex* &h_Homotopy_Sols) {
    //> Read start solutions
    File_Start_Sols.open(File_Name_Start_Sols, std::ios_base::in);
    if (!File_Start_Sols) {
        LOG_FILE_ERROR(File_Name_Start_Sols);
        return false;
    }
    else {
        float s_real, s_imag;
        int d = 0, i = 0; 
        while (File_Start_Sols >> s_real >> s_imag) {
            (h_Start_Sols + i * (num_of_variables+1))[d]     = MAGMA_C_MAKE(s_real, s_imag);
            (h_Homotopy_Sols + i * (num_of_variables+1))[d]  = MAGMA_C_MAKE(s_real, s_imag);
            if (d < num_of_variables-1) d++;
            else {
                d = 0;
                i++;
            }
        }
        for(int k = 0; k < num_of_tracks; k++) {
            (h_Start_Sols + k * (num_of_variables+1))[num_of_variables]    = MAGMA_C_MAKE(1.0, 0.0);
            (h_Homotopy_Sols + k * (num_of_variables+1))[num_of_variables] = MAGMA_C_MAKE(1.0, 0.0);
        }
        return true;
    }
}

bool Data_Reader::Read_Target_Params(magmaFloatComplex* &h_Target_Params) {
    int d = 0;
    File_Target_Params.open(File_Name_Target_Params, std::ios_base::in);
    if (!File_Target_Params) {
        LOG_FILE_ERROR(File_Name_Target_Params);
        return false;
    }
    else {
        float s_real, s_imag;
        while (File_Target_Params >> s_real >> s_imag) {
            h_Target_Params[d] = MAGMA_C_MAKE(s_real, s_imag);
            d++;
        }
        return true;
    }
}

bool Data_Reader::Read_Start_Params(magmaFloatComplex* &h_Start_Params) {
    int d = 0;
    File_Start_Params.open(File_Name_Start_Params, std::ios_base::in);
    if (!File_Start_Params) {
        LOG_FILE_ERROR(File_Name_Start_Params);
        return false;
    }
    else {
        float s_real, s_imag;
        while (File_Start_Params >> s_real >> s_imag) {
            h_Start_Params[d] = MAGMA_C_MAKE(s_real, s_imag);
            d++;
        }
        return true;
    }
}

template< typename T >
bool Data_Reader::Read_dHdx_Indices( T* &h_dHdx_Index ) {
    int index, d = 0;
    File_dHdx_Indices.open(File_Name_dHdx_Indx, std::ios_base::in);
    if (!File_dHdx_Indices) {
        LOG_FILE_ERROR(File_Name_dHdx_Indx);
        return false;
    }
    else {
        while (File_dHdx_Indices >> index) {
            (h_dHdx_Index)[d] = (T)index;
            d++;
        }
#if DATA_READER_DEBUG
    std::cout << "Printing h_dHdx_Index ..." << std::endl;
    for (int i = 0; i < 10; i++) printf("%d\t", (int)h_dHdx_Index[i]);
    std::cout << std::endl;
#endif
        return true;
    }
}

template< typename T >
bool Data_Reader::Read_dHdt_Indices( T* &h_dHdt_Index ) {
    int index, d = 0;
    File_dHdt_Indices.open(File_Name_dHdt_Indx, std::ios_base::in);
    if (!File_dHdt_Indices) {
        LOG_FILE_ERROR(File_Name_dHdt_Indx);
        return false;
    }
    else {
        while (File_dHdt_Indices >> index) {
            (h_dHdt_Index)[d] = (T)index;
            d++;
        }
#if DATA_READER_DEBUG
    std::cout << "Printing h_dHdt_Index ..." << std::endl;
    for (int i = 0; i < 10; i++) printf("%d\t", (int)h_dHdt_Index[i]);
    std::cout << std::endl;
#endif
        return true;
    }
}

#if USE_8BIT_IN_SHARED_MEM == false
template bool Data_Reader::Read_dHdx_Indices< int >( int* & );
template bool Data_Reader::Read_dHdt_Indices< int >( int* & );
#endif

template bool Data_Reader::Read_dHdx_Indices< char >( char* & );
template bool Data_Reader::Read_dHdt_Indices< char >( char* & );

#endif
