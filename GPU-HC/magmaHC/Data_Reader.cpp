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

#include "Data_Reader.hpp"

Data_Reader::Data_Reader(std::string Problem_Filename, std::string RANSAC_Data_File_Path, const int Num_Of_Tracks, const int Num_Of_Vars, const int Num_Of_Params) \
    : num_of_tracks(Num_Of_Tracks), num_of_variables(Num_Of_Vars), num_of_params(Num_Of_Params) {

    //> Define problem file names
    File_Name_Start_Params = Problem_Filename + std::string("/start_params.txt");
    File_Name_Target_Params = Problem_Filename + std::string("/target_params.txt");
    File_Name_Start_Sols = Problem_Filename + std::string("/start_sols.txt");

    File_Name_dHdx_Indx = Problem_Filename + std::string("/dHdx_indx.txt");
    File_Name_dHdt_Indx = Problem_Filename + std::string("/dHdt_indx.txt");

    File_Name_Intrinsic_Matrix = RANSAC_Data_File_Path + std::string("/Synthetic_Intrinsic_Matrix.txt");
    File_Name_Pose21           = RANSAC_Data_File_Path + std::string("/Synthetic_GT_Poses21.txt");
    File_Name_Pose31           = RANSAC_Data_File_Path + std::string("/Synthetic_GT_Poses31.txt");
    File_Name_Triplet_Edgels   = RANSAC_Data_File_Path + std::string("/Synthetic_Triplet_Edgels.txt");
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
        
    }

    //> Copy the start solutions a number of RANSAC iterations times for h_Homotopy_Sols array
    for (int ri = 1; ri < NUM_OF_RANSAC_ITERATIONS; ri++) {
        for (int i = 0; i < (num_of_tracks * (num_of_variables + 1)); i++) {
            (h_Homotopy_Sols + ri * (num_of_tracks * (num_of_variables + 1)))[i] = h_Homotopy_Sols[i];
        }
    }
#if RANSAC_DEBUG
    int copy_id = 2;
    if (NUM_OF_RANSAC_ITERATIONS > copy_id) {
        std::cout << "Printing coopied h_Homotopy_Sols:" << std::endl;
        printf("     Copy 1                       Copy 2\n");
        for (int i = 0; i < num_of_variables+1; i++) {
            printf(" (%.7f, %.7f)        (%.7f, %.7f)\n", MAGMA_C_REAL(h_Homotopy_Sols[i]), MAGMA_C_IMAG(h_Homotopy_Sols[i]), \
                                                        MAGMA_C_REAL((h_Homotopy_Sols + copy_id*(num_of_tracks * (num_of_variables + 1)))[i]), \
                                                        MAGMA_C_IMAG((h_Homotopy_Sols + copy_id*(num_of_tracks * (num_of_variables + 1)))[i]));
        }
    }
#endif
    return true;
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
        h_Target_Params[num_of_params] = MAGMA_C_ONE;
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
        h_Start_Params[num_of_params] = MAGMA_C_ONE;
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

bool Data_Reader::Read_Camera_Matrices( float Pose21[12], float Pose31[12], float K[9] ) {
    //> Intrinsic matrix
    File_Intrinsic_Matrix.open(File_Name_Intrinsic_Matrix, std::ios_base::in);
    if (!File_Intrinsic_Matrix) {
        LOG_FILE_ERROR(File_Name_Intrinsic_Matrix);
        return false;
    }
    else {
        int d = 0;
        float entry = 0.0;
        while (File_Intrinsic_Matrix >> entry) {
            K[d] = entry;
            d++;
        }
    }

    //> Extrinsic matrix - Camera relative pose view 1 & 2
    File_Pose21.open(File_Name_Pose21, std::ios_base::in);
    if (!File_Pose21) {
        LOG_FILE_ERROR(File_Name_Pose21);
        return false;
    }
    else {
        int d = 0;
        float entry = 0.0;
        while (File_Pose21 >> entry) {
            Pose21[d] = entry;
            d++;
        }
    }

    //> Extrinsic matrix - Camera relative pose view 1 & 3
    File_Pose31.open(File_Name_Pose31, std::ios_base::in);
    if (!File_Pose31) {
        LOG_FILE_ERROR(File_Name_Pose31);
        return false;
    }
    else {
        int d = 0;
        float entry = 0.0;
        while (File_Pose31 >> entry) {
            Pose31[d] = entry;
            d++;
        }
    }
#if DATA_READER_DEBUG
    std::cout << std::endl << "Printing K ..." << std::endl;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            printf("%.6f\t", K[i*3 + j]);
        printf("\n");
    }
    std::cout << std::endl;
    std::cout << "Printing Pose21 ..." << std::endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 3; j++)
            printf("%.6f\t", Pose21[i*3 + j]);
        printf("\n");
    }
    std::cout << std::endl;
    std::cout << "Printing Pose31 ..." << std::endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 3; j++)
            printf("%.6f\t", Pose31[i*3 + j]);
        printf("\n");
    }
#endif
    return true;
}

bool Data_Reader::Read_Triplet_Edgels( float* &Triplet_Edge_Locations, float* &Triplet_Edge_Tangents ) {
    File_Triplet_Edgels.open(File_Name_Triplet_Edgels, std::ios_base::in);
    if (!File_Triplet_Edgels) {
        LOG_FILE_ERROR(File_Name_Triplet_Edgels);
        return false;
    }
    else {
        int edge_d = 0, tgt_d = 0;
        float edge1_x, edge1_y, tgt1_x, tgt1_y, edge2_x, edge2_y, tgt2_x, tgt2_y, edge3_x, edge3_y, tgt3_x, tgt3_y;
        while (File_Triplet_Edgels >> edge1_x >> edge1_y >> tgt1_x >> tgt1_y >> \
                                      edge2_x >> edge2_y >> tgt2_x >> tgt2_y >> \
                                      edge3_x >> edge3_y >> tgt3_x >> tgt3_y) {
            //> Edge locations
            Triplet_Edge_Locations(edge_d, 0) = edge1_x;
            Triplet_Edge_Locations(edge_d, 1) = edge1_y;
            Triplet_Edge_Locations(edge_d, 2) = edge2_x;
            Triplet_Edge_Locations(edge_d, 3) = edge2_y;
            Triplet_Edge_Locations(edge_d, 4) = edge3_x;
            Triplet_Edge_Locations(edge_d, 5) = edge3_y;
            edge_d++;

            //> edge tangents
            Triplet_Edge_Tangents(tgt_d, 0) = tgt1_x;
            Triplet_Edge_Tangents(tgt_d, 1) = tgt1_y;
            Triplet_Edge_Tangents(tgt_d, 2) = tgt2_x;
            Triplet_Edge_Tangents(tgt_d, 3) = tgt2_y;
            Triplet_Edge_Tangents(tgt_d, 4) = tgt3_x;
            Triplet_Edge_Tangents(tgt_d, 5) = tgt3_y;
            tgt_d++;
        }
#if DATA_READER_DEBUG
        std::cout << "Number of triplet edgels = " << edge_d << std::endl;
        std::cout << std::endl << "Printing Triplet_Edge_Locations ..." << std::endl;
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 6; j++)
                printf("%.6f\t", Triplet_Edge_Locations[i*6 + j]);
        std::cout << std::endl;
        std::cout << "Printing Triplet_Edge_Tangents ..." << std::endl;
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 6; j++)
                printf("%.6f\t", Triplet_Edge_Tangents[i*6 + j]);
#endif
        return true;
    }
}

void Data_Reader::Print_Out_Target_Params_from_Triplet_Edgels(int sample_index, std::vector<std::array<int,3>> target_params_match_indices, magmaFloatComplex *h_Target_Params) {
    
    std::array<int,3> print_Triplet_Index = target_params_match_indices[sample_index];
    std::cout << std::endl << "Printing triplet edgel indices: ";
    std::cout << print_Triplet_Index[0] << " " << print_Triplet_Index[1] << " " << print_Triplet_Index[2] << std::endl;
    std::cout << "Converting from triplet edgels to target parameters:" << std::endl;
    for (int i = 0; i <= num_of_params; i++) {
        printf("(%.10f, %.10f)\n", MAGMA_C_REAL((h_Target_Params + sample_index*(num_of_params+1))[i]), \
                                   MAGMA_C_IMAG((h_Target_Params + sample_index*(num_of_params+1))[i]));
    }
}

#if USE_8BIT_IN_SHARED_MEM == false
template bool Data_Reader::Read_dHdx_Indices< int >( int* & );
template bool Data_Reader::Read_dHdt_Indices< int >( int* & );
#endif

template bool Data_Reader::Read_dHdx_Indices< char >( char* & );
template bool Data_Reader::Read_dHdt_Indices< char >( char* & );

#endif
