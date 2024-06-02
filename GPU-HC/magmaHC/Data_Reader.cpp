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
#include "./PHC_Coeffs/p2c-5pt_rel_pos_alg_form_quat.h"
#include "./PHC_Coeffs/p2c-5pt_rel_pos_geo_form_quat.h"
#include "./PHC_Coeffs/p2c-5pt_rel_pos_full_form_quat.h"
#include "./PHC_Coeffs/p2c-trifocal_2op1p_30x30.h"
#include "./PHC_Coeffs/p2c-generalized_3views_4pts.h"
#include "./PHC_Coeffs/p2c-numerical_generalized_3views_6lines.h"
#include "./PHC_Coeffs/p2c-uncalibrated_trifocal_rel_pos_CY.h"
#include "./PHC_Coeffs/p2c-uncalibrated_trifocal_rel_pos_CH.h"
#include "./PHC_Coeffs/p2c-optimal_PnP_quat.h"
#include "./PHC_Coeffs/p2c-3view_triangulation.h"
#include "./PHC_Coeffs/p2c-4view_triangulation.h"
#include "./PHC_Coeffs/p2c-6pt_RS_abs_pos.h"
#include "./PHC_Coeffs/p2c-6pt_RS_abs_pos_1lin.h"
#include "./PHC_Coeffs/p2c-dual_reciever_TDOA_5pt.h"
#include "./PHC_Coeffs/p2c-distorted_2view_triangulation.h"
#include "./PHC_Coeffs/p2c-optimal_P4P_abs_pos.h"
#include "./PHC_Coeffs/p2c-3pt_rel_pos_homog.h"
#include "./PHC_Coeffs/p2c-PnP_unkn_principal_pt.h"
#include "./PHC_Coeffs/p2c-rel_pos_quiver.h"
#include "./PHC_Coeffs/p2c-P3P.h"

Data_Reader::Data_Reader(std::string Problem_Filename, std::string RANSAC_Data_File_Path, const int Num_Of_Tracks, const int Num_Of_Vars, const int Num_Of_Params) \
    : num_of_tracks(Num_Of_Tracks), num_of_variables(Num_Of_Vars), num_of_params(Num_Of_Params), RANSAC_Data_Path_(RANSAC_Data_File_Path) {

    //> Define problem file names
    File_Name_Start_Params = Problem_Filename + std::string("/start_params.txt");
    File_Name_Target_Params = Problem_Filename + std::string("/target_params.txt");
    File_Name_Start_Sols = Problem_Filename + std::string("/start_sols.txt");

    File_Name_dHdx_Indx = Problem_Filename + std::string("/dHdx_indx.txt");
    File_Name_dHdt_Indx = Problem_Filename + std::string("/dHdt_indx.txt");

    File_Name_Intrinsic_Matrix = RANSAC_Data_Path_ + std::string("/Intrinsic_Matrix.txt");
}

bool Data_Reader::Construct_Coeffs_From_Params( std::string HC_Problem, \
        magmaFloatComplex* h_Target_Params,     magmaFloatComplex* h_Start_Params, \
        magmaFloatComplex* &h_dHdx_PHC_Coeffs,  magmaFloatComplex* &h_dHdt_PHC_Coeffs ) 
{
    if (HC_Problem == "5pt_rel_pos_geo_form_quat")              magmaHCWrapper::p2c_5pt_rel_pos_geo_form_quat(h_Target_Params, h_Start_Params, h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs);
    else if (HC_Problem == "5pt_rel_pos_alg_form_quat")         magmaHCWrapper::p2c_5pt_rel_pos_alg_form_quat(h_Target_Params, h_Start_Params, h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs);
    else if (HC_Problem == "5pt_rel_pos_full_form_quat")        magmaHCWrapper::p2c_5pt_rel_pos_full_form_quat(h_Target_Params, h_Start_Params, h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs);
    else if (HC_Problem == "trifocal_2op1p_30x30")              magmaHCWrapper::p2c_trifocal_2op1p_30x30(h_Target_Params, h_Start_Params, h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs);
    else if (HC_Problem == "generalized_3views_4pts")           magmaHCWrapper::p2c_generalized_3views_4pts(h_Target_Params, h_Start_Params, h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs);
    else if (HC_Problem == "generalized_3views_6lines")         magmaHCWrapper::p2c_numerical_generalized_3views_6lines(h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs);
    else if (HC_Problem == "uncalibrated_trifocal_rel_pos_CY")  magmaHCWrapper::p2c_uncalibrated_trifocal_rel_pos_CY(h_Target_Params, h_Start_Params, h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs);
    else if (HC_Problem == "uncalibrated_trifocal_rel_pos_CH")  magmaHCWrapper::p2c_uncalibrated_trifocal_rel_pos_CH(h_Target_Params, h_Start_Params, h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs);
    else if (HC_Problem == "optimal_PnP_quat")                  magmaHCWrapper::p2c_optimal_PnP_quat(h_Target_Params, h_Start_Params, h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs);
    else if (HC_Problem == "3view_triangulation")               magmaHCWrapper::p2c_3view_triangulation(h_Target_Params, h_Start_Params, h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs);
    else if (HC_Problem == "4view_triangulation")               magmaHCWrapper::p2c_4view_triangulation(h_Target_Params, h_Start_Params, h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs);
    else if (HC_Problem == "6pt_RS_abs_pos")                    magmaHCWrapper::p2c_6pt_RS_abs_pos(h_Target_Params, h_Start_Params, h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs);
    else if (HC_Problem == "6pt_RS_abs_pos_1lin")               magmaHCWrapper::p2c_6pt_RS_abs_pos_1lin(h_Target_Params, h_Start_Params, h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs);
    else if (HC_Problem == "dual_reciever_TDOA_5pt")            magmaHCWrapper::p2c_dual_reciever_TDOA_5pt(h_Target_Params, h_Start_Params, h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs);
    else if (HC_Problem == "distorted_2view_triangulation")     magmaHCWrapper::p2c_distorted_2view_triangulation(h_Target_Params, h_Start_Params, h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs);
    else if (HC_Problem == "optimal_P4P_abs_pos")               magmaHCWrapper::p2c_optimal_P4P_abs_pos(h_Target_Params, h_Start_Params, h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs);
    else if (HC_Problem == "3pt_rel_pos_homog")                 magmaHCWrapper::p2c_3pt_rel_pos_homog(h_Target_Params, h_Start_Params, h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs);
    else if (HC_Problem == "PnP_unkn_principal_pt")             magmaHCWrapper::p2c_PnP_unkn_principal_pt(h_Target_Params, h_Start_Params, h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs);
    else if (HC_Problem == "rel_pos_quiver")                    magmaHCWrapper::p2c_rel_pos_quiver(h_Target_Params, h_Start_Params, h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs);
    else if (HC_Problem == "P3P")                               magmaHCWrapper::p2c_P3P(h_Target_Params, h_Start_Params, h_dHdx_PHC_Coeffs, h_dHdt_PHC_Coeffs);
    else {
        LOG_ERROR("Invalid HC problem name or P2C function is not included in the Construct_Coeffs_From_Params.");
        return false;
    }
    return true;
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
    }

    // std::cout << "Target parameters:" << std::endl;
    // for (int i = 0; i <= num_of_params; i++)
    //     std::cout << MAGMA_C_REAL(h_Target_Params[i]) << "\t" << MAGMA_C_IMAG(h_Target_Params[i]) << std::endl;
    // std::cout << std::endl;
    return true;
}

bool Data_Reader::Read_Start_Params(magmaFloatComplex* &h_Start_Params) {
    int d = 0;
    File_Start_Params.open(File_Name_Start_Params, std::ios_base::in);
    // LOG_INFOR_MESG("Start params file name: " + File_Name_Start_Params);
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
    }

    // std::cout << "Start parameters:" << std::endl;
    // for (int i = 0; i <= num_of_params; i++)
    //     std::cout << MAGMA_C_REAL(h_Start_Params[i]) << "\t" << MAGMA_C_IMAG(h_Start_Params[i]) << std::endl;
    // std::cout << std::endl;
    return true;
}

bool Data_Reader::Read_dHdx_Indices( int* &h_dHdx_Index ) {
    int index, d = 0;
    File_dHdx_Indices.open(File_Name_dHdx_Indx, std::ios_base::in);
    if (!File_dHdx_Indices) {
        LOG_FILE_ERROR(File_Name_dHdx_Indx);
        return false;
    }
    else {
        while (File_dHdx_Indices >> index) {
            (h_dHdx_Index)[d] = (int)index;
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

bool Data_Reader::Read_dHdt_Indices( int* &h_dHdt_Index ) {
    int index, d = 0;
    File_dHdt_Indices.open(File_Name_dHdt_Indx, std::ios_base::in);
    if (!File_dHdt_Indices) {
        LOG_FILE_ERROR(File_Name_dHdt_Indx);
        return false;
    }
    else {
        while (File_dHdt_Indices >> index) {
            (h_dHdt_Index)[d] = (int)index;
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

bool Data_Reader::Read_Camera_Matrices( float Pose21[12], float Pose31[12], float K[9], int tp_index ) {

    //> Create padded file index
    std::string str_File_Index = std::to_string(tp_index);
    int min_str_length = (3 < str_File_Index.length()) ? 3 : str_File_Index.length();
    auto padded_Index = std::string(3 - min_str_length, '0') + str_File_Index;
    File_Name_Pose21  = RANSAC_Data_Path_ + "/GT_Poses21/GT_Poses21_" + padded_Index + ".txt";
    File_Name_Pose31  = RANSAC_Data_Path_ + "/GT_Poses31/GT_Poses31_" + padded_Index + ".txt";

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

//> Read triplet edgel correspondences file. Return number of triplet edgels.
int Data_Reader::get_Num_Of_Triplet_Edgels( int tp_index ) {

    //> Create padded file index
    std::string str_File_Index = std::to_string(tp_index);
    int min_str_length = (3 < str_File_Index.length()) ? 3 : str_File_Index.length();
    auto padded_Index = std::string(3 - min_str_length, '0') + str_File_Index;

    File_Name_Triplet_Edgels = RANSAC_Data_Path_ + "/Triplet_Edgels/Triplet_Edgels_" + padded_Index + ".txt";
    File_Triplet_Edgels.open(File_Name_Triplet_Edgels, std::ios_base::in);
    if (!File_Triplet_Edgels) {
        LOG_FILE_ERROR(File_Name_Triplet_Edgels);
        return 0;
    }
    else {
        float edge1_x, edge1_y, tgt1_x, tgt1_y, edge2_x, edge2_y, tgt2_x, tgt2_y, edge3_x, edge3_y, tgt3_x, tgt3_y;
        while (File_Triplet_Edgels >> edge1_x >> edge1_y >> tgt1_x >> tgt1_y >> \
                                      edge2_x >> edge2_y >> tgt2_x >> tgt2_y >> \
                                      edge3_x >> edge3_y >> tgt3_x >> tgt3_y) {

            auto edge1 = std::make_pair(edge1_x, edge1_y);
            auto edge2 = std::make_pair(edge2_x, edge2_y);
            auto edge3 = std::make_pair(edge3_x, edge3_y);
            auto tangent1 = std::make_pair(tgt1_x, tgt1_y);
            auto tangent2 = std::make_pair(tgt2_x, tgt2_y);
            auto tangent3 = std::make_pair(tgt3_x, tgt3_y);

            data_read_triplet_edgles_locations.push_back( std::make_tuple(edge1, edge2, edge3) );
            data_read_triplet_edgles_tangents.push_back( std::make_tuple(tangent1, tangent2, tangent3) );
        }
    }

    return data_read_triplet_edgles_locations.size();
}

void Data_Reader::Read_Triplet_Edgels( float* &Triplet_Edge_Locations, float* &Triplet_Edge_Tangents ) {

    //> Assign tuple values to arrays
    for (int ei = 0; ei < data_read_triplet_edgles_locations.size(); ei++) {
        Triplet_Edge_Locations(ei, 0) = std::get<0>(data_read_triplet_edgles_locations[ei]).first;
        Triplet_Edge_Locations(ei, 1) = std::get<0>(data_read_triplet_edgles_locations[ei]).second;
        Triplet_Edge_Locations(ei, 2) = std::get<1>(data_read_triplet_edgles_locations[ei]).first;
        Triplet_Edge_Locations(ei, 3) = std::get<1>(data_read_triplet_edgles_locations[ei]).second;
        Triplet_Edge_Locations(ei, 4) = std::get<2>(data_read_triplet_edgles_locations[ei]).first;
        Triplet_Edge_Locations(ei, 5) = std::get<2>(data_read_triplet_edgles_locations[ei]).second;

        Triplet_Edge_Tangents(ei, 0) = std::get<0>(data_read_triplet_edgles_tangents[ei]).first;
        Triplet_Edge_Tangents(ei, 1) = std::get<0>(data_read_triplet_edgles_tangents[ei]).second;
        Triplet_Edge_Tangents(ei, 2) = std::get<1>(data_read_triplet_edgles_tangents[ei]).first;
        Triplet_Edge_Tangents(ei, 3) = std::get<1>(data_read_triplet_edgles_tangents[ei]).second;
        Triplet_Edge_Tangents(ei, 4) = std::get<2>(data_read_triplet_edgles_tangents[ei]).first;
        Triplet_Edge_Tangents(ei, 5) = std::get<2>(data_read_triplet_edgles_tangents[ei]).second;
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

#endif
