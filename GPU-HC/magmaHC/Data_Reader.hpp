#ifndef DATA_READER_H
#define DATA_READER_H
// ============================================================================
// Data_Reader class: read data from the problem files
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
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <chrono>
#include <cuComplex.h>

#include "definitions.hpp"

class Data_Reader {

public:
    //> Constructor
    Data_Reader(std::string);

    bool Read_Start_Params( cuFloatComplex* &h_Start_Params );
    bool Read_Target_Params( cuFloatComplex* &h_Target_Params );
    bool Read_Start_Sols( cuFloatComplex* &h_Start_Sols, cuFloatComplex* &h_Homotopy_Sols );

    bool Read_dHdx_Indices( int* &h_dHdx_Index );
    bool Read_dHdt_Indices( int* &h_dHdt_Index );

private:
    //> File names
    std::string File_Name_Target_Params;
    std::string File_Name_Start_Params;
    std::string File_Name_Start_Sols;
    std::string File_Name_dHdx_Indx;
    std::string File_Name_dHdt_Indx;

    //> input streams from problem files
    std::fstream File_Start_Params;
    std::fstream File_Target_Params;
    std::fstream File_Start_Sols;
    std::fstream File_dHdx_Indices;
    std::fstream File_dHdt_Indices;
};

#endif
