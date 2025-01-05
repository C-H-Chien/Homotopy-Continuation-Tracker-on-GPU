cmake_minimum_required(VERSION 3.10.2 FATAL_ERROR)
message(STATUS "CMake version: ${CMAKE_VERSION}")

set(CMAKE_CXX_STANDARD "11")
if(NOT CMAKE_CXX_STANDARD_REQUIRED)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

set(MAGMADIR    "/home/ahmad/magma-openblas")
set(OPENBLASDIR "/home/ahmad/openblas")
set(YAMLCPP_DIR "/home/ahmad/yaml-cpp")

# this is for custom options on ICL machines
#set(ICL_BUILD 1)

set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math" CACHE STRING "" FORCE)

project(TrifocalRelativePose)

# -- enable the use of cuda --
option(USE_CUDA "Use CUDA" ON)
enable_language("CUDA")
include(${CMAKE_CURRENT_SOURCE_DIR}/config/cmake/InitializeBuildType.cmake)

# -- Crucial flags for performance --
if(DEFINED ICL_BUILD)
    # for building from the repo, explicitly add the DADD_ option
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -DADD_ -march=native" CACHE STRING "" FORCE)

    # build for volta arch. (compute capability 70)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_70,code=sm_70")

    # build for ampere arch. (compute capability 80)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_80,code=sm_80")

    # build for ampere arch. (compute capability 90)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_90,code=sm_90")

    # add the DADD_ option to nvcc
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -DADD_")

    # this option lets nvcc give some info about register usage/spills
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xptxas -v")

    # this option (currently commented out) can be used to limit the use of registers
    # set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --maxrregcount 128")

else()
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native" CACHE STRING "" FORCE)

    # build for volta arch. (compute capability 70)
    #set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_70,code=sm_70")

    # build for ampere arch. (compute capability 80)
    #set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_80,code=sm_80")

    # build for ampere arch. (compute capability 90)
    #set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_90,code=sm_90")

    # add the DADD_ option to nvcc
    #set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -DADD_")

    # this option lets nvcc give some info about register usage/spills
    #set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xptxas -v")

endif()



#-----------------------------------------------------------------------------
if(NOT COMMAND SETIFEMPTY)
  macro(SETIFEMPTY)
    set(KEY ${ARGV0})
    set(VALUE ${ARGV1})
    if(NOT ${KEY})
      set(${ARGV})
    endif()
  endmacro()
endif()

#-----------------------------------------------------------------------------
SETIFEMPTY(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SETIFEMPTY(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SETIFEMPTY(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

#-----------------------------------------------------------------------------
SETIFEMPTY(CMAKE_INSTALL_LIBRARY_DESTINATION lib)
SETIFEMPTY(CMAKE_INSTALL_ARCHIVE_DESTINATION lib)
SETIFEMPTY(CMAKE_INSTALL_RUNTIME_DESTINATION bin)

SETIFEMPTY(MINUS_INSTALL_RUNTIME_DIR bin)
SETIFEMPTY(MINUS_INSTALL_LIBRARY_DIR lib)
SETIFEMPTY(MINUS_INSTALL_ARCHIVE_DIR lib)
SETIFEMPTY(PROJ_INSTALL_INCLUDE_DIR include/)

# -- CMake support directory --
set(PROJ_ROOT_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR})
set(PROJ_CMAKE_DIR ${CMAKE_CURRENT_LIST_DIR}/config/cmake/Modules)
include(${PROJ_CMAKE_DIR}/StandardOptions.cmake)

# -- Options to add extra compiler and linker flags --
#
# -- These options allow you to specify additional flags without --
# -- affecting the default flags for a particular platform or build type --
# -- This is especially useful for adding extra warning flags --
set( MINUS_EXTRA_CMAKE_C_FLAGS             CACHE STRING "Extra flags appended to CMAKE_C_FLAGS" )
set( MINUS_EXTRA_CMAKE_CXX_FLAGS           CACHE STRING "Extra flags appended to CMAKE_CXX_FLAGS" )
set( MINUS_EXTRA_CMAKE_EXE_LINKER_FLAGS    CACHE STRING "Extra flags appended to CMAKE_EXE_LINKER_FLAGS" )
set( MINUS_EXTRA_CMAKE_MODULE_LINKER_FLAGS CACHE STRING "Extra flags appended to CMAKE_MODULE_LINKER_FLAGS" )
set( MINUS_EXTRA_CMAKE_SHARED_LINKER_FLAGS CACHE STRING "Extra flags appended to CMAKE_SHARED_LINKER_FLAGS" )
set( CMAKE_C_FLAGS             "${CMAKE_C_FLAGS} ${MINUS_EXTRA_CMAKE_C_FLAGS}" )
set( CMAKE_CXX_FLAGS           "${CMAKE_CXX_FLAGS} ${MINUS_EXTRA_CMAKE_CXX_FLAGS}" )
set( CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS} ${MINUS_EXTRA_CMAKE_EXE_LINKER_FLAGS}" )
set( CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${MINUS_EXTRA_CMAKE_MODULE_LINKER_FLAGS}" )
set( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${MINUS_EXTRA_CMAKE_SHARED_LINKER_FLAGS}" )

# -- import openmp to the compiler --
find_package(OpenMP REQUIRED)
#find_package(OpenMP)
#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

#-------------------------------------------------------------------
#-- BUILD CONFIG OPTIONS AND WARNINGS --

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # require at least gcc 4.8
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0)
      message(WARNING "***********************************")
      message(WARNING "GCC must be at least 5.0 for speed!")
      message(WARNING "***********************************")
    endif()
elseif (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
else()
    message(WARNING "You are using an unsupported compiler! Compilation has only been tested with Clang and GCC.")
endif()


# -- Avoid types of path names that could cause havoc with regular expressions --
if( ${PROJECT_BINARY_DIR} MATCHES  \\+ )
  message(SEND_ERROR "You cannot have a + in your binary path")
endif()
if( ${CMAKE_CURRENT_LIST_DIR} MATCHES  \\+ )
  message(SEND_ERROR "You cannot have a + in your source path")
endif()

# include CMakeListsLocal.txt from source directory if it exists
# also include it from the binary dir if different from source dir
if( ${PROJECT_BINARY_DIR} MATCHES ${CMAKE_CURRENT_LIST_DIR} )
  include( ${CMAKE_CURRENT_LIST_DIR}/CMakeListsLocal.txt OPTIONAL )
else()
  include( ${CMAKE_CURRENT_LIST_DIR}/CMakeListsLocal.txt OPTIONAL )
  include( ${PROJECT_BINARY_DIR}/CMakeListsLocal.txt OPTIONAL )
endif()

set(PROJ_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR})    #Source includes
set(PROJ_BINARY_INCLUDE_DIR ${PROJECT_BINARY_DIR}) #Generated includes
set(PROJ_INCLUDE_DIRS ${PROJ_BINARY_INCLUDE_DIR} ${PROJ_INCLUDE_DIR})
set(PROJ_INSTALL_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include/)

# -- include necessary libraries in the commands / executables file --
include_directories(${PROJ_INCLUDE_DIRS})

if(DEFINED ICL_BUILD)
  include_directories( ${MAGMADIR}/include )
  include_directories( ${MAGMADIR}/testing )
  include_directories( ${MAGMADIR}/control )
  include_directories( ${MAGMADIR}/magmablas )
  include_directories( ${YAMLCPP_DIR}/include )
  link_directories(${MAGMADIR}/lib)
  link_directories(${OPENBLASDIR})
  link_directories(${YAMLCPP_DIR}/build)

else()
  include_directories( "/users/cchien3/data/cchien3/magma-2.5.4/include" )
  include_directories( "/users/cchien3/data/cchien3/magma-2.5.4/control" )
  include_directories( "/users/cchien3/data/cchien3/magma-2.5.4/magmablas" )
  include_directories( "/gpfs/data/bkimia/cchien3/yaml-cpp/bin/include" )
  link_directories(/gpfs/data/bkimia/cchien3/yaml-cpp/bin/lib64)
endif()

if(DEFINED ICL_BUILD)
    link_directories("/apps/spacks/2023-11-15/opt/spack/linux-rocky9-x86_64/gcc-11.3.1/gcc-7.5.0-aembyabek5wcmaengfp7yn6c22xe6bq7/lib64/")
endif()

# adds a test for Endianness and a global variable that should be useful for file format programming -MM
include(TestBigEndian)
TEST_BIG_ENDIAN(MY_BIG_ENDIAN)

# main lib
add_subdirectory( magmaHC )
# commands / executables
add_subdirectory(cmd)