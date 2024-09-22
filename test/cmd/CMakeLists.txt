#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> main lib
find_package(Threads REQUIRED)
link_libraries(magmaHC)

add_executable(magmaHC-main magmaHC-main.cpp )

set_source_files_properties(magmaHC-main.cpp PROPERTIES LANGUAGE CUDA)

set_property(TARGET magmaHC-main
             PROPERTY CUDA_SEPARABLE_COMPILATION ON)

target_link_libraries(magmaHC-main Threads::Threads)

