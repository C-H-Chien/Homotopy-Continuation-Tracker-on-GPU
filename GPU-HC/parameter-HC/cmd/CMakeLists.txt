find_package(Threads REQUIRED)
link_libraries(magmaHC)

add_executable(magmaHC-main magmaHC-main.cu)

set_property(TARGET magmaHC-main
             PROPERTY CUDA_SEPARABLE_COMPILATION ON)

target_link_libraries(magmaHC-main Threads::Threads)
