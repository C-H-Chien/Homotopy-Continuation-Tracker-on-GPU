find_package(Threads REQUIRED)
link_libraries(magmaHC)
#link_libraries(magmaHC-DP64)

add_executable(magmaHC-main magmaHC-main.cu)

set_property(TARGET magmaHC-main
             PROPERTY CUDA_SEPARABLE_COMPILATION ON)

target_link_libraries(magmaHC-main Threads::Threads)

#> SP32
#add_executable(magmaHC-main-multiple-params magmaHC-main-multiple-params.cu)

#set_property(TARGET magmaHC-main-multiple-params
#             PROPERTY CUDA_SEPARABLE_COMPILATION ON)

#target_link_libraries(magmaHC-main-multiple-params Threads::Threads)

#> Used for reading multiple target parameters, DP64
#add_executable(magmaHC-main-multiple-params-dp64 magmaHC-main-multiple-params-dp64.cu)

#set_property(TARGET magmaHC-main-multiple-params-dp64
#             PROPERTY CUDA_SEPARABLE_COMPILATION ON)

#target_link_libraries(magmaHC-main-multiple-params-dp64 Threads::Threads)