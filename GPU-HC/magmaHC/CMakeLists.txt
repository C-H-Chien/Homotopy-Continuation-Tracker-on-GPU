
set (control_sources

  GPU_HC_Solver.hpp
  GPU_HC_Solver.cpp

  #> Read problem data
  Data_Reader.hpp
  Data_Reader.cpp

  #> Evaluations on HC solutions
  Evaluations.hpp
  Evaluations.cpp

  #>  magma kernels
  gpu-kernels/magmaHC-kernels.hpp
  gpu-kernels/kernel_GPUHC_trifocal_2op1p_30x30_P2C.cu
  gpu-kernels/kernel_GPUHC_trifocal_2op1p_30x30_PH.cu
  gpu-kernels/kernel_GPUHC_trifocal_2op1p_30x30_PH_CodeOpt.cu
  gpu-kernels/kernel_GPUHC_trifocal_2op1p_30x30_PH_CodeOpt_Volta.cu
  gpu-kernels/kernel_GPUHC_trifocal_2op1p_30x30_PH_CodeOpt_TrunPaths.cu
  gpu-kernels/kernel_GPUHC_trifocal_2op1p_30x30_PH_CodeOpt_TrunPaths_Volta.cu
  gpu-kernels/kernel_GPUHC_trifocal_2op1p_30x30_PH_CodeOpt_TrunPaths_TrunRANSAC.cu
  gpu-kernels/kernel_GPUHC_trifocal_2op1p_30x30_PH_CodeOpt_TrunPaths_TrunRANSAC_Volta.cu

  #>  device function
  dev-cgesv-batched-small.cuh
  dev-cgesv-batched-small-sm.cuh
  dev-get-new-data.cuh
  dev-trifocal_2op1p-eval.cuh
  gpu-idx-evals/dev-eval-indxing-trifocal_2op1p_30x30.cuh
  gpu-idx-evals/dev-eval-indxing-trifocal_2op1p_30x30_P2C.cuh
  gpu-idx-evals/dev-eval-indxing-trifocal_2op1p_30x30_inline_LimUnroll.cuh
  gpu-idx-evals/dev-eval-indxing-trifocal_2op1p_30x30_LimUnroll_L2Cache.cuh
  
  #> P2C
  PHC_Coeffs/p2c-trifocal_2op1p_30x30.h

  #> Utilities
  util.hpp
)

set_source_files_properties(Data_Reader.cpp PROPERTIES LANGUAGE CUDA)
#set_source_files_properties(Data_Reader.hpp PROPERTIES LANGUAGE CUDA)
set_source_files_properties(GPU_HC_Solver.cpp PROPERTIES LANGUAGE CUDA)
#set_source_files_properties(GPU_HC_Solver.hpp PROPERTIES LANGUAGE CUDA)
set_source_files_properties(Evaluations.cpp PROPERTIES LANGUAGE CUDA)
#set_source_files_properties(Evaluations.hpp PROPERTIES LANGUAGE CUDA)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_library( magmaHC ${control_sources} )
target_link_libraries(magmaHC
  -L/apps/spacks/2023-08-14/opt/spack/linux-ubuntu23.04-x86_64/gcc-9.5.0/magma-2.7.1-5hrvjybyve3kqfgx7uv5nxxz7mktnqqw/lib magma
  -L/apps/spacks/2023-08-14/opt/spack/linux-ubuntu23.04-x86_64/gcc-9.5.0/openblas-0.3.23-2tips3wywhkieahh6udqmjfke4flqm3l/lib openblas
  -L/home/chchien/thirdParty/yaml-cpp/yaml-cpp-bin/lib64 yaml-cpp
)
#set_target_properties(magmaHC PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

