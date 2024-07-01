
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
  gpu-kernels/kernel_GPUHC_trifocal_2op1p_30x30_GM32b.cu
  gpu-kernels/kernel_GPUHC_trifocal_2op1p_30x30_GM32b_RKL.cu
  gpu-kernels/kernel_GPUHC_trifocal_2op1p_30x30_GM32b_RKL_inline.cu
  gpu-kernels/kernel_GPUHC_trifocal_2op1p_30x30_GM32b_RKL_inline_LimUnroll.cu
  gpu-kernels/kernel_GPUHC_trifocal_2op1p_30x30_GM32b_RKL_LimUnroll.cu
  gpu-kernels/kernel_GPUHC_trifocal_2op1p_30x30_GM32b_RKL_LimUnroll_TrunPaths.cu
  gpu-kernels/kernel_GPUHC_trifocal_2op1p_30x30_GM32b_RKL_inline_LimUnroll_TrunPaths.cu
  gpu-kernels/kernel_GPUHC_trifocal_2op1p_30x30_GM8b_RKL_inline.cu
  gpu-kernels/kernel_GPUHC_trifocal_2op1p_30x30_SM8b_RKL.cu
  gpu-kernels/kernel_GPUHC_trifocal_2op1p_30x30_SM8b_RKL_inline.cu
  gpu-kernels/kernel_GPUHC_trifocal_2op1p_30x30_SM8b_RKL_LimUnroll.cu

  #>  device function
  dev-cgesv-batched-small.cuh
  dev-cgesv-batched-small-sm.cuh
  dev-get-new-data.cuh
  gpu-idx-evals/dev-eval-indxing-trifocal_2op1p_30x30_8b_GM_inline.cuh
  gpu-idx-evals/dev-eval-indxing-trifocal_2op1p_30x30_8b_SM.cuh
  gpu-idx-evals/dev-eval-indxing-trifocal_2op1p_30x30_8b_SM_inline.cuh
  gpu-idx-evals/dev-eval-indxing-trifocal_2op1p_30x30_8b_SM_LimUnroll.cuh
  gpu-idx-evals/dev-eval-indxing-trifocal_2op1p_30x30_32b.cuh
  gpu-idx-evals/dev-eval-indxing-trifocal_2op1p_30x30_32b_P2C.cuh
  gpu-idx-evals/dev-eval-indxing-trifocal_2op1p_30x30_32b_inline.cuh
  gpu-idx-evals/dev-eval-indxing-trifocal_2op1p_30x30_32b_LimUnroll.cuh
  gpu-idx-evals/dev-eval-indxing-trifocal_2op1p_30x30_32b_inline_LimUnroll.cuh 

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
  #-L/home/ahmad/magma-openblas/lib magma
  #-L/home/ahmad/openblas/lib openblas
        -L/apps/spacks/2023-08-14/opt/spack/linux-ubuntu23.04-x86_64/gcc-9.5.0/magma-2.7.1-5hrvjybyve3kqfgx7uv5nxxz7mktnqqw/lib magma
        -L/apps/spacks/2023-08-14/opt/spack/linux-ubuntu23.04-x86_64/gcc-9.5.0/openblas-0.3.23-2tips3wywhkieahh6udqmjfke4flqm3l/lib openblas
        -L/home/chchien/thirdParty/yaml-cpp/yaml-cpp-bin/lib64 yaml-cpp
)
#set_target_properties(magmaHC PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

