
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
  gpu-kernels/kernel_HC_Solver_5pt_rel_pos_geo_form_quat.cu
  gpu-kernels/kernel_HC_Solver_5pt_rel_pos_alg_form_quat.cu
  gpu-kernels/kernel_HC_Solver_trifocal_2op1p_30x30.cu

  #>  device function
  dev-cgesv-batched-small.cuh
  dev-get-new-data.cuh
  gpu-idx-evals/dev-eval-indxing-5pt_rel_pos_geo_form_quat.cuh
  gpu-idx-evals/dev-eval-indxing-5pt_rel_pos_alg_form_quat.cuh
  gpu-idx-evals/dev-eval-indxing-trifocal_2op1p_30x30.cuh

  #>  p2c constant matrix
  PHC_Coeffs/p2c-5pt_rel_pos_alg_form_quat.h
  PHC_Coeffs/p2c-5pt_rel_pos_geo_form_quat.h
  PHC_Coeffs/p2c-trifocal_2op1p_30x30.h
)

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories( "/users/cchien3/data/cchien3/magma-2.5.4/include" )
include_directories( "/users/cchien3/data/cchien3/magma-2.5.4/testing" )
include_directories( "/users/cchien3/data/cchien3/magma/control" )
include_directories( "/users/cchien3/data/cchien3/magma/magmablas" )

add_library( magmaHC ${control_sources} )
target_link_libraries(magmaHC 
        -L/users/cchien3/data/cchien3/magma-cchien3/lib magma
        -L/gpfs/runtime/opt/openblas/0.3.23-u6k5fey/lib openblas
        )
set_target_properties(magmaHC PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

