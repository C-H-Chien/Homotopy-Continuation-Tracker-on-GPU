
set (control_sources
  #>  input info
  define_params_and_read_files.cu

  #>  main solver code
  magmaHC-problems.cuh
  homotopy_continuation_solver.cu

  #>  magma kernels
  gpu-kernels/magmaHC-kernels.h
  gpu-kernels/kernel_HC_Solver_5pt_rel_pos_geo_form_quat.cu
  gpu-kernels/kernel_HC_Solver_5pt_rel_pos_alg_form_quat.cu

  #>  device function
  dev-cgesv-batched-small.cuh
  dev-get-new-data.cuh
  gpu-idx-evals/dev-eval-indxing-5pt_rel_pos_geo_form_quat.cuh
  gpu-idx-evals/dev-eval-indxing-5pt_rel_pos_alg_form_quat.cuh

  #>  p2c constant matrix
  const-matrices/p2c-5pt_rel_pos_alg_form_quat.h
  const-matrices/p2c-5pt_rel_pos_geo_form_quat.h
)

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories( "/users/cchien3/data/cchien3/magma-2.5.4/include" )
include_directories( "/users/cchien3/data/cchien3/magma-2.5.4/testing" )
include_directories( "/users/cchien3/data/cchien3/magma/control" )
include_directories( "/users/cchien3/data/cchien3/magma/magmablas" )
add_library( magmaHC ${control_sources} )
target_link_libraries(magmaHC 
	PUBLIC OpenMP::OpenMP_CXX
        -L/gpfs/runtime/opt/magma-gpu/2.4.0/lib magma
        -L${CUDAPATH}/lib cublas cudart cusparse
        -L/gpfs/runtime/opt/openblas/0.3.7/lib openblas
        -L/usr/lib64 pthread
        )
set_target_properties(magmaHC PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

#target_compile_options(magmaHC PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
#                                      generate-code arch=compute_75,code=sm_75
#                                      maxrregcount 255
#                                       >)
