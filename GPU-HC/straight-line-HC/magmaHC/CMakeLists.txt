
set (control_sources  
  # -- input info --
  define_params_and_read_files.cu
  
  # -- main gpu code --
  magmaHC-problems.cuh 
  homotopy_continuation_solver.cu
  
  # -- magma kernels --
  gpu-kernels/magmaHC-kernels.h
  gpu-kernels/kernel_HC_Solver_alea6.cu
  gpu-kernels/kernel_HC_Solver_eco12.cu
  gpu-kernels/kernel_HC_Solver_d1.cu
  gpu-kernels/kernel_HC_Solver_cyclic7.cu
  gpu-kernels/kernel_HC_Solver_cyclic8.cu
  gpu-kernels/kernel_HC_Solver_cyclic9.cu
  gpu-kernels/kernel_HC_Solver_katsura6.cu
  gpu-kernels/kernel_HC_Solver_katsura8.cu
  gpu-kernels/kernel_HC_Solver_katsura7.cu
  gpu-kernels/kernel_HC_Solver_katsura9.cu
  gpu-kernels/kernel_HC_Solver_katsura10.cu
  gpu-kernels/kernel_HC_Solver_katsura11.cu
  gpu-kernels/kernel_HC_Solver_katsura12.cu
  gpu-kernels/kernel_HC_Solver_katsura13.cu
  gpu-kernels/kernel_HC_Solver_katsura14.cu
  gpu-kernels/kernel_HC_Solver_katsura15.cu
  gpu-kernels/kernel_HC_Solver_game6two.cu
  gpu-kernels/kernel_HC_Solver_game7two.cu
  #gpu-kernels/kernel_HC_Solver_pole28sys.cu
  
  # -- device function --
  dev-cgesv-batched-small.cuh
  dev-get-new-data.cuh
  gpu-idx-eval/dev-eval-indxing-alea6.cuh
  gpu-idx-eval/dev-eval-indxing-eco12.cuh
  gpu-idx-eval/dev-eval-indxing-d1.cuh
  gpu-idx-eval/dev-eval-indxing-cyclic7.cuh
  gpu-idx-eval/dev-eval-indxing-cyclic8.cuh
  gpu-idx-eval/dev-eval-indxing-cyclic9.cuh
  gpu-idx-eval/dev-eval-indxing-game6two.cuh
  gpu-idx-eval/dev-eval-indxing-game7two.cuh
  gpu-idx-eval/dev-eval-indxing-katsura6.cuh
  gpu-idx-eval/dev-eval-indxing-katsura7.cuh
  gpu-idx-eval/dev-eval-indxing-katsura8.cuh
  gpu-idx-eval/dev-eval-indxing-katsura9.cuh
  gpu-idx-eval/dev-eval-indxing-katsura10.cuh
  gpu-idx-eval/dev-eval-indxing-katsura11.cuh
  gpu-idx-eval/dev-eval-indxing-katsura12.cuh
  gpu-idx-eval/dev-eval-indxing-katsura13.cuh
  gpu-idx-eval/dev-eval-indxing-katsura14.cuh
  gpu-idx-eval/dev-eval-indxing-katsura15.cuh
  #gpu-idx-eval/dev-eval-indxing-pole28sys.cuh
  
  # -- CPU-HC --
  # -- 1) solver --
  cpu-compute/cpu-compute.h
  cpu-compute/cpu-hc-solver/cpu-hc-solver-alea6.cpp
  cpu-compute/cpu-hc-solver/cpu-hc-solver-eco12.cpp
  cpu-compute/cpu-hc-solver/cpu-hc-solver-cyclic7.cpp
  cpu-compute/cpu-hc-solver/cpu-hc-solver-cyclic8.cpp
  cpu-compute/cpu-hc-solver/cpu-hc-solver-cyclic9.cpp
  cpu-compute/cpu-hc-solver/cpu-hc-solver-d1.cpp
  cpu-compute/cpu-hc-solver/cpu-hc-solver-katsura6.cpp
  cpu-compute/cpu-hc-solver/cpu-hc-solver-katsura7.cpp
  cpu-compute/cpu-hc-solver/cpu-hc-solver-katsura8.cpp
  cpu-compute/cpu-hc-solver/cpu-hc-solver-katsura9.cpp
  cpu-compute/cpu-hc-solver/cpu-hc-solver-katsura10.cpp
  cpu-compute/cpu-hc-solver/cpu-hc-solver-katsura11.cpp
  cpu-compute/cpu-hc-solver/cpu-hc-solver-katsura12.cpp
  cpu-compute/cpu-hc-solver/cpu-hc-solver-katsura13.cpp
  cpu-compute/cpu-hc-solver/cpu-hc-solver-katsura14.cpp
  cpu-compute/cpu-hc-solver/cpu-hc-solver-katsura15.cpp
  cpu-compute/cpu-hc-solver/cpu-hc-solver-game6two.cpp
  cpu-compute/cpu-hc-solver/cpu-hc-solver-game7two.cpp
  cpu-compute/cpu-hc-solver/cpu-hc-solver-pole28sys.cpp
  # -- 2) HxHt evaluations --
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-alea6.h
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-cyclic7.h
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-cyclic8.h
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-cyclic9.h
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-d1.h
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-eco12.h
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-katsura6.h
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-katsura7.h
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-katsura8.h
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-katsura9.h
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-katsura10.h
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-katsura11.h
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-katsura12.h
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-katsura13.h
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-katsura14.h
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-katsura15.h
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-game6two.h
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-game7two.h
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-pole28sys.h
  # -- 3) HxH evaluations --
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-alea6.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-d1.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-eco12.h 
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-cyclic7.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-cyclic8.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-cyclic9.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-katsura6.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-katsura7.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-katsura8.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-katsura9.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-katsura10.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-katsura11.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-katsura12.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-katsura13.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-katsura14.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-katsura15.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-game6two.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-game7two.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-pole28sys.h
  
  # -- const matrices definition --
  const-matrices/define-const-matrices.h
  const-matrices/const-matrices-allocations.h
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
#                                       --generate-code arch=compute_75,code=sm_75
#                                       --maxrregcount 255
#                                       >)

