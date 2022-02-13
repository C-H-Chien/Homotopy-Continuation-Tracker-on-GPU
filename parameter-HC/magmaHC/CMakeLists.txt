
set (control_sources
  # -- input info --
  define_params_and_read_files.cu

  # -- main gpu code --
  magmaHC-problems.cuh
  homotopy_continuation_solver.cu

  # -- magma kernels --
  gpu-kernels/magmaHC-kernels.h
  gpu-kernels/kernel_HC_Solver_3view_unknownf_pHC.cu
  gpu-kernels/kernel_HC_Solver_3vTrg.cu
  gpu-kernels/kernel_HC_Solver_3vTrg_relax.cu
  gpu-kernels/kernel_HC_Solver_4vTrg.cu
  gpu-kernels/kernel_HC_Solver_5pt_rel_pose_w_depth_recon.cu
  gpu-kernels/kernel_HC_Solver_optimalPnP_w_quaternion.cu
  gpu-kernels/kernel_HC_Solver_3pt_rel_pose_w_homo_constraint.cu
  gpu-kernels/kernel_HC_Solver_r6p.cu
  gpu-kernels/kernel_HC_Solver_refractive_p5p.cu
  gpu-kernels/kernel_HC_Solver_refractive_p6p.cu

  # -- device function --
  dev-cgesv-batched-small.cuh
  dev-get-new-data.cuh
  gpu-idx-evals/dev-eval-indxing-3view_unknownf_pHC.cuh
  gpu-idx-evals/dev-eval-indxing-3vTrg.cuh
  gpu-idx-evals/dev-eval-indxing-3vTrg_relax.cuh
  gpu-idx-evals/dev-eval-indxing-4vTrg.cuh
  gpu-idx-evals/dev-eval-indxing-5pt_rel_pose_w_depth_recon.cuh
  gpu-idx-evals/dev-eval-indxing-optimalPnP_w_quaternion.cuh
  gpu-idx-evals/dev-eval-indxing-3pt_rel_pose_w_homo_constraint.cuh
  gpu-idx-evals/dev-eval-indxing-r6p.cuh
  gpu-idx-evals/dev-eval-indxing-refractive_p5p.cuh
  gpu-idx-evals/dev-eval-indxing-refractive_p6p.cuh

  # -- cpu computation --
  cpu-compute/cpu-compute.h
  cpu-compute/cpu-hc-solver/cpu-hc-solver-3view_unknownf_pHC.cpp
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-3view_unknownf_pHC.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-3view_unknownf_pHC.h
  cpu-compute/cpu-hc-solver/cpu-hc-solver-3vTrg.cpp
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-3vTrg.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-3vTrg.h
  cpu-compute/cpu-hc-solver/cpu-hc-solver-3vTrg_relax.cpp
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-3vTrg_relax.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-3vTrg_relax.h
  cpu-compute/cpu-hc-solver/cpu-hc-solver-4vTrg.cpp
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-4vTrg.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-4vTrg.h
  cpu-compute/cpu-hc-solver/cpu-hc-solver-5pt_rel_pose_w_depth_recon.cpp
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-5pt_rel_pose_w_depth_recon.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-5pt_rel_pose_w_depth_recon.h
  cpu-compute/cpu-hc-solver/cpu-hc-solver-optimalPnP_w_quaternion.cpp
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-optimalPnP_w_quaternion.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-optimalPnP_w_quaternion.h
  cpu-compute/cpu-hc-solver/cpu-hc-solver-3pt_rel_pose_w_homo_constraint.cpp
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-3pt_rel_pose_w_homo_constraint.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-3pt_rel_pose_w_homo_constraint.h
  cpu-compute/cpu-hc-solver/cpu-hc-solver-r6p.cpp
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-r6p.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-r6p.h
  cpu-compute/cpu-hc-solver/cpu-hc-solver-refractive_p5p.cpp
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-refractive_p5p.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-refractive_p5p.h
  cpu-compute/cpu-hc-solver/cpu-hc-solver-refractive_p6p.cpp
  cpu-compute/cpu-eval-HxHt/cpu-eval-HxHt-refractive_p6p.h
  cpu-compute/cpu-eval-HxH/cpu-eval-HxH-refractive_p6p.h

  # -- const matrices definition --
  const-matrices/define-const-matrices.h
  const-matrices/const-matrices-allocations.h

  # -- p2c constant matrix --
  const-matrices/p2c-3view_unknownf_pHC.h
  const-matrices/p2c-3vTrg.h
  const-matrices/p2c-3vTrg_relax.h
  const-matrices/p2c-4vTrg.h
  const-matrices/p2c-5pt_rel_pose_w_depth_recon.h
  const-matrices/p2c-optimalPnP_w_quaternion.h
  const-matrices/p2c-3pt_rel_pose_w_homo_constraint.h
  const-matrices/p2c-r6p.h
  const-matrices/p2c-refractive_p5p.h
  const-matrices/p2c-refractive_p6p.h
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
