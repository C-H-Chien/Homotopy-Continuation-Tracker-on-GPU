
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
  gpu-kernels/kernel_HC_Solver_5pt_rel_pos_full_form_quat.cu
  gpu-kernels/kernel_HC_Solver_trifocal_2op1p_30x30.cu
  gpu-kernels/kernel_HC_Solver_generalized_3views_4pts.cu
  gpu-kernels/kernel_HC_Solver_generalized_3views_6lines.cu
  gpu-kernels/kernel_HC_Solver_uncalibrated_trifocal_rel_pos_CY.cu
  gpu-kernels/kernel_HC_Solver_uncalibrated_trifocal_rel_pos_CH.cu
  gpu-kernels/kernel_HC_Solver_optimal_PnP_quat.cu
  gpu-kernels/kernel_HC_Solver_3view_triangulation.cu
  gpu-kernels/kernel_HC_Solver_4view_triangulation.cu
  gpu-kernels/kernel_HC_Solver_6pt_RS_abs_pos.cu
  gpu-kernels/kernel_HC_Solver_6pt_RS_abs_pos_1lin.cu
  gpu-kernels/kernel_HC_Solver_dual_reciever_TDOA_5pt.cu
  gpu-kernels/kernel_HC_Solver_distorted_2view_triangulation.cu
  gpu-kernels/kernel_HC_Solver_optimal_P4P_abs_pos.cu
  gpu-kernels/kernel_HC_Solver_3pt_rel_pos_homog.cu
  gpu-kernels/kernel_HC_Solver_PnP_unkn_principal_pt.cu
  gpu-kernels/kernel_HC_Solver_rel_pos_quiver.cu
  gpu-kernels/kernel_HC_Solver_P3P.cu

  #>  device function
  dev-cgesv-batched-small.cuh
  dev-cgesv-batched-small-sm.cuh
  dev-get-new-data.cuh
  gpu-idx-evals/dev-eval-indxing-5pt_rel_pos_geo_form_quat.cuh
  gpu-idx-evals/dev-eval-indxing-5pt_rel_pos_alg_form_quat.cuh
  gpu-idx-evals/dev-eval-indxing-5pt_rel_pos_full_form_quat.cuh
  gpu-idx-evals/dev-eval-indxing-trifocal_2op1p_30x30.cuh
  gpu-idx-evals/dev-eval-indxing-generalized_3views_4pts.cuh
  gpu-idx-evals/dev-eval-indxing-generalized_3views_6lines.cuh
  gpu-idx-evals/dev-eval-indxing-uncalibrated_trifocal_rel_pos_CY.cuh
  gpu-idx-evals/dev-eval-indxing-uncalibrated_trifocal_rel_pos_CH.cuh
  gpu-idx-evals/dev-eval-indxing-optimal_PnP_quat.cuh
  gpu-idx-evals/dev-eval-indxing-3view_triangulation.cuh
  gpu-idx-evals/dev-eval-indxing-4view_triangulation.cuh
  gpu-idx-evals/dev-eval-indxing-6pt_RS_abs_pos.cuh
  gpu-idx-evals/dev-eval-indxing-6pt_RS_abs_pos_1lin.cuh
  gpu-idx-evals/dev-eval-indxing-dual_reciever_TDOA_5pt.cuh
  gpu-idx-evals/dev-eval-indxing-distorted_2view_triangulation.cuh
  gpu-idx-evals/dev-eval-indxing-optimal_P4P_abs_pos.cuh
  gpu-idx-evals/dev-eval-indxing-3pt_rel_pos_homog.cuh
  gpu-idx-evals/dev-eval-indxing-PnP_unkn_principal_pt.cuh
  gpu-idx-evals/dev-eval-indxing-rel_pos_quiver.cuh
  gpu-idx-evals/dev-eval-indxing-P3P.cuh

  #>  P2C
  PHC_Coeffs/p2c-5pt_rel_pos_alg_form_quat.h
  PHC_Coeffs/p2c-5pt_rel_pos_geo_form_quat.h
  PHC_Coeffs/p2c-5pt_rel_pos_full_form_quat.h
  PHC_Coeffs/p2c-trifocal_2op1p_30x30.h
  PHC_Coeffs/p2c-generalized_3views_4pts.h
  PHC_Coeffs/p2c-numerical_generalized_3views_6lines.h
  PHC_Coeffs/p2c-uncalibrated_trifocal_rel_pos_CY.h
  PHC_Coeffs/p2c-uncalibrated_trifocal_rel_pos_CH.h
  PHC_Coeffs/p2c-optimal_PnP_quat.h
  PHC_Coeffs/p2c-3view_triangulation.h
  PHC_Coeffs/p2c-4view_triangulation.h
  PHC_Coeffs/p2c-6pt_RS_abs_pos.h
  PHC_Coeffs/p2c-6pt_RS_abs_pos_1lin.h
  PHC_Coeffs/p2c-dual_reciever_TDOA_5pt.h
  PHC_Coeffs/p2c-distorted_2view_triangulation.h
  PHC_Coeffs/p2c-optimal_P4P_abs_pos.h
  PHC_Coeffs/p2c-3pt_rel_pos_homog.h
  PHC_Coeffs/p2c-PnP_unkn_principal_pt.h
  PHC_Coeffs/p2c-rel_pos_quiver.h
  PHC_Coeffs/p2c-P3P.h

  #> Utilities
  util.hpp
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
        -L/gpfs/data/bkimia/cchien3/yaml-cpp/bin/lib64 yaml-cpp
        )
set_target_properties(magmaHC PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

