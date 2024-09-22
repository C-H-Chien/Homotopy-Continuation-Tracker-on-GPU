
set (control_sources
  any_class.cpp
  
)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_library( magmaHC ${control_sources} )
target_link_libraries(magmaHC
  -L/apps/spacks/2023-08-14/opt/spack/linux-ubuntu23.04-x86_64/gcc-9.5.0/magma-2.7.1-5hrvjybyve3kqfgx7uv5nxxz7mktnqqw/lib magma
  -L/apps/spacks/2023-08-14/opt/spack/linux-ubuntu23.04-x86_64/gcc-9.5.0/openblas-0.3.23-2tips3wywhkieahh6udqmjfke4flqm3l/lib openblas
  -L/home/chchien/thirdParty/yaml-cpp/yaml-cpp-bin/lib64 yaml-cpp
)
#set_target_properties(magmaHC PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

