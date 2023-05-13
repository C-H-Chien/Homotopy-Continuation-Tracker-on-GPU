Path Traking for Homotopy Continuation in GPU - Parameter Homotopy Continuation applied to Computer Vision Applications

# 1. After build and compile

## Input arguments and the associated command:
```
<input-argument>      <command>
       -p             <problem-name>         # (or --problem) : name of the polynomial problem
```
For example, run 5-point relative pose minimal problem (geometric form) ```./magmaHC-main -p 5pt_rel_pos_geo_form_quat```

If not specified, the help message will be automatically shown.


# 2. Important Changes in Local Machine:
```
(1) ``include_directories" of MAGMA libraries in CMakeLists.txt
    ``include_directories" of MAGMA libraries and ``target_link_libraries" of dependencies in magmaHC/CMakeLists.txt
(2) The global repo directory in cmd/magmaHC-main.cu
```
