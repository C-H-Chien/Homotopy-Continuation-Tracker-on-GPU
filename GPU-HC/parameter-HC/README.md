Path Traking for Homotopy Continuation in GPU - Parameter Homotopy Continuation applied to Computer Vision Applications

# 1. After build and compile

## Input arguments and the associated command:
```
<input-argument>      <command>
       -p             <problem-name>         # (or --problem) : name of the polynomial problem
```
For example, run 3-view Triangulation problem by ```./magmaHC-main -p 3vTrg```

If not specified, the help message will be automatically shown.


# 2. Important Changes in Local Machine:
```
(1) ``include_directories" of MAGMA package  in CMakeLists.txt
    ``include_directories" of MAGMA package and ``target_link_libraries" of dependencies  in magmaHC/CMakeLists.txt
(3) The cblas.h file directory in line 18 of magmaHC/cpu-compute/cpu-compute.h
(4) The absolute directories in line 27 in cmd/magmaHC-main.cu
```
