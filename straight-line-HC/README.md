Path Traking for Homotopy Continuation in GPU: Straight-line Homotopy Continuation

# 1. After build and compile

## Input arguments and the associated commands are listed below:
```
<input-argument>      <command>
       -p             <problem-name>         # (or --problem) : name of the polynomial problem
```
For example, run eco12 benchmark polynomial problem by ```./magmaHC-main -p eco12```
If not specified, the help message will be automatically shown.


# 2. Important Changes in Local Machine:
```
(1) ``include_directories" of MAGMA package  in CMakeLists.txt
    ``include_directories" of MAGMA package and ``target_link_libraries" of dependencies  in magmaHC/CMakeLists.txt
(3) The cblas.h file directory in line 18 of magmaHC/cpu-compute/cpu-compute.h
(4) The absolute directories in line 27 in cmd/magmaHC-main.cu
```

# 3. Some Notices:
The sizes of start solutions of katsura20 and katsura21 problems are large, and thus they are not included in this repo. For those who want to run these problems, please create them via Julia.
