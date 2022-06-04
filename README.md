GPU-HC: Path Traking for Homotopy Continuation in GPU for Benchmark Polynomial Problems

# 1. Contents
This repository primarily contains two folders:
``GPU-HC'': GPU implementation of both straight-line and parameter homotopy continuation algorithms
``auto-gen-tools'': software tools used to automatically generate necessary data and code for GPU-HC to solve new problems

# 1. Dependencies:
(1) CMake 3.14 or higher

(2) MAGMA 2.6.1 or higher

(3) CUDA 9.0 or higher

(4) cuBlas

(5) openBlas

(6) pthread

# 2. How to build and compile the code
(1) clone the repo
```bash
git clone https://github.com/C-H-Chien/Homotopy-Continuation-Tracker-on-GPU.git
```
(2) under the repo folder, create a 'build' directory
```bash
mkdir build
```
(3) enter the build folder
```bash
cd build
```

(4) make sure to change routes in CMakeFiles according to your local machine. See README.md of either straight-line or parameter HC.

(5) create a make file
```bash
cmake ..
```
(5) compile the entire code
```bash
make -j
```

# 3. How to run the execution file after successful compilation
(1) enter the bin foler under build folder
```bash
cd bin
```

(2) run the code by specifically typing input arguments. See README.md of either straight-line or parameter HC
```bash
./magmaHC-main <input-argument> <command>
```

# 4. How to use the code to run a new minimal problem
This part will be updated soon.

# 5. Reference
The paper that uses this code will be released as the reference once it is accepted for publication.

