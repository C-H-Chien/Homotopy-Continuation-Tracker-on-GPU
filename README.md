Path Traking for Homotopy Continuation in GPU for Benchmark Polynomial Problems

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

# 4. Reference
This code is made by the work proposed in the following paper:

Chiang-Heng Chien, Hongyi Fan, Elias Tsigaridas, Ahmad Abdelfattah,
Stanimire Tomov and Benjamin Kimia, "Parallel Path Tracking for Ho-
motopy Continuation using GPU," submitted to International Symposium
on Symbolic and Algebraic Computation (ISSAC’22), 2022.

