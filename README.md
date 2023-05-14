# GPU-HC: Homotopy Continuation Path Tracking in GPU
### Research @ LEMS, Brown University (CVPR 2022)
## Introduction
GPU-HC, as its name suggests, is a GPU implementation of Homotopy Continuation Solver. It is a general tool for finding roots of a polynomial systems (check our papers for more details). Our research aims at applying GPU-HC for computer vision problems, especially multiview geometry problems, where timings are crucial as the computation is typically performed under a RANSAC loop. Please refer and cite the following two papers if you intend to use them in your work. Also, please do not hesitate to contact chiang-heng_chien@brown.edu if you have any questions on using GPU-HC. Currently, the only restriction of GPU-HC is that the polynomial system of interest must have the number of variables lesser than 32. We are planning to develop a solver enabling polynomials of unknowns greater than 32.

This is the source code of the following two papers: <br />
1. ``Chien, Chiang-Heng, Hongyi Fan, Ahmad Abdelfattah, Elias Tsigaridas, Stanimire Tomov, and Benjamin Kimia. "GPU-Based Homotopy Continuation for Minimal Problems in Computer Vision." In Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition, pp. 15765-15776. 2022.`` [Paper link](https://openaccess.thecvf.com/content/CVPR2022/html/Chien_GPU-Based_Homotopy_Continuation_for_Minimal_Problems_in_Computer_Vision_CVPR_2022_paper.html) <br />

2. ``Chien, Chiang-Heng, Hongyi Fan, Elias Tsigaridas, Ahmad Abdelfattah, Stanimire Tomov, and Benjamin Kimia. "Parallel Path Tracking for Homotopy Continuation using GPU." In Proceedings of the International Symposium on Symbolic and Algebraic Computation. 2022.`` [Paper link](https://par.nsf.gov/biblio/10333125) <br /> <br />

## New Updates
1. 2023.05.14 GPU-HC has a new release!
2. 2023.03.12 GPU-HC has now been used in solving three-view relative pose problem of a generalized camera! <br />

## Contents
This repository primarily contains three folders: <br />
``GPU-HC``: main code of GPU-HC solver.
``auto-gen-tools``: software tools used to automatically generate necessary data for GPU-HC to solve new problems <br />
``problem-data-generation``: example polynomial system data. <br />

## Dependencies:
(1) CMake 3.14 or higher <br />
(2) MAGMA 2.6.1 or higher <br />
(3) CUDA 9.0 or higher <br />
(4) cuBlas <br />
(5) openBlas <br />
(6) pthread <br />
(7) if you want to solve a new polynomial problem, you will need Matlab 2019 or higher and [Julia](https://julialang.org/downloads/) in order to generate start parameters and start solutions.

## How to build and compile the code
(1) After cloning the repo, cd to the repo folder and create a 'build' directory and enter it
```bash
mkdir build && cd build
```
(2) make sure to change routes in CMakeFiles according to your local machine. See README.md of either straight-line or parameter HC. <br />
(3) create a make file and compile the entire code <br />
(There are a lot of problems run in this repo. It may be time consuming if all GPU kernels are compiled. To save your time, feel free to comment out the kernel in the CMakeLists.txt under ``/GPU-HC/straight-line-HC/magmaHC/`` or ``/GPU-HC/parameter-HC/magmaHC/``)
```bash
cmake .. && make -j
```

## How to run the execution file after successful compilation
(1) enter the bin foler under build folder of the repo dorectory
```bash
cd bin
```
(2) run the code by specifically typing input arguments. See README.md in folders of either straight-line or parameter HC under ``GPU-HC`` folder.
```bash
./magmaHC-main <input-argument> <command>
```

## How to use GPU-HC to solve a new polynomial problem
