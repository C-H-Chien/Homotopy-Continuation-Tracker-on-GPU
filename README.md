# GPU-HC: Homotopy Continuation Path Tracking in GPU
### Research @ LEMS, Brown University (CVPR 2022)
## Introduction
GPU-HC, as its name suggests, is a GPU implementation of Homotopy Continuation Solver. It is a general tool for finding roots of a polynomial systems (check our papers for more details). Our research aims at applying GPU-HC for computer vision problems, especially multiview geometry problems, where timings are crucial as the computation is typically performed under a RANSAC loop. Please refer and cite the following two papers if you intend to use them in your work. Also, please do not hesitate to contact chiang-heng_chien@brown.edu if you have any questions on using GPU-HC. Currently, the only restriction of GPU-HC is that the polynomial system of interest must have the number of variables lesser than 32. We are planning to develop a solver enabling polynomials of unknowns greater than 32. <br />

1. ``Chien, Chiang-Heng, Hongyi Fan, Ahmad Abdelfattah, Elias Tsigaridas, Stanimire Tomov, and Benjamin Kimia. "GPU-Based Homotopy Continuation for Minimal Problems in Computer Vision." In Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition, pp. 15765-15776. 2022.`` [[Paper link](https://openaccess.thecvf.com/content/CVPR2022/html/Chien_GPU-Based_Homotopy_Continuation_for_Minimal_Problems_in_Computer_Vision_CVPR_2022_paper.html)] <br />

2. ``Chien, Chiang-Heng, Hongyi Fan, Elias Tsigaridas, Ahmad Abdelfattah, Stanimire Tomov, and Benjamin Kimia. "Parallel Path Tracking for Homotopy Continuation using GPU." In Proceedings of the International Symposium on Symbolic and Algebraic Computation. 2022.`` [[Paper link](https://par.nsf.gov/biblio/10333125)] <br /> <br />

## New Updates
1. 2023.05.14 GPU-HC has a new release! <br />
2. 2023.03.12 GPU-HC has now been used in solving generalized three-view relative pose problem. Checkout this [ICCV paper](https://openaccess.thecvf.com/content/ICCV2023/papers/Ding_Minimal_Solutions_to_Generalized_Three-View_Relative_Pose_Problem_ICCV_2023_paper.pdf) and [GitHub page](https://github.com/C-H-Chien/Three_View_Generalized_Camera). <br />

## Contents
This repository primarily contains three folders: <br />
- [x] ``GPU-HC``: main code of GPU-HC solver. <br />
- [x] ``auto-data-gen-tools``: software tools used to automatically generate necessary data for GPU-HC to solve new problems <br />
- [x] ``problem-data-generation``: example polynomial system data. <br />

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
(2) Make sure to change paths of dependencies in CMakeFiles according to your local machine. <br />
(3) Create a make file and compile the entire code <br />
```bash
cmake .. && make -j
```

## How to run the execution file after successful compilation
(1) Enter the bin foler under build folder of the repo dorectory
```bash
cd bin
```
(2) Run the code by specifically typing input arguments.
```bash
./magmaHC-main <input-argument> <command>
```
As an example, to run 5-point relative pose problem in geometric form, type
```bash
./magmaHC-main -p 5pt_rel_pos_geo_form_quat
```

## How to use GPU-HC to solve a new polynomial problem
There are two example problems provided in this repository: 5-point relative pose problem in (i) geometric and (ii) algebraic forms. Their polynomials can be generated from either the Julia script or the MATLAB script under ``problem-data-generation/`` <br />
- **Step 1. Create the polynomial system:** Formulate your _explicit_ polynomial system as a matlab script. An example can be found in ``auto-data-gen-tools/sys_5pt_rel_pos_geo_form_quat.m`` where ``p`` is the system parameters and ``x`` are the unknowns. 
- **Step 2. Generate start parameters and solutions:** Create a Julia script which perfoms a Monodromy solver that finds the solutions of a start system. An example can be found under ``problem-data-generation/``, _e.g._, for 5-point relative pose estimation of geometric form, refer to ``Julia_Monodromy_Solver_Examples/5-Point-Relative-Pose-Geometric-Form.jl``. Remember to write the start parameters and solutions to files.
- **Step 3. Reformulate the start parameters and solutions:** <br />
	- **Start solutions:** Run ``auto-data-gen-tools/reformatStartSolsFromJulia.m`` to reformulate the start solutions. Make sure to specify the file path and names. <br />
	- **Start parameters:** Reformulation is trivial. Manually reformulate it as the one in ``GPU-HC/problems/5pt_rel_pos_geo_form_quat/start_params.txt``. The first column is the real part while the second is the imaginary part. <br />


<br />
(To be updated...) <br />

## Limitations
Several limitations of the GPU-HC solver: <br />
(1) The system size must not exceed 32x32. <br />
(2) GPU-HC is unable to solve an over-determined system. One possible way to tackle this is to choose only partial polynomials to make the system well-determined. This might not guarantee to find the solution of interest, but sometimes it works. <br />
(3) There is no _generic_ (_e.g._, circular arc) gamma trick applied in the solver. This will however be introduced in the future.
