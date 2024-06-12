# GPU-HC: Homotopy Continuation Path Tracking in GPU
### Research @ LEMS, Brown University (CVPR 2022)
## Introduction
GPU-HC, as its name suggests, is a GPU implementation of Homotopy Continuation Solver. It is a general tool for finding roots of a polynomial systems (check our papers for more details). Our research aims at applying GPU-HC for computer vision problems, especially multiview geometry problems, where timings are crucial as the computation is typically performed under a RANSAC loop. Please refer and cite the following two papers if you intend to use them in your work. Also, please do not hesitate to contact chiang-heng_chien@brown.edu if you have any questions on using GPU-HC. <br />

## New Updates
1. 2023.05.14 GPU-HC has a new release! <br />
2. 2023.03.12 GPU-HC has now been used in solving generalized three-view relative pose problem. Checkout this [ICCV paper](https://openaccess.thecvf.com/content/ICCV2023/papers/Ding_Minimal_Solutions_to_Generalized_Three-View_Relative_Pose_Problem_ICCV_2023_paper.pdf) and [GitHub page](https://github.com/C-H-Chien/Three_View_Generalized_Camera). <br />
3. 2024.05.29 Major change in the structure of GPU-HC. The GPU code is also optimized, which provides 1.3-6.0x speedup over the older version. This optimized version will become the second release. <br />

## Contents
This repository primarily contains three folders: <br />
- [x] ``GPU-HC``: main code of GPU-HC solver. <br />
- [x] ``auto-data-gen-tools``: software tools used to automatically generate necessary data for GPU-HC to solve new problems <br />
- [x] ``problem-data-generation``: example polynomial system data. <br />

## Dependencies:
(We are going to release a new version independent of MAGMA library) <br />
(1) CMake 3.14 or higher <br />
(2) MAGMA 2.6.1 or higher <br />
(3) CUDA 9.0 or higher <br />
(4) openBlas 0.3.7 or higher <br />
(5) YAML-CPP, can be built from its [official repo](https://github.com/jbeder/yaml-cpp). <br />
(7) if you want to solve a new polynomial problem, you will need Matlab 2019 or higher and [Julia](https://julialang.org/downloads/) in order to generate start parameters and start solutions.

## Setup
Make sure to change paths of dependencies in CMakeFiles according to your local machine. Follow the standard steps for build and compile, _i.e._, <br />
```bash
mkdir build && cd build
cmake .. && make -j
```
The executable file should appear under ``/buid/bin/`` <br />

## How to run the execution file after successful compilation
(1) Enter the bin foler under build folder of the repo dorectory
```bash
cd bin
```
(2) Run the code by specifically typing the minimal problem name with a ``-p`` tag, _e.g._,
```bash
./magmaHC-main -p generalized_3views_4pts
```
which solves the generalized three-view relative pose problem using 4 points.

## Limitations
Several limitations of the GPU-HC solver: <br />
(1) The system size must not exceed 32x32. We are planning to enable solving such large system size in the future. <br />
(2) GPU-HC is unable to solve an over-determined system. One possible way to tackle this is to choose only partial polynomials to make the system well-determined. This might not guarantee to find the solution of interest, but sometimes it works. <br />

## References
```BibTeX
@InProceedings{@inproceedings{chien2022gpu,
  title={{GPU}-based homotopy continuation for minimal problems in computer vision},
  author={Chien, Chiang-Heng and Fan, Hongyi and Abdelfattah, Ahmad and Tsigaridas, Elias and Tomov, Stanimire and Kimia, Benjamin},
  booktitle={Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition},
  pages={15765--15776},
  year={2022}
}
```
```BibTeX
@inproceedings{chien2022parallel,
  title={Parallel path tracking for homotopy continuation using {GPU}},
  author={Chien, Chiang-Heng and Fan, Hongyi and Abdelfattah, Ahmad and Tsigaridas, Elias and Tomov, Stanimire and Kimia, Benjamin},
  booktitle={Proceedings of the International Symposium on Symbolic and Algebraic Computation},
  year={2022}
}
```
