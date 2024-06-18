# GPU-HC: Homotopy Continuation Path Tracking in GPU
### Research @ LEMS, Brown University (CVPR 2022)
## Introduction
GPU-HC, as its name suggests, is a GPU implementation of Homotopy Continuation Solver. It is a general tool for finding roots of a polynomial systems (check our papers for more details). Our research aims at applying GPU-HC for computer vision problems, especially multiview geometry problems, where timings are crucial as the computation is typically performed under a RANSAC loop. Please refer and cite the reference papers if you intend to use them in your work. Also, please do not hesitate to contact chiang-heng_chien@brown.edu if you have any questions on using GPU-HC. <br />

## :arrow_heading_up: New Updates
1. 2023.05.14 GPU-HC has a new release! <br />
2. 2023.03.12 GPU-HC has now been used in solving generalized three-view relative pose problem. Checkout this [ICCV paper](https://openaccess.thecvf.com/content/ICCV2023/papers/Ding_Minimal_Solutions_to_Generalized_Three-View_Relative_Pose_Problem_ICCV_2023_paper.pdf) and [GitHub page](https://github.com/C-H-Chien/Three_View_Generalized_Camera). <br />
3. 2024.05.29 Major change in the structure of GPU-HC. The GPU code is also optimized, which provides 1.3-6.0x speedup over the older version. This optimized version will become the second release. Refer to the speed in [the table](https://github.com/C-H-Chien/Homotopy-Continuation-Tracker-on-GPU/tree/main?tab=readme-ov-file#timer_clock-speed-in-milliseconds-for-selected-minimal-problems) which is much faster than what was reported in the papers. <br />
4. 2024.06.14 GPU-HC has a second release! Faster in speed and more friendly to new users who aim to build their own GPU-HC solver. <br />

## :floppy_disk: What do we have in this repo?
This repository primarily contains three folders: <br />
:file_folder: ``GPU-HC``: main code of GPU-HC solver. <br />
:file_folder: ``auto-data-gen-tools``: software tools used to automatically generate necessary data for GPU-HC to solve new problems <br />
:file_folder: ``start-sys-gen``: generating start systems by Julia example code. <br />

## :dependabot: Dependencies:
(1) CMake 3.14 or higher <br />
(2) MAGMA 2.6.1 or higher <br />
(3) CUDA 9.0 or higher <br />
(4) openBlas 0.3.7 or higher <br />
(5) YAML-CPP, can be built from its [official repo](https://github.com/jbeder/yaml-cpp). <br />
(7) if you want to solve a new polynomial problem, you will need Matlab 2019 or higher and [Julia](https://julialang.org/downloads/) in order to generate start parameters and start solutions.

## :hammer_and_wrench: Setup and run
Make sure to change paths of dependencies in CMakeFiles according to your local machine. Follow the standard steps for build and compile, _i.e._, <br />
```bash
mkdir build && cd build
cmake .. && make -j
```
The executable file should appear under ``/buid/bin/`` <br />
```bash
cd bin
```
Run the code by specifically typing the minimal problem name with a ``-p`` tag, _e.g._,
```bash
./magmaHC-main -p generalized_3views_4pts
```
which solves the generalized three-view relative pose problem using 4 points.

## :beginner: How to solve a new minimal problem
Refer to [Solve_New_Problems](https://github.com/C-H-Chien/Homotopy-Continuation-Tracker-on-GPU/blob/main/Solve_New_Problems.md) document and follow a step-by-step instruction.

## :timer_clock: Speed in milliseconds for selected minimal problems
E.T. means Elimination Template which has an out of memory issue ("X" in the table) for large, hard minimal problems (the first 7 rows). The speed of CPU-HC is chosen from the fastest among Julia, [MiNuS](https://github.com/rfabbri/minus), and my own implementation, running on 8-core multi-threading. Refer to my [another repository](https://github.com/C-H-Chien/Minimal-Problem-Solver-on-CPU) for how you could use elimination template, Julia, and my CPU-HC for solving minimal problems. The selected problems below are provided in this repository under ``GPU-HC/problems/`` with references given in respective YAML files.
| Problem           | # of Unkn. | # of Sol. | E.T. | CPU-HC | GPU-HC |
| :---------------: | :------: | :----: | :-------: | :-------: | :-------: |
| trifocal relative pose, unkown focal length              | 12 |  180 | X | 152.83 | **15.17** | 
| trifocal relative pose from lines at points (30x30)      | 30 |  312 | X | 234.86 | **54.44** |
| generalized 3-view relative pose from 4 points           | 12 |  583 | X | 428.26 |  **4.77** |
| generalized 3-view relative pose from 6 lines            |  6 |  600 | X | 1103   |  **5.11** |
| 5-point relative pose with depths                        | 16 |   40 | X |  27.46 |  **7.26** |
| 6-point rolling-shutter absolute pose (linearized)       | 18 |   40 | X | 145.79 |  **8.92** |
| 4-view triangulation                                     | 11 |  142 | X |  48.47 |  **1.74** |
| 3-view triangulation                                     |  9 |  142 | 612.43 | 17.84 | **0.88** |
| optimal PnP (quaternion)                                 |  4 |   80 |  36.33 | 55.91 | **0.77** |
| dual-receiver TDOA 5-points                              |  5 |   24 |   9.03 | 22.04 | **3.17** |
| distorted 2-view triangulation                           |  5 |   28 |   5.92 |  7.83 | **0.83** |
| 5-point relative pose (reduced/algebraic form)           |  6 |   40 |   2.97 | 24.69 | **1.77** |
| optimal P4P absolute pose                                |  5 |   32 |   1.86 |  3.92 | **0.43** |
| 6-point rolling shutter absolute pose                    |  6 |   20 |   1.58 | 37.30 | **0.65** |
| 3-point relative pose with homography constraint         |  8 |    8 |   1.47 |  1.72 | **0.84** |
| PnP with unknown principal point                         | 10 |    8 |   **1.47** |  6.40 | 3.45 |
| relative pose using quiver, unknown focal length         |  4 |   20 |   1.08 |  2.97 | **0.75** |
| P3P, absolute pose                                       |  3 |    8 |   **0.06** |  1.18 | 0.15 |
| 5-point relative pose (right null-space)                 |  3 |    8 |   **0.04** | 25.76 | 0.61 |

## Limitations
Several limitations of the GPU-HC solver: <br />
(1) The system size must not exceed 32x32. We are planning to enable solving such large system size in the future. <br />
(2) GPU-HC is unable to solve an over-determined system. One possible way to tackle this is to choose only partial polynomials to make the system well-determined. This might not guarantee to find the solution of interest, but sometimes it works. <br />

## :bookmark: References
```BibTeX
@InProceedings{chien2022gpu,
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
