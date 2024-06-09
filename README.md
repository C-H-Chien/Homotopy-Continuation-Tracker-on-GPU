# GPU-HC: Homotopy Continuation Path Tracking in GPU
### Research @ LEMS, Brown University (CVPR 2022)
## Introduction
GPU-HC, as its name suggests, is a GPU implementation of Homotopy Continuation Solver. It is a general tool for finding roots of a polynomial systems (check our papers for more details). Our research aims at applying GPU-HC for computer vision problems, especially multiview geometry problems, where timings are crucial as the computation is typically performed under a RANSAC loop. Please refer and cite the following two papers if you intend to use them in your work. Also, please do not hesitate to contact chiang-heng_chien@brown.edu if you have any questions on using GPU-HC. <br />

## New Updates
1. 2023.05.14 GPU-HC has a new release! <br />
2. 2023.03.12 GPU-HC has now been used in solving generalized three-view relative pose problem. Checkout this [ICCV paper](https://openaccess.thecvf.com/content/ICCV2023/papers/Ding_Minimal_Solutions_to_Generalized_Three-View_Relative_Pose_Problem_ICCV_2023_paper.pdf) and [GitHub page](https://github.com/C-H-Chien/Three_View_Generalized_Camera). <br />
3. 2023.12.29 Major change in the GPU-HC main code which homogenize the code for all minimal problems. Also adding gamma trick for circular-arc homotopy. <br />

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
(4) openBlas <br />
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
(2) Run the code by specifically typing input arguments.
```bash
./magmaHC-main <input-argument> <command>
```
As an example, to run 5-point relative pose problem in geometric form, type
```bash
./magmaHC-main -p 5pt_rel_pos_geo_form_quat
```

## How to use GPU-HC to solve a new polynomial problem
There are two example problems provided in this repository: 5-point relative pose problem in __(i)__ geometric and __(ii)__ algebraic forms. Their polynomials can be generated from either the Julia script or the MATLAB script under ``problem-data-generation/``. Let your polynomial system name be <problem-name>. <br />
- **Step 1. Create Your Polynomial System:** Formulate your _explicit_ polynomial system as a matlab script. An example can be found in ``auto-data-gen-tools/sys_5pt_rel_pos_geo_form_quat.m`` where ``p`` is the system parameters and ``x`` are the unknowns. 
- **Step 2. Generate Start Parameters and Solutions:** 
	- Create a Julia script which perfoms a Monodromy solver that finds the solutions of a start system. An example can be found under ``problem-data-generation/``, _e.g._, for 5-point relative pose estimation of geometric form, refer to ``Julia_Monodromy_Solver_Examples/5-Point-Relative-Pose-Geometric-Form.jl``. <br />
	- Remember to write the start parameters and solutions to files when running Julia. <br />
- **Step 3. Reformulate the Start Parameters and Solutions:** <br />
	- **Start solutions:** Run ``auto-data-gen-tools/reformatStartSolsFromJulia.m`` to reformulate the start solutions. Make sure to specify the file path and names. <br />
	- **Start parameters:** Reformulation is trivial. Manually reformulate it as the one in ``GPU-HC/problems/5pt_rel_pos_geo_form_quat/start_params.txt``. The first column is the real part while the second is the imaginary part. <br />
- **Step 4. Create Jacobian Matrices Evaluation Indices Automatically:** <br />
	- Add your problem in the ``auto-data-gen-tools/params2coeffs.m`` (refer to lines 60-64). This needs the matlab script you prepared in Step 1. <br />
	- Run the ``auto-data-gen-tools/params2coeffs.m`` script. You should expect several output files, but the following files are necessary for GPU-HC: <br />
		- ``Hx_idx.txt``: Jacobian matrix $\frac{\partial H}{\partial x}$ evaluation indices. <br />
		- ``Ht_idx.txt``: Jacobian matrix $\frac{\partial H}{\partial t}$ evaluation indices. <br />
		- ``P2C_script.txt``: Parameters to coefficient conversion.
- **Step 5. Deploy your problem in GPU-HC:** 
	- Make a folder named <problem-name> under ``GPU-HC/problems/`` and put __(i)__ start system solutions and parameters made in Step 3, __(ii)__ ``Hx_idx.txt`` and ``Ht_idx.txt`` made in Step 4, and __(iii)__ your target parameters for which the formulation is identical to the start parameters.
	- Create the following files for your polynomial problem:
		- ``GPU-HC/magmaHC/const-matrices/<problem-name>.h``. Paste the content in ``P2C_script.txt`` into it and make necessary revisions.
		- ``GPU-HC/gpu-kernels/<problem-name>.cu``. You can make it as a copy from the example and make necessary revisions.
		- ``GPU-HC/gpu-idx-evals/dev-eval-indxing-<problem-name>.cuh``. Again, you can make it as a copy from the example and make necessary revisions.
	- Insert your polynomial problem:
		- ``GPU-HC/magmaHC/define_params_and_read_files.cu``. Numbers for Hx_maximal_terms, etc. are shown after running ``params2coeffs.m`` in Step 4.
		- ``GPU-HC/magmaHC/homotopy_continuation_solver.cu``, lines 181-192.
		- ``GPU-HC/cmd/magmaHC-main.cu``, lines 223-224.
	- Edit the following scripts:
		- ``GPU-HC/gpu-kernels/<problem-name>.cu``: input arguments when launching the kernel are: <br />

		- ``GPU-HC/gpu-idx-evals/dev-eval-indxing-<problem-name>.cuh``:
			- Function ``eval_parameter_homotopy``: Check the comments from the example script ``dev-eval-indxing-5pt_rel_pos_geo_form_quat.cuh`` to make necessary changes.
			- Function ``eval_Jacobian_Hx``: the number of ``s_track`` is the value ``Hx_maximal_parts`` minus 1.
			- Function ``eval_Jacobian_Ht``: the number of ``s_track`` is the value ``Ht_maximal_parts`` minus 1.
			- Function ``eval_Homotopy``: the same as ``eval_Jacobian_Ht``.

<br />
(To be updated...) <br />

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
