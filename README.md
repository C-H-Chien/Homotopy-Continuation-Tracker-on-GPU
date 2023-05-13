# GPU-HC: Homotopy Continuation Path Tracking in GPU

# 0. References
This is the source code of the following two papers: <br />
1. ``Chien, Chiang-Heng, Hongyi Fan, Ahmad Abdelfattah, Elias Tsigaridas, Stanimire Tomov, and Benjamin Kimia. "GPU-Based Homotopy Continuation for Minimal Problems in Computer Vision." In Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition, pp. 15765-15776. 2022.`` [Paper link](https://openaccess.thecvf.com/content/CVPR2022/html/Chien_GPU-Based_Homotopy_Continuation_for_Minimal_Problems_in_Computer_Vision_CVPR_2022_paper.html) <br /> <br />

2. ``Chien, Chiang-Heng, Hongyi Fan, Elias Tsigaridas, Ahmad Abdelfattah, Stanimire Tomov, and Benjamin Kimia. "Parallel Path Tracking for Homotopy Continuation using GPU." In Proceedings of the International Symposium on Symbolic and Algebraic Computation. 2022.`` [Paper link](https://par.nsf.gov/biblio/10333125) <br /> <br />

# 1. Contents
This repository primarily contains three folders: <br />
``GPU-HC``: GPU implementation of both straight-line and parameter homotopy continuation algorithms <br />
``auto-gen-tools``: software tools used to automatically generate necessary data and code for GPU-HC to solve new problems <br />
``example-polynomial-data``: polynomial problems solved by our GPU-HC <br />

# 2. Dependencies:
(1) CMake 3.14 or higher <br />
(2) MAGMA 2.6.1 or higher <br />
(3) CUDA 9.0 or higher <br />
(4) cuBlas <br />
(5) openBlas <br />
(6) pthread <br />
(7) if you want to solve a new polynomial problem not listed in the ``polynomial-problems`` folder, you will need Matlab 2019 or higher as well as [Macaulay2](http://www2.macaulay2.com/Macaulay2/Downloads/) and [Julia](https://julialang.org/downloads/) in order to generate symbolic Jacobian expressions and start solutions.

# 3. Constraints
The only restriction of GPU-HC is that the polynomial system of interest must have the number of variables lesser than 32. For those having more than 32 variables, you may need to reduce the size of the problem before using GPU-HC. <br />

# 4. How to build and compile the code
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


# 5. How to run the execution file after successful compilation
(1) enter the bin foler under build folder of the repo dorectory
```bash
cd bin
```
(2) run the code by specifically typing input arguments. See README.md in folders of either straight-line or parameter HC under ``GPU-HC`` folder.
```bash
./magmaHC-main <input-argument> <command>
```

# 6. How to use software tools and edit the code to solve a new polynomial problem
First, create a folder of your problem-name where all data can be placed. All scripts and materials mentioned below are given under the ``polynoamil-data-preparations`` folder. <br />
**STAIGHT-LINE HC** <br /> <br />
In this instruction, we take alea6 as an example. <br /><br />

**STEP 1: PREPARE YOUR TARGET SYSTEM**<br />
**(1) Create the problem:** create a matlab script of your target system, see ``/auto-gen-tools/polynomial-problems/sys_alea6.m`` for an example. <br />
**(2) Change symbols:** Edit your problem name, directory names, etc. in ``/auto-gen-tools/polynomial-problems/poly_coeffs_rep.m`` and run the script in Matlab. Specify the directory as the problem folder you created, and you will see three output files: <br />
``rep_problem.txt``: the taget system where coefficients are replaced by variables *a_i*. Same coefficients will be replaced by the same variables. <br />
``target_coeffs.txt``: the target coefficients are represented by complex numbers. Each row is a coefficient, while the first and second columns are real and imaginary numbers, respecitvely. <br />
``rep_coeffs.txt``: the target coefficients associated with the replaced variables *a_i*. This is provided for reference but might not be necessarily used. <br /><br />

**STEP 2: PREPARE YOUR START SYSTEM**<br />
**(1) Create your start system:** The start system has to have the same formulation as the target system. <br />
**(2) Find start solutions and start parameters:** Use Julia's monodromy solver to find the start solutions and start parameters of your start system. Please refer to ``/example-polynomial-data/julia_monodromySolver.jl`` for more details. <br />
**(3) Write to files:** After Julia's monodromy solver is done, write your start solutions into a file named ``julia-start-sols-raw``. For start coefficients, write to a file named ``start_coeffs.txt``. Examples are provided under ``/example-polynomial-data/alea6/``. <br />
**(4) Reformation:** Use ``/auto-gen-tools/reformateStartSolsFromJulia.m`` matlab script to reformat start solutions. Edit the input/output file directory and the output file named ``start-sols.txt`` will be generated. For start coefficients, manually reformat them as they appear in the example. <br /><br />

**STEP 3: GENERATE SYMBOLIC EXPRESSIONS OF JACOBIANS**<br />
**This is for the CPU-HC solver. Skip this step if you do not need a CPU solver.** <br />
**(1) Generate raw format from M2:** Use Macaulay2 to generate the C code of the symbolic Jacobian evaluations. An example is provided in ``/example-polynomial-data/alea6/M2-Jacobian-Evluation-Ccode``. <br />
**(2) Copy M2 generated C code to files:** Copy the generated C code to files (Please refer to the example files under ``/example-polynomial-data/alea6/``): <br /> 
``M2-HxHt-only-G-raw``: Only the G's of the HxHt symbolic evaluations. <br />
``M2-HxHt-only-y-raw``: Only the y's of the HxHt symbolic evaluations. <br />
``M2-HxH-only-G-raw``: Only the G's of the HxH symbolic evaluations. <br />
``M2-HxH-only-y-raw``: Only the y's of the HxH symbolic evaluations. <br />
**(3) Reformations:** Use ``/auto-gen-tools/reformatEval.m`` to generate files that change the Macaulay2 C code of ``M2-HxHt-only-G-raw`` and ``M2-HxH-only-G-raw`` into MAGMA data type so that we can use them directly in the c++ code. Be sure to remove the constants ``C0 = ...``, ``C1 = ...``, etc in ``M2-HxHt-only-G-raw`` and ``M2-HxH-only-G-raw``. Likewise, use ``/auto-gen-tools/reformat_y.m`` to generate files that change the Macaulay2 C code of ``M2-HxHt-only-y-raw`` and ``M2-HxH-only-y-raw``. <br /><br />

**STEP 4: GENERATE JACOBIAN EVALUATION INDICES FOR GPU-HC**<br />
**(1) Create your problem for index generator:** under ``/auto-gen-tools/idx-matrices-generator/``, create your problem matlab script. Add your problem in ``/auto-gen-tools/idx-matrices-generator/const_matrix_generator.m``. Make sure that the indices start from 0 and follow the format as they are in the provided scripts. <br />
**(2) Run index generator:** Run ``/auto-gen-tools/idx-matrices-generator/const_matrix_generator.m``. Change the file write directory if necessary. Two files, ``Hx.txt`` and ``Ht.txt`` will be created. <br /><br />

**STEP 5: CONSTRUCT YOUR HC SOLVER** <br />
**(1) Create a problem folder:** Under ``/GPU-HC/straight-line-HC/problems/``, create a problem folder containing ``Hx.txt``, ``Ht.txt``, your start coefficients ``start_coefs.txt``, your start solutions ``start_sols.txt``, and your target coefficients ``target_coefs.txt``. <br />
**(2) Define your problem parameters:** Go to ``/GPU-HC/straight-line-HC/define_params_and_read_files.cu`` and add your problem. <br />
**(3) Direct the solver to your problem:** Go to ``/GPU-HC/straight-line-HC/homotopy_continuation_solver.cu`` and put in your problem. Make sure to include your problem function in ``/GPU-HC/straight-line-HC/magmaHC/cpu-compute/cpu-compute.h`` for a CPU-HC solver, and ``/GPU-HC/straight-line-HC/magmaHC/gpu-kernels/magmaHC-kernels.h`` for GPU-HC solver. <br />
**(4) Create a CPU-HC solver:** Create Jacobian evaluation functions for your problem under ``/GPU-HC/straight-line-HC/magmaHC/cpu-compute/cpu-eval-HxHt.h`` and ``/GPU-HC/straight-line-HC/magmaHC/cpu-compute/cpu-eval-HxH.h``. Copy and paste the code from the files generated in STEP 3 (3). Create the solver under ``GPU-HC/straight-line-HC/magmaHC/cpu-compute/cpu-hc-solver/``. Change the values of ``C0``, ``C1``, etc if necessary (see STEP 3 (3)). <br />
**(5) Create a GPU-HC solver:** Create a CUDA kernel under ``/GPU-HC/straight-line-HC/magmaHC/gpu-kernels/`` and a Jacobian evaluation device function under ``/GPU-HC/straight-line-HC/magmaHC/gpu-idx-eval/``. <br /> For the CUDA kernel, please refer to the alea6 problem for reference. Make sure to put in the problem parameters as the templates for the kernels. For the jacobian evaluation device function, there are four functions: <br />
``eval_cdt_<problem_name>``: the evaluation of the coefficient homotopy ``(target_coefficient)*t + (start_coefficient)*(1-t)``. Also please refer to the device function in the alea6 problem. Change the loop number according to the number of coefficients and the number of your problem unknowns if necessary. <br />
``eval_Hx_<problem_name>``: the evaluation of the Jacobian Hx. <br />
``eval_Ht_<problem_name>``: the evaluation of the Jacobian Ht. <br />
``eval_H_<problem_name>``: the evaluation of the Homotopy H. <br />
Change the number of unknowns in each term of your problem for the above three device functions if necessary. <br /><br />

**STEP 6: INCLUDE NEWLY CREATED CODE FILES IN CMAKELISTS** <br />
Finally, add the code files you just created in the ``CMakeLists.txt`` under ``/GPU-HC/straight-line-HC/magmaHC/``. <br /><BR />

**STEP 7: NOW YOU SHOULD BE ABLE TO RUN YOUR HC SOLVER!** <br />

# 7. Important Update Notice (Oct. 9th 2022)
(1) In the published papers, start solutions are generated using Julia's homotopy continuation package. We found out that among all the generated start solutions, some of them are almost identical, creating redundant HC paths. Therefore, we have made a change on this where Julia's monodromy solver is used to generate start solutons. An example Julia script is given in ``/example-polynomial-data/``. <br />
(2) Straight-line HC CUDA kernel code have been optimized a bit, and the tolerance rate is strict in this repo. The timings for each problem could be a bit different from what we published in the papers. <br />
(3) We have discovered some numerical instable issues in the parameter HC. We are now improving this now and will update the repo as soon as possible. <br />

# 8. Known issues
When using ``#pragma unroll`` for some of the for-loops, an error message from nvcc comes out. Commenting out ``#pragma unroll`` solves the issue. <br />

