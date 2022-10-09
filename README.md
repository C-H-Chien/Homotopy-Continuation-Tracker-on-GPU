GPU-HC: Path Traking for Homotopy Continuation in GPU for Benchmark Polynomial Problems

# 1. Contents
This repository primarily contains three folders: <br />
``GPU-HC``: GPU implementation of both straight-line and parameter homotopy continuation algorithms <br />
``auto-gen-tools``: software tools used to automatically generate necessary data and code for GPU-HC to solve new problems <br />
``polynomial-problems``: polynomial problems solved by our GPU-HC
``polynomial-data-preparations``: example data files used for the software tools

# 2. Dependencies:
(1) CMake 3.14 or higher <br />
(2) MAGMA 2.6.1 or higher <br />
(3) CUDA 9.0 or higher <br />
(4) cuBlas <br />
(5) openBlas <br />
(6) pthread <br />
(7) if you want to solve a new polynomial problem not listed in the ``polynomial-problems`` folder, you will need Matlab 2019 or higher as well as [Macaulay2](http://www2.macaulay2.com/Macaulay2/Downloads/) and [Julia](https://julialang.org/downloads/) in order to generate symbolic Jacobian expressions and start solutions.


# 3. How to build and compile the code
(1) After cloning the repo, cd to the repo folder and create a 'build' directory and enter it
```bash
mkdir build && cd build
```
(2) make sure to change routes in CMakeFiles according to your local machine. See README.md of either straight-line or parameter HC. <br />
(3) create a make file and compile the entire code
```bash
cmake .. && make -j
```

# 4. How to run the execution file after successful compilation
(1) enter the bin foler under build folder of the repo dorectory
```bash
cd bin
```
(2) run the code by specifically typing input arguments. See README.md in folders of either straight-line or parameter HC under ``GPU-HC`` folder.
```bash
./magmaHC-main <input-argument> <command>
```

# 5. How to use software tools and edit the code to solve a new polynomial problem
First, create a folder of your problem-name where all data can be placed. All scripts and materials mentioned below are given under the ``polynoamil-data-preparations`` folder. <br />
**STAIGHT-LINE HC** <br /> <br />
In this instruction, we take alea6 as an example. <br /><br />

**STEP 1: PREPARE YOUR TARGET SYSTEM**<br />
(1) Once you have the formulation of the polynomial problem you need to solve (i.e., target system), create a matlab script of the system. See ``sys_alea6.m`` under ``auto-gen-tools/polynomial-problems/``. <br />
(2) Edit the problem name, directory names, etc. in ``poly_coeffs_rep.m`` and run the script in Matlab. Specify the directory as the problem folder you created, and you will see three output files: <br />
``rep_problem.txt``: the taget system where coefficients are replaced by variables *a_i*. Same coefficients will be replaced by the same variables. <br />
``target_coeffs.txt``: the target coefficients are represented by complex numbers. Each row is a coefficient, while the first and second columns are real and imaginary numbers, respecitvely. <br />
``rep_coeffs.txt``: the target coefficients associated with the replaced variables *a_i*. This is provided for reference but might not be necessary used. <br /><br />

**STEP 2: PREPARE YOUR START SYSTEM**<br />
**(1) Create your start system:** The start system has to have the same formulation as the target system. <br />
**(2) Find start solutions and start parameters:** Use Julia's monodromy solver to find the start solutions and start parameters of your start system. Please refer to ``/example-polynomial-data/julia_monodromySolver.jl`` for more details. <br />
**(3) Write to files:** After Julia's monodromy solver is done, write your start solutions into a file named ``julia-start-sols-raw``. For start coefficients, write to a file named ``start_coeffs.txt``. Examples are provided under ``/example-polynomial-data/alea6/``.
**(4) Reformation:** Use ``/auto-gen-tools/reformateStartSolsFromJulia.m`` matlab script to reformat start solutions. Edit the input/output file directory and the output file named ``start-sols.txt`` will be generated. For start coefficients, manually reformat them as they appear in the example. <br /><br />

**STEP 3: GENERATE SYMBOLIC EXPRESSIONS OF JACOBIANS**<br />
**(1) Generate raw format from M2:** Use Macaulay2 to generate the C code of the symbolic Jacobian evaluations. An example is provided in ``/example-polynomial-data/alea6/M2-Jacobian-Evluation-Ccode``.
**(2) Copy M2 generated C code to files:** Copy the generated C code to files. TODO 

# 6. Important Update Notice (Oct. 9th 2022)
(1) In the published papers, start solutions are generated using Julia's homotopy continuation package. We found out that among all the generated start solutions, some of them are almost identical, creating redundant HC paths. Therefore, we have made a change on this where Julia's monodromy solver is used to generate start solutons. An example Julia script is given in ``/example-polynomial-data/``. <br />
(2) Straight-line HC CUDA kernel code have been optimized a bit. The timings for each problem could be a bit different from what we published in the papers.
(3) We have discovered some numerical instable issues in the parameter HC. We are now improving this and will update the repo as soon as possible.

# 7. Reference
Please cite the following papers if you use this code: <br />
Straight-line HC: <br />
``Chien, Chiang-Heng, Hongyi Fan, Ahmad Abdelfattah, Elias Tsigaridas, Stanimire Tomov, and Benjamin Kimia. "GPU-Based Homotopy Continuation for Minimal Problems in Computer Vision." In Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition, pp. 15765-15776. 2022.`` <br />
Parameter HC: <br />
``Chien, Chiang-Heng, Hongyi Fan, Ahmad Abdelfattah, Elias Tsigaridas, Stanimire Tomov, and Benjamin Kimia. "Parallel Path Tracking for Homotopy Continuation using GPU." In Proceedings of the International Symposium on Symbolic and Algebraic Computation. 2022.``

