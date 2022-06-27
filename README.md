GPU-HC: Path Traking for Homotopy Continuation in GPU for Benchmark Polynomial Problems

# 1. Contents
This repository primarily contains two folders: <br />
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
(7) if you want to solve a new polynomial problem not listed in the ``polynomial-problems`` folder, you will need Matlab 2019 or higher as well as [MacaulAy2](http://www2.macaulay2.com/Macaulay2/Downloads/) and [Julia](https://julialang.org/downloads/) in order to generate symbolic Jacobian expressions and start solutions.

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
**STEP 1: PREPARE TARGET SYSTEM**<br />
(1) Once you have the formulation of the polynomial problem you need to solve (i.e., target system), create a matlab script of the system. See ``sys_alea6.m`` under ``auto-gen-tools/polynomial-problems/``.
(2) Edit the problem name, directory names, etc. in ``poly_coeffs_rep.m`` and run the script in Matlab. Specify the directory as the problem folder you created, and you will see four output files: <br />
``rep_problem.txt``: the taget system where coefficients are replaced by variables *a_i*. Same coefficients will be replaced by the same variables. <br />
``target_coeffs.txt``: the target coefficients are represented by complex numbers. Each row is a coefficient, while the first and second columns are real and imaginary numbers, respecitvely. <br />
``rep_coeffs.txt``: the target coefficients associated with the replaced variables *a_i*. This is provided for reference but might not be necessary used. <br /><br />
**STEP 2: PREPARE START SYSTEM**<br />
(1) Because the start system has to have the same formate as the target system, we can use the target problem representation in ``rep_problem.txt`` to create a start system and generate start solutions. To do so, edit variables, start coefficients, start system formulation from the file ``jl_start_sols.jl`` under the ``polynomial-data-preparations``. <br />
(2) Use julia to run ``jl_start_sols.jl``. In the end, the start solutions will be put in the ``julia-start-sols-raw`` file which Julia writes to.<br />
(3) To reformate the start solutions Julia created, use ``reformateStartSolsFromJulia.m`` matlab script. Edit the input/output file directory and a file named ``start-sols.txtx`` will be generated.<br /><br />
**STEP 3: GENERATE SYMBOLIC EXPRESSIONS OF JACOBIANS**<br />
This part will be updated soon.

# 6. Reference
The straigh-line HC was proposed in the paper <br />
``Chien, Chiang-Heng, Hongyi Fan, Ahmad Abdelfattah, Elias Tsigaridas, Stanimire Tomov, and Benjamin Kimia. "GPU-Based Homotopy Continuation for Minimal Problems in Computer Vision." In Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition, pp. 15765-15776. 2022.`` <br />
and the parameter HC was proposed in the paper <br />
``Chien, Chiang-Heng, Hongyi Fan, Ahmad Abdelfattah, Elias Tsigaridas, Stanimire Tomov, and Benjamin Kimia. "Parallel Path Tracking for Homotopy Continuation using GPU." In Proceedings of the International Symposium on Symbolic and Algebraic Computation. 2022.``

