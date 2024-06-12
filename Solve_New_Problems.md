## How to solve a new (minimal) problem

Below we will use ``<problem_name>`` as the name of your minimal problem. <br />

- **Step 1. Create Your Polynomial System:** Formulate your _explicit_ polynomial system as a matlab script. An example can be found in ``auto-data-gen-tools/sys_*.m`` where ``p`` is the system parameters and ``x`` are the unknowns. Be sure to expand your polynomial equations. 
- **Step 2. Generate Data to be fed to GPU-HC:**
    - Generate Start Parameters and Solutions
	    - Create a Julia script running a Monodromy solver that finds the solutions of a start system. Examples can be found under ``start-sys-gen/<problem_name>_monodromySolve.jl``. Remember to write the start parameters and solutions to files when complete running Julia. <br />
	    - Run ``auto-data-gen-tools/reformatStartSolsFromJulia.m`` to reformulate the start solutions. Make sure to specify the file path and names. <br />
        - Reformulating start parameters is trivial. Manually reformulate it to the one as in ``GPU-HC/problems/<problem_name>/start_params.txt``. The first column is the real part while the second is the imaginary part. <br />
    - Create dH/dx, dH/dt, and H Evaluation Indices
        - Add your problem in the ``auto-data-gen-tools/autogen_gpuhc.m``. This needs the matlab script you prepared in Step 1. <br />
	    - Run ``auto-data-gen-tools/autogen_gpuhc.m``. You should expect several output files, but only the following files are necessary for GPU-HC: <br />
		    - ``dHdx_indx.txt``: Jacobian matrix $\frac{\partial H}{\partial x}$ evaluation indices. <br />
		    - ``dHdt_indx.txt``: Jacobian matrix $\frac{\partial H}{\partial t}$ evaluation indices. <br />
		    - ``p2c-<problem_name>.h``: Parameters to coefficient conversion C++ header code. <br />
			- ``dev-eval-indxing-<problem_name>.cuh``: device function CUH code. <br />
			- ``kernel_HC_Solver_<problem_name>.cu``: GPU kernel CUDA code. <br /> 
- **Step 3. Deploy your problem in GPU-HC:** 
	- Make a folder named <problem_name> under ``GPU-HC/problems/`` and put __(i)__ start system solutions and parameters, __(ii)__ ``dHdx_indx.txt`` and ``dHdt_indx.txt``, __(iii)__ your target parameters for which the formulation is identical to the start parameters, and __(iv)__ a ``gpuhc_settings.yaml`` file containing all necessary information of the problem.  <br />
		- In the YAML file, ``Num_Of_Tracks`` is the number of solutions which is given by the Julia monodromy solver. ``dHdx_Max_Terms``, ``dHdx_Max_Parts``, ``dHdt_Max_Terms``, ``dHdt_Max_Parts``, ``Max_Order_Of_T``, and ``Num_Of_Coeffs_From_Params`` can be found from the printed output after running ``auto-data-gen-tools/autogen_gpuhc.m`` MATLAB code. <br />
	- Move generated files in Step 2: <br />
		- Move ``p2c-<problem_name>.h`` under ``GPU-HC/magmaHC/PHC_Coeffs/``. <br />
		- Move ``dev-eval-indxing-<problem_name>.cuh`` under ``GPU-HC/magmaHC/gpu-idx-evals/``. <br />
		- Move ``kernel_HC_Solver_<problem_name>.cu`` under ``GPU-HC/magmaHC/gpu-kernels/``. Hard-coded the number of solutions by replacing XX in the line 
```cpp
const int num_of_tracks             = XX;
```
	- Add your minimal problem in ``GPU-HC/magmaHC/gpu-kernels/magmaHC-kernels.hpp`` and ``GPU-HC/magmaHC/GPU_HC_Solver.cpp`` (member function ``Solve_by_GPU_HC()``). <br />
	- Add your newly added files to the ``GPU-HC/magmaHC/CMakeLists.txt``. <br />
- **Step 4. Build and Compile:**
	- Now you should be able to build and compile and run GPU-HC. The converged solutions are written in a file under a created folder ``Output_Write_Files``. <br />
	- You can tune the HC hyper-parameters in the ``<problem_name>/gpuhc_settings.yaml``, _e.g._, ``GPUHC_Max_Steps``, ``GPUHC_Max_Correction_Steps``, and ``GPUHC_Num_Of_Steps_to_Increase_Delta_t`` to gain the fastest speed. <br />

<br />