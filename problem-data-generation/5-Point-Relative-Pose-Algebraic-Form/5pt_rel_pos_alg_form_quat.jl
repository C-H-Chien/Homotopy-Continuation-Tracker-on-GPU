#> 5-point relative pose problem (Algebraic form, quaternion parametrization)

using Base
using HomotopyContinuation
using LinearAlgebra

#> Rotation parametrization: Quaternion
function quat2R(x,y,z,w)
    M = [1-2*y*y-2*z*z    2*x*y-2*w*z       2*x*z+2*w*y;
         2*x*y+2*w*z      1-2*x*x-2*z*z     2*(y*z-w*x);
         2*(x*z-w*y)      2*(y*z+w*x)       1-2*x*x-2*y*y];
    return M;
end

#> Set rotations
@var r[1:4]
R = quat2R(r[1], r[2], r[3], r[4]);


#> Set translations (make the third dimension as 1 for normalization)
@var t[1:2];

#> Set parameters (pairs of point correspondences on two images)
@var x[1:5, 1:2, 1:2];

#> Essential Matrix
t_x = [0     -1    t[2]; 
       1     0    -t[1]; 
       -t[2] t[1] 0     ];
E = t_x * R;

#> Point Equations
#> 5 pairs of correspondences
Eq1 = [x[1,2,:];1]' * E * [x[1,1,:];1];
Eq2 = [x[2,2,:];1]' * E * [x[2,1,:];1];
Eq3 = [x[3,2,:];1]' * E * [x[3,1,:];1];
Eq4 = [x[4,2,:];1]' * E * [x[4,1,:];1];
Eq5 = [x[5,2,:];1]' * E * [x[5,1,:];1];

#> Quaternion Normalization Constraint Eqs
quat_nom_R = r[1]^2 + r[2]^2 + r[3]^2 + r[4]^2 - 1 

#> All Equations, 16 equations in total
Eqs = [Eq1; Eq2; Eq3; Eq4; Eq5; quat_nom_R];
variables_list = collect(Iterators.flatten([t, r]));
parameters_list = collect(Iterators.flatten([permutedims(x,[3,2,1])]))
F = System(Eqs;variables=variables_list, parameters =parameters_list);

S = monodromy_solve(F)
start_solutions = solutions(S);
start_params = S.parameters;

#Try to solve with random Parameters
#target_params = rand(33) + rand(33) * im

#@time for i = 1:10
#solve(F, start_solutions; start_parameters=start_params, target_parameters=target_params)
#end

#io = open("/home/chchien/hcOutput", "w");
#using DelimitedFiles
#writedlm(io, start_solutions);
#close(io)