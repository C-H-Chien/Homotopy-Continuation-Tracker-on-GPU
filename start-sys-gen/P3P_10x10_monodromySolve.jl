
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

#> Set rotation
@var r[1:4]
R = quat2R(r[1], r[2], r[3], r[4]);

#> Set translation
@var t[1:3];

#> Set depth
@var d[1:3]

#> Set parameters (2D image points)
@var x[1:3, 1:2];
@var X[1:3, 1:3];

Eq1 = d[1] * R * [x[1,:];1] + t - X[1,:];
Eq2 = d[2] * R * [x[2,:];1] + t - X[2,:];
Eq3 = d[3] * R * [x[3,:];1] + t - X[3,:];

#> Quaternion Normalization Constraint Eqs
quat_nom_R = r[1]^2 + r[2]^2 + r[3]^2 + r[4]^2 - 1 

#> All Equations, 16 equations in total
Eqs = [Eq1; Eq2; Eq3; quat_nom_R];
variables_list = collect(Iterators.flatten([r, t, d]));
parameters_list = collect(Iterators.flatten([permutedims(x,[2,1]), permutedims(X,[2,1])]))
F = System(Eqs;variables=variables_list, parameters =parameters_list);

S = monodromy_solve(F)
start_solutions = solutions(S);
start_params = S.parameters;

