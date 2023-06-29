using Base
using HomotopyContinuation
using LinearAlgebra

#> Rotation parametrization: Caylay
function cay2R(x,y,z)
    M = [1+x*x-(y*y+z*z)  2*(x*y-z)         2*(x*z+y);
         2*(x*y+z)        1+y^2-(x*x+z*z)   2*(y*z-x);
         2*(x*z-y)        2*(y*z+x)         1+z*z-(x*x+y*y)];
    return M;
end

#> Set rotations
@var r2[1:3] r3[1:3]
R2 = cay2R(r2[1], r2[2], r2[3]);
R3 = cay2R(r3[1], r3[2], r3[3]);

#> Set translations
@var t2[1:3];
@var t3[1:3];

#> Set parameters
@var x[1:3, 1:3, 1:2];
@var d[1:2, 1:3, 1:2];

#> Set unknowns
@var a[1:3, 1:3];
@var e[1:2, 1:3];
@var u[1:2, 1:3];

#> Point Equations
p = 1
pointEquations2 = a[p,2] * [x[p,2,:];1] - (R2 * (a[p,1] * [x[p,1,:];1]) + t2)
for p = 2:3
    Eq = a[p,2] * [x[p,2,:];1] - (R2 * (a[p,1] * [x[p,1,:];1]) + t2)
    for n = 1:3
        push!(pointEquations2, Eq[n]);
    end
end

p = 1
pointEquations3 = a[p,3] * [x[p,3,:];1] - (R3 * (a[p,1] * [x[p,1,:];1]) + t3)
for p = 2:3
    Eq = a[p,3] * [x[p,3,:];1] - (R3 * (a[p,1] * [x[p,1,:];1]) + t3)
    for n = 1:3
        push!(pointEquations3, Eq[n]);
    end
end
pointEquations = [pointEquations2; pointEquations3];

#> Tangent Equations
p = 1
tangentEquations2 = (e[p,2] * [x[p,2,:];1] + u[p,2] * [d[p,2,:];0]) - (R2 * (e[p,1] * [x[p,1,:];1] + u[p,1] * [d[p,1,:];0]));
for p = 2:2
    Eq = (e[p,2] * [x[p,2,:];1] + u[p,2] * [d[p,2,:];0]) - (R2 * (e[p,1] * [x[p,1,:];1] + u[p,1] * [d[p,1,:];0]))
    for n = 1:3
        push!(tangentEquations2, Eq[n]);
    end
end
p = 1
tangentEquations3 = (e[p,3] * [x[p,3,:];1] + u[p,3] * [d[p,3,:];0]) - (R3 * (e[p,1] * [x[p,1,:];1] + u[p,1] * [d[p,1,:];0]));
for p = 2:2
    Eq = (e[p,3] * [x[p,3,:];1] + u[p,3] * [d[p,3,:];0]) - (R3 * (e[p,1] * [x[p,1,:];1] + u[p,1] * [d[p,1,:];0]))
    for n = 1:3
        push!(tangentEquations3, Eq[n]);
    end
end
tangentEquations = [tangentEquations2; tangentEquations3];

#> Valid unknowns
a_unknowns = [a[2,1], a[3,1], a[1,2], a[2,2], a[3,2], a[1,3], a[2,3], a[3,3]];
e_unknowns = [e[1,2], e[1,3], e[2,2], e[2,3]];

#> All Equations, 30 equations in total
Eqs = [pointEquations; tangentEquations];
variables_list = collect(Iterators.flatten([transpose(a_unknowns), transpose(e_unknowns), transpose(u), t2, t3, r2, r3]));
parameters_list = collect(Iterators.flatten([permutedims(x,[3,2,1]), permutedims(d,[3,2,1]), a[1,1], e[1,1], e[2,1]]))
F = System(Eqs;variables=variables_list, parameters =parameters_list);

S = monodromy_solve(F)
start_solutions = solutions(S);
start_params = S.parameters;

#Try to solve with random Parameters
#target_params = rand(33) + rand(33) * im

#@time for i = 1:10
#solve(F, start_solutions; start_parameters=start_params, target_parameters=target_params)
#end

#io = open("/path/to/folder/", "w");
#using DelimitedFiles
#writedlm(io, start_solutions);
#close(io)