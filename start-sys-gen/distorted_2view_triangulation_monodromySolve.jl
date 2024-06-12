using Base
using HomotopyContinuation
using LinearAlgebra

@var x1 x2 x3 x4 x5
@var p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15

f1 = p13 + x1*p7 + x2*p10 + x3*p11 + x4*p12 + x1^2*p13*p14 + x2^2*p13*p14 + x3^2*p13*p15 + x4^2*p13*p15 + x1*x3*p5 + x1*x4*p6 + x2*x3*p8 + x2*x4*p9 + x1*x3^2*p7*p15 + x1*x4^2*p7*p15 + x1^2*x3*p11*p14 + x2*x3^2*p10*p15 + x2^2*x3*p11*p14 + x2*x4^2*p10*p15 + x1^2*x4*p12*p14 + x2^2*x4*p12*p14 + x1^2*x3^2*p13*p14*p15 + x1^2*x4^2*p13*p14*p15 + x2^2*x3^2*p13*p14*p15 + x2^2*x4^2*p13*p14*p15;
f2 = 2*x1 - 2*p1 + 2*x5*p7 + 2*x3*x5*p5 + 2*x4*x5*p6 + 4*x1*x5*p13*p14 + 2*x3^2*x5*p7*p15 + 2*x4^2*x5*p7*p15 + 4*x1*x3*x5*p11*p14 + 4*x1*x4*x5*p12*p14 + 4*x1*x3^2*x5*p13*p14*p15 + 4*x1*x4^2*x5*p13*p14*p15;
f3 = 2*x2 - 2*p2 + 2*x5*p10 + 2*x3*x5*p8 + 2*x4*x5*p9 + 4*x2*x5*p13*p14 + 2*x3^2*x5*p10*p15 + 2*x4^2*x5*p10*p15 + 4*x2*x3*x5*p11*p14 + 4*x2*x4*x5*p12*p14 + 4*x2*x3^2*x5*p13*p14*p15 + 4*x2*x4^2*x5*p13*p14*p15;
f4 = 2*x3 - 2*p3 + 2*x5*p7 + 2*x1*x5*p5 + 2*x2*x5*p6 + 4*x3*x5*p13*p15 + 2*x1^2*x5*p7*p14 + 2*x2^2*x5*p7*p14 + 4*x1*x3*x5*p11*p15 + 4*x2*x3*x5*p12*p15 + 4*x1^2*x3*x5*p13*p14*p15 + 4*x2^2*x3*x5*p13*p14*p15;
f5 = 2*x4 - 2*p4 + 2*x5*p10 + 2*x1*x5*p8 + 2*x2*x5*p9 + 4*x4*x5*p13*p15 + 2*x1^2*x5*p10*p14 + 2*x2^2*x5*p10*p14 + 4*x1*x4*x5*p11*p15 + 4*x2*x4*x5*p12*p15 + 4*x1^2*x4*x5*p13*p14*p15 + 4*x2^2*x4*x5*p13*p14*p15;

F = System([f1, f2, f3, f4, f5]; variables=[x1, x2, x3, x4, x5], parameters = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15]);

S = monodromy_solve(F)
start_solutions = solutions(S);
start_params = S.parameters;

target_params = rand(15) + rand(15) * im;
@time solve(F, start_solutions; start_parameters=start_params, target_parameters=target_params)


#> write target solutions to a file
#io = open("/users/cchien3/data/cchien3/hcOutput", "w");
#using DelimitedFiles
#writedlm(io, solutions(R));
#close(io)
