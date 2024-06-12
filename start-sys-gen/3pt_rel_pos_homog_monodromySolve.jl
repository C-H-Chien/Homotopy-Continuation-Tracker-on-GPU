using Base
using HomotopyContinuation
using LinearAlgebra

@var x1 x2 x3 x4 x5 x6 x7 x8;
@var p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18;


f1  = x1*x6*x8*p1*p16 - x1*x6*p1*p13 + x2*x6*x8*p4*p16 - x2*x6*p4*p13 + x3*x6*x8*p7*p16 - x3*x6*p7*p13 - x4*p4*p16 - x5*p1*p16 + p7*p13;
f2 = -x1*x6*x7*p1*p16 + x1*x6*p1*p10 - x2*x6*x7*p4*p16 + x2*x6*p4*p10 - x3*x6*x7*p7*p16 + x3*x6*p7*p10 + x4*p1*p16 - x5*p4*p16 - p7*p10;
f3 = x1*x6*x8*p2*p17 - x1*x6*p2*p14 + x2*x6*x8*p5*p17 - x2*x6*p5*p14 + x3*x6*x8*p8*p17 - x3*x6*p8*p14 - x4*p5*p17 - x5*p2*p17 + p8*p14;
f4 = -x1*x6*x7*p2*p17 + x1*x6*p2*p11 - x2*x6*x7*p5*p17 + x2*x6*p5*p11 - x3*x6*x7*p8*p17 + x3*x6*p8*p11 + x4*p2*p17 - x5*p5*p17 - p8*p11;
f5 = x1*x6*x8*p3*p18 - x1*x6*p3*p15 + x2*x6*x8*p6*p18 - x2*x6*p6*p15 + x3*x6*x8*p9*p18 - x3*x6*p9*p15 - x4*p6*p18 - x5*p3*p18 + p9*p15;
f6 = -x1*x6*x7*p3*p18 + x1*x6*p3*p12 - x2*x6*x7*p6*p18 + x2*x6*p6*p12 - x3*x6*x7*p9*p18 + x3*x6*p9*p12 + x4*p3*p18 - x5*p6*p18 - p9*p12;
f7 = x1^2 + x2^2 + x3^2 - 1;
f8 = x4^2 + x5^2 - 1;

Eqs = [f1; f2; f3; f4; f5; f6; f7; f8];
variables_list = [x1; x2; x3; x4; x5; x6; x7; x8];
parameters_list = [p1; p2; p3; p4; p5; p6; p7; p8; p9; p10; p11; p12; p13; p14; p15; p16; p17; p18];
F = System(Eqs;variables=variables_list, parameters =parameters_list);

S = monodromy_solve(F)
start_solutions = solutions(S);
start_params = S.parameters;

#Try to solve with random Parameters
target_params = rand(18) + rand(18) * im

@time solve(F, start_solutions; start_parameters=start_params, target_parameters=target_params)


#io = open("/home/chchien/hcOutput", "w");
#using DelimitedFiles
#writedlm(io, start_solutions);
#close(io)
