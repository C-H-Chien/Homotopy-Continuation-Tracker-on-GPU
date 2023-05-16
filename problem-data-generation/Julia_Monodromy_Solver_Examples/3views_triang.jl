using Base
using HomotopyContinuation
using LinearAlgebra

@var x1 x2 x3 x4 x5 x6 x7 x8;
@var p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20 p21 p22 p23 p24;


f1  = 2*x1 + p1*x3*x7 + p2*x4*x7 + p3*x7 - 2*p19*1;
f2  = 2*x2 + p4*x3*x7 + p5*x4*x7 + p6*x7 - 2*p20*1;
f3  = p1*x1*x7 + p4*x2*x7 + 2*x3 + p10*x5*x8 + p11*x6*x8 + p7*x7 + p12*x8 - 2*p21*1;
f4  = p2*x1*x7 + p5*x2*x7 + 2*x4 + p13*x5*x8 + p14*x6*x8 + p8*x7 + p15*x8 - 2*p22*1;
f5  = p10*x3*x8 + p13*x4*x8 + 2*x5 + p16*x8 - 2*p23*1;
f6  = p11*x3*x8 + p14*x4*x8 + 2*x6 + p17*x8 - 2*p24*1;
f7  = p1*x1*x3 + p2*x1*x4 + p3*x1 + p4*x2*x3 + p5*x2*x4 + p6*x2 + p7*x3 + p8*x4 + p9*1;
f8  = p10*x3*x5 + p11*x3*x6 + p12*x3 + p13*x4*x5 + p14*x4*x6 + p15*x4 + p16*x5 + p17*x6 + p18*1;

Eqs = [f1; f2; f3; f4; f5; f6; f7; f8];
variables_list = [x1; x2; x3; x4; x5; x6; x7; x8];
parameters_list = [p1; p2; p3; p4; p5; p6; p7; p8; p9; p10; p11; p12; p13; p14; p15; p16; p17; p18; p19; p20; p21; p22; p23; p24];
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