using Base
using HomotopyContinuation
using LinearAlgebra

@var x1 x2 x3 x4 x5 x6 x7 x8 x9 x10
@var p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20 p21 p22 p23 p24 p25 p26 p27 p28 p29 p30 p31 p32 p33 p34

f1 = 2*x1 + 2*x5*x7*p1 - 2*x5*x7*p4 + 2*x5*x8*p1 - 2*x5*x8*p7 + 2*x5*x9*p1 - 2*x5*x9*p10 + 2*x5*x10*p1 - 2*x5*x10*p13;
f2 = 2*x2 + 2*x5*x7*p2 - 2*x5*x7*p5 + 2*x5*x8*p2 - 2*x5*x8*p8 + 2*x5*x9*p2 - 2*x5*x9*p11 + 2*x5*x10*p2 - 2*x5*x10*p14;
f3 = 2*x3 - 2*x6*x7*p16 + 2*x6*x7*p19 - 2*x6*x8*p16 + 2*x6*x8*p22 - 2*x6*x9*p16 + 2*x6*x9*p25 - 2*x6*x10*p16 + 2*x6*x10*p28;
f4 = 2*x4 - 2*x6*x7*p17 + 2*x6*x7*p20 - 2*x6*x8*p17 + 2*x6*x8*p23 - 2*x6*x9*p17 + 2*x6*x9*p26 - 2*x6*x10*p17 + 2*x6*x10*p29;
f5 = 2*x1*x7*p1 - 2*x1*x7*p4 + 2*x1*x8*p1 - 2*x1*x8*p7 + 2*x1*x9*p1 - 2*x1*x9*p10 + 2*x1*x10*p1 - 2*x1*x10*p13 + 2*x2*x7*p2 - 2*x2*x7*p5 + 2*x2*x8*p2 - 2*x2*x8*p8 + 2*x2*x9*p2 - 2*x2*x9*p11 + 2*x2*x10*p2 - 2*x2*x10*p14 - 2*x7*p1*p31 - 2*x7*p2*p32 - 2*x7*p3 + 2*x7*p4*p31 + 2*x7*p5*p32 + 2*x7*p6 - 2*x8*p1*p31 - 2*x8*p2*p32 - 2*x8*p3 + 2*x8*p7*p31 + 2*x8*p8*p32 + 2*x8*p9 - 2*x9*p1*p31 - 2*x9*p2*p32 - 2*x9*p3 + 2*x9*p10*p31 + 2*x9*p11*p32 + 2*x9*p12 - 2*x10*p1*p31 - 2*x10*p2*p32 - 2*x10*p3 + 2*x10*p13*p31 + 2*x10*p14*p32 + 2*x10*p15;
f6 = -2*x3*x7*p16 + 2*x3*x7*p19 - 2*x3*x8*p16 + 2*x3*x8*p22 - 2*x3*x9*p16 + 2*x3*x9*p25 - 2*x3*x10*p16 + 2*x3*x10*p28 - 2*x4*x7*p17 + 2*x4*x7*p20 - 2*x4*x8*p17 + 2*x4*x8*p23 - 2*x4*x9*p17 + 2*x4*x9*p26 - 2*x4*x10*p17 + 2*x4*x10*p29 + 2*x7*p16*p33 + 2*x7*p17*p34 + 2*x7*p18 - 2*x7*p19*p33 - 2*x7*p20*p34 - 2*x7*p21 + 2*x8*p16*p33 + 2*x8*p17*p34 + 2*x8*p18 - 2*x8*p22*p33 - 2*x8*p23*p34 - 2*x8*p24 + 2*x9*p16*p33 + 2*x9*p17*p34 + 2*x9*p18 - 2*x9*p25*p33 - 2*x9*p26*p34 - 2*x9*p27 + 2*x10*p16*p33 + 2*x10*p17*p34 + 2*x10*p18 - 2*x10*p28*p33 - 2*x10*p29*p34 - 2*x10*p30;
f7 = 2*x1*x5*p1 - 2*x1*x5*p4 + 2*x2*x5*p2 - 2*x2*x5*p5 - 2*x3*x6*p16 + 2*x3*x6*p19 - 2*x4*x6*p17 + 2*x4*x6*p20 - 2*x5*p1*p31 - 2*x5*p2*p32 - 2*x5*p3 + 2*x5*p4*p31 + 2*x5*p5*p32 + 2*x5*p6 + 2*x6*p16*p33 + 2*x6*p17*p34 + 2*x6*p18 - 2*x6*p19*p33 - 2*x6*p20*p34 - 2*x6*p21 + p1^2 + p2^2 + p3^2 - p4^2 - p5^2 - p6^2 - p16^2 - p17^2 - p18^2 + p19^2 + p20^2 + p21^2;
f8 = 2*x1*x5*p1 - 2*x1*x5*p7 + 2*x2*x5*p2 - 2*x2*x5*p8 - 2*x3*x6*p16 + 2*x3*x6*p22 - 2*x4*x6*p17 + 2*x4*x6*p23 - 2*x5*p1*p31 - 2*x5*p2*p32 - 2*x5*p3 + 2*x5*p7*p31 + 2*x5*p8*p32 + 2*x5*p9 + 2*x6*p16*p33 + 2*x6*p17*p34 + 2*x6*p18 - 2*x6*p22*p33 - 2*x6*p23*p34 - 2*x6*p24 + p1^2 + p2^2 + p3^2 - p7^2 - p8^2 - p9^2 - p16^2 - p17^2 - p18^2 + p22^2 + p23^2 + p24^2;
f9 = 2*x1*x5*p1 - 2*x1*x5*p10 + 2*x2*x5*p2 - 2*x2*x5*p11 - 2*x3*x6*p16 + 2*x3*x6*p25 - 2*x4*x6*p17 + 2*x4*x6*p26 - 2*x5*p1*p31 - 2*x5*p2*p32 - 2*x5*p3 + 2*x5*p10*p31 + 2*x5*p11*p32 + 2*x5*p12 + 2*x6*p16*p33 + 2*x6*p17*p34 + 2*x6*p18 - 2*x6*p25*p33 - 2*x6*p26*p34 - 2*x6*p27 + p1^2 + p2^2 + p3^2 - p10^2 - p11^2 - p12^2 - p16^2 - p17^2 - p18^2 + p25^2 + p26^2 + p27^2;
f10 = 2*x1*x5*p1 - 2*x1*x5*p13 + 2*x2*x5*p2 - 2*x2*x5*p14 - 2*x3*x6*p16 + 2*x3*x6*p28 - 2*x4*x6*p17 + 2*x4*x6*p29 - 2*x5*p1*p31 - 2*x5*p2*p32 - 2*x5*p3 + 2*x5*p13*p31 + 2*x5*p14*p32 + 2*x5*p15 + 2*x6*p16*p33 + 2*x6*p17*p34 + 2*x6*p18 - 2*x6*p28*p33 - 2*x6*p29*p34 - 2*x6*p30 + p1^2 + p2^2 + p3^2 - p13^2 - p14^2 - p15^2 - p16^2 - p17^2 - p18^2 + p28^2 + p29^2 + p30^2;

F = System([f1, f2, f3, f4, f5, f6, f7, f8, f9, f10]; variables=[x1, x2, x3, x4, x5, x6, x7, x8, x9, x10], parameters = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34]);

S = monodromy_solve(F)
start_solutions = solutions(S);
start_params = S.parameters;

target_params = rand(34) + rand(34) * im;
@time solve(F, start_solutions; start_parameters=start_params, target_parameters=target_params)


#> write target solutions to a file
#io = open("/home/chchien/hcOutput", "w");
#using DelimitedFiles
#writedlm(io, solutions(R));
#close(io)
