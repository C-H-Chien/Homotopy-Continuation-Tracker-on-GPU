using Base
using HomotopyContinuation
using LinearAlgebra

@var x1 x2 x3 x4 x5
@var p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20 p21 p22 p23 p24 p25 p26 p27 p28 p29 p30 p31 p32 p33 p34 p35 p36 p37 p38 p39 p40 p41 p42 p43 p44 p45 p46 p47 p48

f1 = x1^2*p13 + x1*x2*p14 + x1*x2*p16 + x1*x3*p15 + x1*x3*p19 + x1*x4*p2*p14 + x1*x4*p2*p16 + x1*x4*p3*p15 + x1*x4*p3*p19 + x1*x5*p2*p15 + x1*x5*p2*p19 - x1*x5*p3*p14 - x1*x5*p3*p16 + 2*x1*p1*p13 + x2^2*p17 + x2*x3*p18 + x2*x3*p20 + 2*x2*x4*p2*p17 + x2*x4*p3*p18 + x2*x4*p3*p20 + x2*x5*p2*p18 + x2*x5*p2*p20 - 2*x2*x5*p3*p17 + x2*p1*p14 + x2*p1*p16 + x3^2*p21 + x3*x4*p2*p18 + x3*x4*p2*p20 + 2*x3*x4*p3*p21 + 2*x3*x5*p2*p21 - x3*x5*p3*p18 - x3*x5*p3*p20 + x3*p1*p15 + x3*p1*p19 + x4^2*p2^2*p17 + x4^2*p2*p3*p18 + x4^2*p2*p3*p20 + x4^2*p3^2*p21 + x4*x5*p2^2*p18 + x4*x5*p2^2*p20 - 2*x4*x5*p2*p3*p17 + 2*x4*x5*p2*p3*p21 - x4*x5*p3^2*p18 - x4*x5*p3^2*p20 + x4*p1*p2*p14 + x4*p1*p2*p16 + x4*p1*p3*p15 + x4*p1*p3*p19 + x5^2*p2^2*p21 - x5^2*p2*p3*p18 - x5^2*p2*p3*p20 + x5^2*p3^2*p17 + x5*p1*p2*p15 + x5*p1*p2*p19 - x5*p1*p3*p14 - x5*p1*p3*p16 + p1^2*p13;
f2 = x1^2*p22 + x1*x2*p23 + x1*x2*p25 + x1*x3*p24 + x1*x3*p28 + x1*x4*p5*p23 + x1*x4*p5*p25 + x1*x4*p6*p24 + x1*x4*p6*p28 + x1*x5*p5*p24 + x1*x5*p5*p28 - x1*x5*p6*p23 - x1*x5*p6*p25 + 2*x1*p4*p22 + x2^2*p26 + x2*x3*p27 + x2*x3*p29 + 2*x2*x4*p5*p26 + x2*x4*p6*p27 + x2*x4*p6*p29 + x2*x5*p5*p27 + x2*x5*p5*p29 - 2*x2*x5*p6*p26 + x2*p4*p23 + x2*p4*p25 + x3^2*p30 + x3*x4*p5*p27 + x3*x4*p5*p29 + 2*x3*x4*p6*p30 + 2*x3*x5*p5*p30 - x3*x5*p6*p27 - x3*x5*p6*p29 + x3*p4*p24 + x3*p4*p28 + x4^2*p5^2*p26 + x4^2*p5*p6*p27 + x4^2*p5*p6*p29 + x4^2*p6^2*p30 + x4*x5*p5^2*p27 + x4*x5*p5^2*p29 - 2*x4*x5*p5*p6*p26 + 2*x4*x5*p5*p6*p30 - x4*x5*p6^2*p27 - x4*x5*p6^2*p29 + x4*p4*p5*p23 + x4*p4*p5*p25 + x4*p4*p6*p24 + x4*p4*p6*p28 + x5^2*p5^2*p30 - x5^2*p5*p6*p27 - x5^2*p5*p6*p29 + x5^2*p6^2*p26 + x5*p4*p5*p24 + x5*p4*p5*p28 - x5*p4*p6*p23 - x5*p4*p6*p25 + p4^2*p22;
f3 = x1^2*p31 + x1*x2*p32 + x1*x2*p34 + x1*x3*p33 + x1*x3*p37 + x1*x4*p8*p32 + x1*x4*p8*p34 + x1*x4*p9*p33 + x1*x4*p9*p37 + x1*x5*p8*p33 + x1*x5*p8*p37 - x1*x5*p9*p32 - x1*x5*p9*p34 + 2*x1*p7*p31 + x2^2*p35 + x2*x3*p36 + x2*x3*p38 + 2*x2*x4*p8*p35 + x2*x4*p9*p36 + x2*x4*p9*p38 + x2*x5*p8*p36 + x2*x5*p8*p38 - 2*x2*x5*p9*p35 + x2*p7*p32 + x2*p7*p34 + x3^2*p39 + x3*x4*p8*p36 + x3*x4*p8*p38 + 2*x3*x4*p9*p39 + 2*x3*x5*p8*p39 - x3*x5*p9*p36 - x3*x5*p9*p38 + x3*p7*p33 + x3*p7*p37 + x4^2*p8^2*p35 + x4^2*p8*p9*p36 + x4^2*p8*p9*p38 + x4^2*p9^2*p39 + x4*x5*p8^2*p36 + x4*x5*p8^2*p38 - 2*x4*x5*p8*p9*p35 + 2*x4*x5*p8*p9*p39 - x4*x5*p9^2*p36 - x4*x5*p9^2*p38 + x4*p7*p8*p32 + x4*p7*p8*p34 + x4*p7*p9*p33 + x4*p7*p9*p37 + x5^2*p8^2*p39 - x5^2*p8*p9*p36 - x5^2*p8*p9*p38 + x5^2*p9^2*p35 + x5*p7*p8*p33 + x5*p7*p8*p37 - x5*p7*p9*p32 - x5*p7*p9*p34 + p7^2*p31;
f4 = x1^2*p40 + x1*x2*p41 + x1*x2*p43 + x1*x3*p42 + x1*x3*p46 + x1*x4*p11*p41 + x1*x4*p11*p43 + x1*x4*p12*p42 + x1*x4*p12*p46 + x1*x5*p11*p42 + x1*x5*p11*p46 - x1*x5*p12*p41 - x1*x5*p12*p43 + 2*x1*p10*p40 + x2^2*p44 + x2*x3*p45 + x2*x3*p47 + 2*x2*x4*p11*p44 + x2*x4*p12*p45 + x2*x4*p12*p47 + x2*x5*p11*p45 + x2*x5*p11*p47 - 2*x2*x5*p12*p44 + x2*p10*p41 + x2*p10*p43 + x3^2*p48 + x3*x4*p11*p45 + x3*x4*p11*p47 + 2*x3*x4*p12*p48 + 2*x3*x5*p11*p48 - x3*x5*p12*p45 - x3*x5*p12*p47 + x3*p10*p42 + x3*p10*p46 + x4^2*p11^2*p44 + x4^2*p11*p12*p45 + x4^2*p11*p12*p47 + x4^2*p12^2*p48 + x4*x5*p11^2*p45 + x4*x5*p11^2*p47 - 2*x4*x5*p11*p12*p44 + 2*x4*x5*p11*p12*p48 - x4*x5*p12^2*p45 - x4*x5*p12^2*p47 + x4*p10*p11*p41 + x4*p10*p11*p43 + x4*p10*p12*p42 + x4*p10*p12*p46 + x5^2*p11^2*p48 - x5^2*p11*p12*p45 - x5^2*p11*p12*p47 + x5^2*p12^2*p44 + x5*p10*p11*p42 + x5*p10*p11*p46 - x5*p10*p12*p41 - x5*p10*p12*p43 + p10^2*p40;
f5 = x4^2 + x5^2 - 1;
    
F = System([f1, f2, f3, f4, f5]; variables=[x1, x2, x3, x4, x5], parameters = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48]);

S = monodromy_solve(F)
start_solutions = solutions(S);
start_params = S.parameters;

target_params = rand(48) + rand(48) * im;
@time solve(F, start_solutions; start_parameters=start_params, target_parameters=target_params)


#> write target solutions to a file
#io = open("/users/cchien3/data/cchien3/hcOutput", "w");
#using DelimitedFiles
#writedlm(io, solutions(R));
#close(io)
