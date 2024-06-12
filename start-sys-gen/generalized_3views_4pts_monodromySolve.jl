####################################################################################
#> Code Description
#    This code generates solutions from a set of polynomials using Julia's 
#    monodromy solver. At the end of the code, it is optional for the user to 
#    write the start solutions to a file in a user-defined directory.
#
#> How to run this code
#    1. In terminal, launch Julia
#    2. use the command $ include("/path/to/this/code/six_lines_16.jl");
#
#> (c) LEMS, Brown University
#> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
#> Last modified:  Nov. 22nd, 2022
####################################################################################

using Base
using HomotopyContinuation
using LinearAlgebra

#> 4-point problem has 12 unknowns and 45 coefficients

@var x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12;
@var p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20 p21 p22 p23 p24 p25 p26 p27 p28 p29 p30 p31 p32 p33 p34 p35 p36 p37 p38 p39 p40 p41 p42 p43 p44 p45;

#> polynomial equations
f1 = p1^2*x1^2 - 2*p1*p4*x1*x2 + 2*p37*p1*x1 + p2^2*x1^2 - 2*p2*p5*x1*x2 + 2*p38*p2*x1 + p3^2*x1^2 - 2*p3*p6*x1*x2 + 2*p39*p3*x1 + p4^2*x2^2 - 2*p37*p4*x2 + p5^2*x2^2 - 2*p38*p5*x2 + p6^2*x2^2 - 2*p39*p6*x2 - p13^2*x5^2 + 2*p13*p16*x5*x6 - 2*p37*p13*x5 - p14^2*x5^2 + 2*p14*p17*x5*x6 - 2*p38*p14*x5 - p15^2*x5^2 + 2*p15*p18*x5*x6 - 2*p39*p15*x5 - p16^2*x6^2 + 2*p37*p16*x6 - p17^2*x6^2 + 2*p38*p17*x6 - p18^2*x6^2 + 2*p39*p18*x6;
f2 = p1^2*x1^2 - 2*p1*p7*x1*x3 + 2*p40*p1*x1 + p2^2*x1^2 - 2*p2*p8*x1*x3 + 2*p41*p2*x1 + p3^2*x1^2 - 2*p3*p9*x1*x3 + 2*p42*p3*x1 + p7^2*x3^2 - 2*p40*p7*x3 + p8^2*x3^2 - 2*p41*p8*x3 + p9^2*x3^2 - 2*p42*p9*x3 - p13^2*x5^2 + 2*p13*p19*x5*x7 - 2*p40*p13*x5 - p14^2*x5^2 + 2*p14*p20*x5*x7 - 2*p41*p14*x5 - p15^2*x5^2 + 2*p15*p21*x5*x7 - 2*p42*p15*x5 - p19^2*x7^2 + 2*p40*p19*x7 - p20^2*x7^2 + 2*p41*p20*x7 - p21^2*x7^2 + 2*p42*p21*x7;
f3 = p1^2*x1^2 - 2*p1*p10*x1*x4 + 2*p43*p1*x1 + p2^2*x1^2 - 2*p2*p11*x1*x4 + 2*p44*p2*x1 + p3^2*x1^2 - 2*p3*p12*x1*x4 + 2*p45*p3*x1 + p10^2*x4^2 - 2*p43*p10*x4 + p11^2*x4^2 - 2*p44*p11*x4 + p12^2*x4^2 - 2*p45*p12*x4 - p13^2*x5^2 + 2*p13*p22*x5*x8 - 2*p43*p13*x5 - p14^2*x5^2 + 2*p14*p23*x5*x8 - 2*p44*p14*x5 - p15^2*x5^2 + 2*p15*p24*x5*x8 - 2*p45*p15*x5 - p22^2*x8^2 + 2*p43*p22*x8 - p23^2*x8^2 + 2*p44*p23*x8 - p24^2*x8^2 + 2*p45*p24*x8;
f4 = p4^2*x2^2 + p5^2*x2^2 + p6^2*x2^2 + p7^2*x3^2 + p8^2*x3^2 + p9^2*x3^2 - p16^2*x6^2 - p17^2*x6^2 - p18^2*x6^2 - p19^2*x7^2 - p20^2*x7^2 - p21^2*x7^2 - 2*p4*p37*x2 - 2*p5*p38*x2 + 2*p4*p40*x2 - 2*p6*p39*x2 + 2*p7*p37*x3 + 2*p5*p41*x2 + 2*p8*p38*x3 + 2*p6*p42*x2 - 2*p7*p40*x3 + 2*p9*p39*x3 - 2*p8*p41*x3 - 2*p9*p42*x3 + 2*p16*p37*x6 + 2*p17*p38*x6 - 2*p16*p40*x6 + 2*p18*p39*x6 - 2*p19*p37*x7 - 2*p17*p41*x6 - 2*p20*p38*x7 - 2*p18*p42*x6 + 2*p19*p40*x7 - 2*p21*p39*x7 + 2*p20*p41*x7 + 2*p21*p42*x7 - 2*p4*p7*x2*x3 - 2*p5*p8*x2*x3 - 2*p6*p9*x2*x3 + 2*p16*p19*x6*x7 + 2*p17*p20*x6*x7 + 2*p18*p21*x6*x7;
f5 = p4^2*x2^2 + p5^2*x2^2 + p6^2*x2^2 + p10^2*x4^2 + p11^2*x4^2 + p12^2*x4^2 - p16^2*x6^2 - p17^2*x6^2 - p18^2*x6^2 - p22^2*x8^2 - p23^2*x8^2 - p24^2*x8^2 - 2*p4*p37*x2 - 2*p5*p38*x2 - 2*p6*p39*x2 + 2*p4*p43*x2 + 2*p5*p44*x2 + 2*p10*p37*x4 + 2*p6*p45*x2 + 2*p11*p38*x4 + 2*p12*p39*x4 - 2*p10*p43*x4 - 2*p11*p44*x4 + 2*p16*p37*x6 - 2*p12*p45*x4 + 2*p17*p38*x6 + 2*p18*p39*x6 - 2*p16*p43*x6 - 2*p17*p44*x6 - 2*p22*p37*x8 - 2*p18*p45*x6 - 2*p23*p38*x8 - 2*p24*p39*x8 + 2*p22*p43*x8 + 2*p23*p44*x8 + 2*p24*p45*x8 - 2*p4*p10*x2*x4 - 2*p5*p11*x2*x4 - 2*p6*p12*x2*x4 + 2*p16*p22*x6*x8 + 2*p17*p23*x6*x8 + 2*p18*p24*x6*x8;
f6 = p7^2*x3^2 + p8^2*x3^2 + p9^2*x3^2 + p10^2*x4^2 + p11^2*x4^2 + p12^2*x4^2 - p19^2*x7^2 - p20^2*x7^2 - p21^2*x7^2 - p22^2*x8^2 - p23^2*x8^2 - p24^2*x8^2 - 2*p7*p40*x3 - 2*p8*p41*x3 + 2*p7*p43*x3 - 2*p9*p42*x3 + 2*p10*p40*x4 + 2*p8*p44*x3 + 2*p11*p41*x4 + 2*p9*p45*x3 - 2*p10*p43*x4 + 2*p12*p42*x4 - 2*p11*p44*x4 - 2*p12*p45*x4 + 2*p19*p40*x7 + 2*p20*p41*x7 - 2*p19*p43*x7 + 2*p21*p42*x7 - 2*p22*p40*x8 - 2*p20*p44*x7 - 2*p23*p41*x8 - 2*p21*p45*x7 + 2*p22*p43*x8 - 2*p24*p42*x8 + 2*p23*p44*x8 + 2*p24*p45*x8 - 2*p7*p10*x3*x4 - 2*p8*p11*x3*x4 - 2*p9*p12*x3*x4 + 2*p19*p22*x7*x8 + 2*p20*p23*x7*x8 + 2*p21*p24*x7*x8;
f7 = p1^2*x1^2 - 2*p1*p4*x1*x2 + 2*p37*p1*x1 + p2^2*x1^2 - 2*p2*p5*x1*x2 + 2*p38*p2*x1 + p3^2*x1^2 - 2*p3*p6*x1*x2 + 2*p39*p3*x1 + p4^2*x2^2 - 2*p37*p4*x2 + p5^2*x2^2 - 2*p38*p5*x2 + p6^2*x2^2 - 2*p39*p6*x2 - p25^2*x9^2 + 2*p25*p28*x9*x10 - 2*p37*p25*x9 - p26^2*x9^2 + 2*p26*p29*x9*x10 - 2*p38*p26*x9 - p27^2*x9^2 + 2*p27*p30*x9*x10 - 2*p39*p27*x9 - p28^2*x10^2 + 2*p37*p28*x10 - p29^2*x10^2 + 2*p38*p29*x10 - p30^2*x10^2 + 2*p39*p30*x10;
f8 = p1^2*x1^2 - 2*p1*p7*x1*x3 + 2*p40*p1*x1 + p2^2*x1^2 - 2*p2*p8*x1*x3 + 2*p41*p2*x1 + p3^2*x1^2 - 2*p3*p9*x1*x3 + 2*p42*p3*x1 + p7^2*x3^2 - 2*p40*p7*x3 + p8^2*x3^2 - 2*p41*p8*x3 + p9^2*x3^2 - 2*p42*p9*x3 - p25^2*x9^2 + 2*p25*p31*x9*x11 - 2*p40*p25*x9 - p26^2*x9^2 + 2*p26*p32*x9*x11 - 2*p41*p26*x9 - p27^2*x9^2 + 2*p27*p33*x9*x11 - 2*p42*p27*x9 - p31^2*x11^2 + 2*p40*p31*x11 - p32^2*x11^2 + 2*p41*p32*x11 - p33^2*x11^2 + 2*p42*p33*x11;
f9 = p1^2*x1^2 - 2*p1*p10*x1*x4 + 2*p43*p1*x1 + p2^2*x1^2 - 2*p2*p11*x1*x4 + 2*p44*p2*x1 + p3^2*x1^2 - 2*p3*p12*x1*x4 + 2*p45*p3*x1 + p10^2*x4^2 - 2*p43*p10*x4 + p11^2*x4^2 - 2*p44*p11*x4 + p12^2*x4^2 - 2*p45*p12*x4 - p25^2*x9^2 + 2*p25*p34*x9*x12 - 2*p43*p25*x9 - p26^2*x9^2 + 2*p26*p35*x9*x12 - 2*p44*p26*x9 - p27^2*x9^2 + 2*p27*p36*x9*x12 - 2*p45*p27*x9 - p34^2*x12^2 + 2*p43*p34*x12 - p35^2*x12^2 + 2*p44*p35*x12 - p36^2*x12^2 + 2*p45*p36*x12;
f10 = p4^2*x2^2 + p5^2*x2^2 + p6^2*x2^2 + p7^2*x3^2 + p8^2*x3^2 + p9^2*x3^2 - p28^2*x10^2 - p29^2*x10^2 - p30^2*x10^2 - p31^2*x11^2 - p32^2*x11^2 - p33^2*x11^2 - 2*p4*p37*x2 - 2*p5*p38*x2 + 2*p4*p40*x2 - 2*p6*p39*x2 + 2*p7*p37*x3 + 2*p5*p41*x2 + 2*p8*p38*x3 + 2*p6*p42*x2 - 2*p7*p40*x3 + 2*p9*p39*x3 - 2*p8*p41*x3 - 2*p9*p42*x3 + 2*p28*p37*x10 + 2*p29*p38*x10 - 2*p28*p40*x10 + 2*p30*p39*x10 - 2*p31*p37*x11 - 2*p29*p41*x10 - 2*p32*p38*x11 - 2*p30*p42*x10 + 2*p31*p40*x11 - 2*p33*p39*x11 + 2*p32*p41*x11 + 2*p33*p42*x11 - 2*p4*p7*x2*x3 - 2*p5*p8*x2*x3 - 2*p6*p9*x2*x3 + 2*p28*p31*x10*x11 + 2*p29*p32*x10*x11 + 2*p30*p33*x10*x11;
f11 = p4^2*x2^2 + p5^2*x2^2 + p6^2*x2^2 + p10^2*x4^2 + p11^2*x4^2 + p12^2*x4^2 - p28^2*x10^2 - p29^2*x10^2 - p30^2*x10^2 - p34^2*x12^2 - p35^2*x12^2 - p36^2*x12^2 - 2*p4*p37*x2 - 2*p5*p38*x2 - 2*p6*p39*x2 + 2*p4*p43*x2 + 2*p5*p44*x2 + 2*p10*p37*x4 + 2*p6*p45*x2 + 2*p11*p38*x4 + 2*p12*p39*x4 - 2*p10*p43*x4 - 2*p11*p44*x4 - 2*p12*p45*x4 + 2*p28*p37*x10 + 2*p29*p38*x10 + 2*p30*p39*x10 - 2*p28*p43*x10 - 2*p29*p44*x10 - 2*p34*p37*x12 - 2*p30*p45*x10 - 2*p35*p38*x12 - 2*p36*p39*x12 + 2*p34*p43*x12 + 2*p35*p44*x12 + 2*p36*p45*x12 - 2*p4*p10*x2*x4 - 2*p5*p11*x2*x4 - 2*p6*p12*x2*x4 + 2*p28*p34*x10*x12 + 2*p29*p35*x10*x12 + 2*p30*p36*x10*x12;
f12 = p7^2*x3^2 + p8^2*x3^2 + p9^2*x3^2 + p10^2*x4^2 + p11^2*x4^2 + p12^2*x4^2 - p31^2*x11^2 - p32^2*x11^2 - p33^2*x11^2 - p34^2*x12^2 - p35^2*x12^2 - p36^2*x12^2 - 2*p7*p40*x3 - 2*p8*p41*x3 + 2*p7*p43*x3 - 2*p9*p42*x3 + 2*p10*p40*x4 + 2*p8*p44*x3 + 2*p11*p41*x4 + 2*p9*p45*x3 - 2*p10*p43*x4 + 2*p12*p42*x4 - 2*p11*p44*x4 - 2*p12*p45*x4 + 2*p31*p40*x11 + 2*p32*p41*x11 - 2*p31*p43*x11 + 2*p33*p42*x11 - 2*p34*p40*x12 - 2*p32*p44*x11 - 2*p35*p41*x12 - 2*p33*p45*x11 + 2*p34*p43*x12 - 2*p36*p42*x12 + 2*p35*p44*x12 + 2*p36*p45*x12 - 2*p7*p10*x3*x4 - 2*p8*p11*x3*x4 - 2*p9*p12*x3*x4 + 2*p31*p34*x11*x12 + 2*p32*p35*x11*x12 + 2*p33*p36*x11*x12;

#> stack all polynomials
Eqs = [f1; f2; f3; f4; f5; f6; f7; f8; f9; f10; f11; f12]

F = System(Eqs; variables=[x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12], parameters=[p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45]);
S = monodromy_solve(F)
start_solutions = solutions(S);
start_params = S.parameters;

#> write the resultant start solutions to a file and store under the user-defined directory
#io = open("/your/start-sols/dir/julia-start-sols-raw", "w");
#using DelimitedFiles
#writedlm(io, start_solutions);
#close(io)
