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

#> six-line problem has 14 unknowns and 36 coefficients

@var x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16;
@var p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20 p21 p22 p23 p24 p25 p26 p27 p28 p29 p30 p31 p32 p33 p34 p35 p36;

#> polynomial equations
f1 = p13*x13 - p25*x5 + x2*x13 - x5*x10 + p2*p25*x8 - p13*p25*x3 - p14*p25*x4 - p2*p13*x16 + p13*p25*x11 + p13*p26*x12 - p2*x2*x16 + p2*x8*x10 + p1*x5*x16 - p1*x8*x13 - p13*x3*x10 + p14*x1*x13 - p14*x4*x10 + p25*x2*x11 + p26*x2*x12 - p26*x5*x9 + p2*p13*p25*x6 + p2*p14*p25*x7 - p2*p13*p25*x14 - p2*p13*p26*x15 + p2*p13*x6*x10 + p1*p13*x3*x16 - p1*p13*x6*x13 - p2*p14*x1*x16 + p2*p14*x7*x10 + p1*p14*x4*x16 - p1*p14*x7*x13 - p2*p25*x2*x14 + p1*p25*x5*x14 - p1*p25*x8*x11 - p2*p26*x2*x15 + p2*p26*x8*x9 + p1*p26*x5*x15 - p1*p26*x8*x12 - p13*p26*x3*x9 + p14*p25*x1*x11 + p14*p26*x1*x12 - p14*p26*x4*x9 + p1*p13*p25*x3*x14 - p1*p13*p25*x6*x11 + p2*p13*p26*x6*x9 - p2*p14*p25*x1*x14 + p1*p13*p26*x3*x15 - p1*p13*p26*x6*x12 + p1*p14*p25*x4*x14 - p1*p14*p25*x7*x11 - p2*p14*p26*x1*x15 + p2*p14*p26*x7*x9 + p1*p14*p26*x4*x15 - p1*p14*p26*x7*x12;
f2 = p15*x13 - p27*x5 + x2*x13 - x5*x10 + p4*p27*x8 - p15*p27*x3 - p16*p27*x4 - p4*p15*x16 + p15*p27*x11 + p15*p28*x12 - p4*x2*x16 + p4*x8*x10 + p3*x5*x16 - p3*x8*x13 - p15*x3*x10 + p16*x1*x13 - p16*x4*x10 + p27*x2*x11 + p28*x2*x12 - p28*x5*x9 + p4*p15*p27*x6 + p4*p16*p27*x7 - p4*p15*p27*x14 - p4*p15*p28*x15 + p4*p15*x6*x10 + p3*p15*x3*x16 - p3*p15*x6*x13 - p4*p16*x1*x16 + p4*p16*x7*x10 + p3*p16*x4*x16 - p3*p16*x7*x13 - p4*p27*x2*x14 + p3*p27*x5*x14 - p3*p27*x8*x11 - p4*p28*x2*x15 + p4*p28*x8*x9 + p3*p28*x5*x15 - p3*p28*x8*x12 - p15*p28*x3*x9 + p16*p27*x1*x11 + p16*p28*x1*x12 - p16*p28*x4*x9 + p3*p15*p27*x3*x14 - p3*p15*p27*x6*x11 + p4*p15*p28*x6*x9 - p4*p16*p27*x1*x14 + p3*p15*p28*x3*x15 - p3*p15*p28*x6*x12 + p3*p16*p27*x4*x14 - p3*p16*p27*x7*x11 - p4*p16*p28*x1*x15 + p4*p16*p28*x7*x9 + p3*p16*p28*x4*x15 - p3*p16*p28*x7*x12;
f3 = p17*x13 - p29*x5 + x2*x13 - x5*x10 + p6*p29*x8 - p17*p29*x3 - p18*p29*x4 - p6*p17*x16 + p17*p29*x11 + p17*p30*x12 - p6*x2*x16 + p6*x8*x10 + p5*x5*x16 - p5*x8*x13 - p17*x3*x10 + p18*x1*x13 - p18*x4*x10 + p29*x2*x11 + p30*x2*x12 - p30*x5*x9 + p6*p17*p29*x6 + p6*p18*p29*x7 - p6*p17*p29*x14 - p6*p17*p30*x15 + p6*p17*x6*x10 + p5*p17*x3*x16 - p5*p17*x6*x13 - p6*p18*x1*x16 + p6*p18*x7*x10 + p5*p18*x4*x16 - p5*p18*x7*x13 - p6*p29*x2*x14 + p5*p29*x5*x14 - p5*p29*x8*x11 - p6*p30*x2*x15 + p6*p30*x8*x9 + p5*p30*x5*x15 - p5*p30*x8*x12 - p17*p30*x3*x9 + p18*p29*x1*x11 + p18*p30*x1*x12 - p18*p30*x4*x9 + p5*p17*p29*x3*x14 - p5*p17*p29*x6*x11 + p6*p17*p30*x6*x9 - p6*p18*p29*x1*x14 + p5*p17*p30*x3*x15 - p5*p17*p30*x6*x12 + p5*p18*p29*x4*x14 - p5*p18*p29*x7*x11 - p6*p18*p30*x1*x15 + p6*p18*p30*x7*x9 + p5*p18*p30*x4*x15 - p5*p18*p30*x7*x12;
f4 = p19*x13 - p31*x5 + x2*x13 - x5*x10 + p8*p31*x8 - p19*p31*x3 - p20*p31*x4 - p8*p19*x16 + p19*p31*x11 + p19*p32*x12 - p8*x2*x16 + p8*x8*x10 + p7*x5*x16 - p7*x8*x13 - p19*x3*x10 + p20*x1*x13 - p20*x4*x10 + p31*x2*x11 + p32*x2*x12 - p32*x5*x9 + p8*p19*p31*x6 + p8*p20*p31*x7 - p8*p19*p31*x14 - p8*p19*p32*x15 + p8*p19*x6*x10 + p7*p19*x3*x16 - p7*p19*x6*x13 - p8*p20*x1*x16 + p8*p20*x7*x10 + p7*p20*x4*x16 - p7*p20*x7*x13 - p8*p31*x2*x14 + p7*p31*x5*x14 - p7*p31*x8*x11 - p8*p32*x2*x15 + p8*p32*x8*x9 + p7*p32*x5*x15 - p7*p32*x8*x12 - p19*p32*x3*x9 + p20*p31*x1*x11 + p20*p32*x1*x12 - p20*p32*x4*x9 + p7*p19*p31*x3*x14 - p7*p19*p31*x6*x11 + p8*p19*p32*x6*x9 - p8*p20*p31*x1*x14 + p7*p19*p32*x3*x15 - p7*p19*p32*x6*x12 + p7*p20*p31*x4*x14 - p7*p20*p31*x7*x11 - p8*p20*p32*x1*x15 + p8*p20*p32*x7*x9 + p7*p20*p32*x4*x15 - p7*p20*p32*x7*x12;
f5 = p21*x13 - p33*x5 + x2*x13 - x5*x10 + p10*p33*x8 - p21*p33*x3 - p22*p33*x4 - p10*p21*x16 + p21*p33*x11 + p21*p34*x12 - p10*x2*x16 + p10*x8*x10 + p9*x5*x16 - p9*x8*x13 - p21*x3*x10 + p22*x1*x13 - p22*x4*x10 + p33*x2*x11 + p34*x2*x12 - p34*x5*x9 + p10*p21*p33*x6 + p10*p22*p33*x7 - p10*p21*p33*x14 - p10*p21*p34*x15 + p10*p21*x6*x10 + p9*p21*x3*x16 - p9*p21*x6*x13 - p10*p22*x1*x16 + p10*p22*x7*x10 + p9*p22*x4*x16 - p9*p22*x7*x13 - p10*p33*x2*x14 + p9*p33*x5*x14 - p9*p33*x8*x11 - p10*p34*x2*x15 + p10*p34*x8*x9 + p9*p34*x5*x15 - p9*p34*x8*x12 - p21*p34*x3*x9 + p22*p33*x1*x11 + p22*p34*x1*x12 - p22*p34*x4*x9 + p9*p21*p33*x3*x14 - p9*p21*p33*x6*x11 + p10*p21*p34*x6*x9 - p10*p22*p33*x1*x14 + p9*p21*p34*x3*x15 - p9*p21*p34*x6*x12 + p9*p22*p33*x4*x14 - p9*p22*p33*x7*x11 - p10*p22*p34*x1*x15 + p10*p22*p34*x7*x9 + p9*p22*p34*x4*x15 - p9*p22*p34*x7*x12;
f6 = p23*x13 - p35*x5 + x2*x13 - x5*x10 + p12*p35*x8 - p23*p35*x3 - p24*p35*x4 - p12*p23*x16 + p23*p35*x11 + p23*p36*x12 - p12*x2*x16 + p12*x8*x10 + p11*x5*x16 - p11*x8*x13 - p23*x3*x10 + p24*x1*x13 - p24*x4*x10 + p35*x2*x11 + p36*x2*x12 - p36*x5*x9 + p12*p23*p35*x6 + p12*p24*p35*x7 - p12*p23*p35*x14 - p12*p23*p36*x15 + p12*p23*x6*x10 + p11*p23*x3*x16 - p11*p23*x6*x13 - p12*p24*x1*x16 + p12*p24*x7*x10 + p11*p24*x4*x16 - p11*p24*x7*x13 - p12*p35*x2*x14 + p11*p35*x5*x14 - p11*p35*x8*x11 - p12*p36*x2*x15 + p12*p36*x8*x9 + p11*p36*x5*x15 - p11*p36*x8*x12 - p23*p36*x3*x9 + p24*p35*x1*x11 + p24*p36*x1*x12 - p24*p36*x4*x9 + p11*p23*p35*x3*x14 - p11*p23*p35*x6*x11 + p12*p23*p36*x6*x9 - p12*p24*p35*x1*x14 + p11*p23*p36*x3*x15 - p11*p23*p36*x6*x12 + p11*p24*p35*x4*x14 - p11*p24*p35*x7*x11 - p12*p24*p36*x1*x15 + p12*p24*p36*x7*x9 + p11*p24*p36*x4*x15 - p11*p24*p36*x7*x12;
f7 = x1*x7+x2*x8+x6;
f8 = x1*x4+x2*x5+x3;
f9 = x3*x5+x6*x8+x2;
f10 = x4^2+x5^2-x6^2-1;
f11 = x3*x4+x6*x7+x1;
f12 = x9*x15+x10*x16+x14;
f13 = x9*x12+x10*x13+x11;
f14 = x11*x13+x14*x16+x10;
f15 = x12^2+x13^2-x14^2-1;
f16 = x11*x12+x14*x15+x9;


#> stack all polynomials
Eqs = [f1; f2; f3; f4; f5; f6; f7; f8; f9; f10; f11; f12; f13; f14; f15; f16]

#> make a System
F = System(Eqs; variables=[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16], parameters=[p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36]);
S = monodromy_solve(F)
start_solutions = solutions(S);
start_params = S.parameters;

#> write the resultant start solutions to a file and store under the user-defined directory
#io = open("/your/start-sols/dir/julia-start-sols-raw", "w");
#using DelimitedFiles
#writedlm(io, start_solutions);
#close(io)
