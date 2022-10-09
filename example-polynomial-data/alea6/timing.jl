
# -- load packages --
using HomotopyContinuation
using Base
using LinearAlgebra
using BenchmarkTools

# -- decalre variables --
@var a[1:29]
@var x1, x2, x3, x4, x5, x6

f1 = a[1]*x1^2*x4 + a[2]*x2*x4*x5 + a[3]*x2*x4*x6 + a[4]*x4*x6 + a[5]*x5*x6
f2 = a[6]*x1*x2*x6 + a[7]*x2^2*x5 + a[8]*x2*x3*x5 + a[9]*x2*x5^2 + a[10]*x3^2 + a[11]*x4*x5*x6
f3 = a[12]*x1^2*x4 + a[13]*x1^2 + a[14]*x1*x4*x6 + a[13]*x2^2*x5 + a[3]*x2*x4^2 + a[15]*x6^3
f4 = a[16]*x1*x4^2 + a[14]*x2*x4 + a[17]*x2*x5^2 + a[18]*x3*x4 + a[19]*x3*x5*x6 + a[20]*x4*x5^2
f5 = a[21]*x1^2*x3 + a[22]*x1*x2*x3 + a[23]*x1*x4*x5 + a[24]*x1*x4 + a[25]*x2^2*x5
f6 = a[26]*x1*x3*x4 + a[19]*x3^2*x4 + a[27]*x3^2*x6 + a[28]*x3 + a[29]*x4^3 + a[25]*x5

Eq = [f1, f2, f3, f4, f5, f6]

start_coeffs = (rand(29)+rand(29)*im)

Eq_start = subs(Eq, a => start_coeffs);

F_start = System(Eq_start;variables=[x1, x2, x3, x4, x5, x6]);

start_solution = solve(F_start);

target_coeff = [5+0*im, 37+0*im, 32+0*im, 21+0*im, 55+0*im, 39+0*im, 23+0*im, 57+0*im, 56+0*im, 10+0*im, 52+0*im, 33+0*im, 51+0*im, 42+0*im, 1+0*im, 44+0*im, 47+0*im, 12+0*im, 2+0*im, 43+0*im, 49+0*im, 11+0*im, 83+0*im, 54+0*im, 45+0*im, 48+0*im, 59+0*im, 17+0*im, 36+0*im]

F = System(Eq;variables=[x1, x2, x3, x4, x5, x6], parameters=a)

sol_target = solve(F, solutions(start_solution); start_parameters=start_coeffs, target_parameters=target_coeff)
