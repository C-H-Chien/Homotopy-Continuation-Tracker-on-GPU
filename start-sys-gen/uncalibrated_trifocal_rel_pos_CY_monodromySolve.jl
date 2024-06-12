using Base
using HomotopyContinuation
using LinearAlgebra

@var s a2 a3 a4 b1 b2 b3 b4 c1 c2 c3 c4 k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11 k12 k13 k14 k15 k16 k17 k18 k19 k20 k21 k22 k23 k24

f1 = a2^2*s + a2^2*k2^2 + a2^2*k6^2 - 2*a2*s - 2*a2*k1*k2 - 2*a2*k5*k6 - b1^2*s - b1^2*k9^2 - b1^2*k13^2 + 2*b1*b2*s + 2*b1*b2*k9*k10 + 2*b1*b2*k13*k14 - b2^2*s - b2^2*k10^2 - b2^2*k14^2 + s + k1^2 + k5^2;
f2 = a3^2*s + a3^2*k3^2 + a3^2*k7^2 - 2*a3*s - 2*a3*k1*k3 - 2*a3*k5*k7 - b1^2*s - b1^2*k9^2 - b1^2*k13^2 + 2*b1*b3*s + 2*b1*b3*k9*k11 + 2*b1*b3*k13*k15 - b3^2*s - b3^2*k11^2 - b3^2*k15^2 + s + k1^2 + k5^2;
f3 = a2^2*s + a2^2*k2^2 + a2^2*k6^2 - 2*a2*a3*s - 2*a2*a3*k2*k3 - 2*a2*a3*k6*k7 + a3^2*s + a3^2*k3^2 + a3^2*k7^2 - b2^2*s - b2^2*k10^2 - b2^2*k14^2 + 2*b2*b3*s + 2*b2*b3*k10*k11 + 2*b2*b3*k14*k15 - b3^2*s - b3^2*k11^2 - b3^2*k15^2;
f4 = a2^2*s + a2^2*k2^2 + a2^2*k6^2 - 2*a2*a4*s - 2*a2*a4*k2*k4 - 2*a2*a4*k6*k8 + a4^2*s + a4^2*k4^2 + a4^2*k8^2 - b2^2*s - b2^2*k10^2 - b2^2*k14^2 + 2*b2*b4*s + 2*b2*b4*k10*k12 + 2*b2*b4*k14*k16 - b4^2*s - b4^2*k12^2 - b4^2*k16^2;
f5 = a4^2*s + a4^2*k4^2 + a4^2*k8^2 - 2*a4*s - 2*a4*k1*k4 - 2*a4*k5*k8 - b1^2*s - b1^2*k9^2 - b1^2*k13^2 + 2*b1*b4*s + 2*b1*b4*k9*k12 + 2*b1*b4*k13*k16 - b4^2*s - b4^2*k12^2 - b4^2*k16^2 + s + k1^2 + k5^2;
f6 = a3^2*s + a3^2*k3^2 + a3^2*k7^2 - 2*a3*a4*s - 2*a3*a4*k3*k4 - 2*a3*a4*k7*k8 + a4^2*s + a4^2*k4^2 + a4^2*k8^2 - b3^2*s - b3^2*k11^2 - b3^2*k15^2 + 2*b3*b4*s + 2*b3*b4*k11*k12 + 2*b3*b4*k15*k16 - b4^2*s - b4^2*k12^2 - b4^2*k16^2;
f7 = a2^2*s + a2^2*k2^2 + a2^2*k6^2 - 2*a2*s - 2*a2*k1*k2 - 2*a2*k5*k6 - c1^2*s - c1^2*k17^2 - c1^2*k21^2 + 2*c1*c2*s + 2*c1*c2*k17*k18 + 2*c1*c2*k21*k22 - c2^2*s - c2^2*k18^2 - c2^2*k22^2 + s + k1^2 + k5^2;
f8 = a3^2*s + a3^2*k3^2 + a3^2*k7^2 - 2*a3*s - 2*a3*k1*k3 - 2*a3*k5*k7 - c1^2*s - c1^2*k17^2 - c1^2*k21^2 + 2*c1*c3*s + 2*c1*c3*k17*k19 + 2*c1*c3*k21*k23 - c3^2*s - c3^2*k19^2 - c3^2*k23^2 + s + k1^2 + k5^2;
f9 = a2^2*s + a2^2*k2^2 + a2^2*k6^2 - 2*a2*a3*s - 2*a2*a3*k2*k3 - 2*a2*a3*k6*k7 + a3^2*s + a3^2*k3^2 + a3^2*k7^2 - c2^2*s - c2^2*k18^2 - c2^2*k22^2 + 2*c2*c3*s + 2*c2*c3*k18*k19 + 2*c2*c3*k22*k23 - c3^2*s - c3^2*k19^2 - c3^2*k23^2;
f10 = a2^2*s + a2^2*k2^2 + a2^2*k6^2 - 2*a2*a4*s - 2*a2*a4*k2*k4 - 2*a2*a4*k6*k8 + a4^2*s + a4^2*k4^2 + a4^2*k8^2 - c2^2*s - c2^2*k18^2 - c2^2*k22^2 + 2*c2*c4*s + 2*c2*c4*k18*k20 + 2*c2*c4*k22*k24 - c4^2*s - c4^2*k20^2 - c4^2*k24^2;
f11 = a4^2*s + a4^2*k4^2 + a4^2*k8^2 - 2*a4*s - 2*a4*k1*k4 - 2*a4*k5*k8 - c1^2*s - c1^2*k17^2 - c1^2*k21^2 + 2*c1*c4*s + 2*c1*c4*k17*k20 + 2*c1*c4*k21*k24 - c4^2*s - c4^2*k20^2 - c4^2*k24^2 + s + k1^2 + k5^2;
f12 = a3^2*s + a3^2*k3^2 + a3^2*k7^2 - 2*a3*a4*s - 2*a3*a4*k3*k4 - 2*a3*a4*k7*k8 + a4^2*s + a4^2*k4^2 + a4^2*k8^2 - c3^2*s - c3^2*k19^2 - c3^2*k23^2 + 2*c3*c4*s + 2*c3*c4*k19*k20 + 2*c3*c4*k23*k24 - c4^2*s - c4^2*k20^2 - c4^2*k24^2;

F = System([f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12]; variables=[s, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4], parameters = [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, k18, k19, k20, k21, k22, k23, k24]);

S = monodromy_solve(F)
start_solutions = solutions(S);
start_params = S.parameters;

target_params = rand(24) + rand(24) * im;
@time for i = 1:10
solve(F, start_solutions; start_parameters=start_params, target_parameters=target_params)
end

#> write target solutions to a file
#io = open("/home/chchien/hcOutput", "w");
#using DelimitedFiles
#writedlm(io, solutions(R));
#close(io)
