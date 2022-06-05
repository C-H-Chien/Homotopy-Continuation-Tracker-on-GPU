function [f, numOfVars] = sys_refractive_p5p()
    % -- define systems --
    % -- variables --
    syms x1 x2 x3 x4 x5
    % -- parameters --
    syms c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 c16 c17 c18 c19 c20 c21 c22 c23 c24 c25 c26 c27 c28 c29
    
    numOfVars = 5;

    f(1) = -x1* x3*x3*c2 + 2*x1*x3*x4*c1 + 2*x1*x3*c3 + x1* x4*x4*c2 + 2*x1*x4*x5*c3 - x1* x5*x5*c2 - 2*x1*x5*c1 + x1*c2 - x2* x3*x3*c1 - 2*x2*x3*x4*c2 - 2*x2*x3*x5*c3 + x2* x4*x4*c1 + 2*x2*x4*c3 + x2* x5*x5*c1 - 2*x2*x5*c2 - x2*c1 + x3*x3*c1*c17 + x3*x3*c2*c16 - 2*x3*x4*c1*c16 + 2*x3*x4*c2*c17 + 2*x3*x5*c3*c17 - 2*x3*c3*c16 - x4*x4*c1*c17 - x4*x4*c2*c16 - 2*x4*x5*c3*c16 - 2*x4*c3*c17 - x5*x5*c1*c17 + x5*x5*c2*c16 + 2*x5*c1*c16 + 2*x5*c2*c17 +c1*c17 -c2*c16;
    f(2) = -x1* x3*x3*c5 + 2*x1*x3*x4*c4 + 2*x1*x3*c6 + x1* x4*x4*c5 + 2*x1*x4*x5*c6 - x1* x5*x5*c5 - 2*x1*x5*c4 + x1*c5 - x2* x3*x3*c4 - 2*x2*x3*x4*c5 - 2*x2*x3*x5*c6 + x2* x4*x4*c4 + 2*x2*x4*c6 + x2* x5*x5*c4 - 2*x2*x5*c5 - x2*c4 + x3*x3*c4*c20 + x3*x3*c5*c19 - 2*x3*x4*c4*c19 + 2*x3*x4*c5*c20 + 2*x3*x5*c6*c20 - 2*x3*c6*c19 - x4*x4*c4*c20 - x4*x4*c5*c19 - 2*x4*x5*c6*c19 - 2*x4*c6*c20 - x5*x5*c4*c20 + x5*x5*c5*c19 + 2*x5*c4*c19 + 2*x5*c5*c20 +c4*c20 -c5*c19;
    f(3) = -x1* x3*x3*c8 + 2*x1*x3*x4*c7 + 2*x1*x3*c9 + x1* x4*x4*c8 + 2*x1*x4*x5*c9 - x1* x5*x5*c8 - 2*x1*x5*c7 + x1*c8 - x2* x3*x3*c7 - 2*x2*x3*x4*c8 - 2*x2*x3*x5*c9 + x2* x4*x4*c7 + 2*x2*x4*c9 + x2* x5*x5*c7 - 2*x2*x5*c8 - x2*c7 + x3*x3*c7*c23 + x3*x3*c8*c22 - 2*x3*x4*c7*c22 + 2*x3*x4*c8*c23 + 2*x3*x5*c9*c23 - 2*x3*c9*c22 - x4*x4*c7*c23 - x4*x4*c8*c22 - 2*x4*x5*c9*c22 - 2*x4*c9*c23 - x5*x5*c7*c23 + x5*x5*c8*c22 + 2*x5*c7*c22 + 2*x5*c8*c23 +c7*c23 -c8*c22;
    f(4) = -x1* x3*x3*c11 + 2*x1*x3*x4*c10 + 2*x1*x3*c12 + x1* x4*x4*c11 + 2*x1*x4*x5*c12 - x1* x5*x5*c11 - 2*x1*x5*c10 + x1*c11 - x2* x3*x3*c10 - 2*x2*x3*x4*c11 - 2*x2*x3*x5*c12 + x2* x4*x4*c10 + 2*x2*x4*c12 + x2* x5*x5*c10 - 2*x2*x5*c11 - x2*c10 + x3*x3*c10*c26 + x3*x3*c11*c25 - 2*x3*x4*c10*c25 + 2*x3*x4*c11*c26 + 2*x3*x5*c12*c26 - 2*x3*c12*c25 - x4*x4*c10*c26 - x4*x4*c11*c25 - 2*x4*x5*c12*c25 - 2*x4*c12*c26 - x5*x5*c10*c26 + x5*x5*c11*c25 + 2*x5*c10*c25 + 2*x5*c11*c26 +c10*c26 -c11*c25;
    f(5) = -x1* x3*x3*c14 + 2*x1*x3*x4*c13 + 2*x1*x3*c15 + x1* x4*x4*c14 + 2*x1*x4*x5*c15 - x1* x5*x5*c14 - 2*x1*x5*c13 + x1*c14 - x2* x3*x3*c13 - 2*x2*x3*x4*c14 - 2*x2*x3*x5*c15 + x2* x4*x4*c13 + 2*x2*x4*c15 + x2* x5*x5*c13 - 2*x2*x5*c14 - x2*c13 + x3*x3*c13*c29 + x3*x3*c14*c28 - 2*x3*x4*c13*c28 + 2*x3*x4*c14*c29 + 2*x3*x5*c15*c29 - 2*x3*c15*c28 - x4*x4*c13*c29 - x4*x4*c14*c28 - 2*x4*x5*c15*c28 - 2*x4*c15*c29 - x5*x5*c13*c29 + x5*x5*c14*c28 + 2*x5*c13*c28 + 2*x5*c14*c29 +c13*c29 -c14*c28;
end
