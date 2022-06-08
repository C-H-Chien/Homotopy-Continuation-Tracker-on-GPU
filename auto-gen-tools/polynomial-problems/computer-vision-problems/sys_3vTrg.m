function [f, numOfVars] = sys_3vTrg()
    % -- define systems --
    % -- variables --
    syms x1 x2 x3 x4 x5 x6 x7 x8 x9
    % -- parameters --
    syms c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 c16 c17 c18 c19 c20 c21 c22 c23 c24 c25 c26 c27 c28 c29 c30 c31 c32 c33
    
    numOfVars = 9;

    f(1) = 2*x1 + x3*x7*c1 + x4*x7*c2 + x5*x9*c19 + x6*x9*c20 + x7*c3 + x9*c21 - 2*c28;
    f(2) = 2*x2 + x3*x7*c4 + x4*x7*c5 + x5*x9*c22 + x6*x9*c23 + x7*c6 + x9*c24 - 2*c29;
    f(3) = x1*x7*c1 + x2*x7*c4 + 2*x3 + x5*x8*c10 + x6*x8*c11 + x7*c7 + x8*c12 - 2*c30;
    f(4) = x1*x7*c2 + x2*x7*c5 + 2*x4 + x5*x8*c13 + x6*x8*c14 + x7*c8 + x8*c15 - 2*c31;
    f(5) = x1*x9*c19 + x2*x9*c22 + x3*x8*c10 + x4*x8*c13 + 2*x5 + x8*c16 + x9*c25 - 2*c32;
    f(6) = x1*x9*c20 + x2*x9*c23 + x3*x8*c11 + x4*x8*c14 + 2*x6 + x8*c17 + x9*c26 - 2*c33;
    f(7) = x1*x3*c1 + x1*x4*c2 + x1*c3 + x2*x3*c4 + x2*x4*c5 + x2*c6 + x3*c7 + x4*c8 +c9;
    f(8) = x3*x5*c10 + x3*x6*c11 + x3*c12 + x4*x5*c13 + x4*x6*c14 + x4*c15 + x5*c16 + x6*c17 +c18;
    f(9) = x1*x5*c19 + x1*x6*c20 + x1*c21 + x2*x5*c22 + x2*x6*c23 + x2*c24 + x5*c25 + x6*c26 +c27;

end
