function [f, numOfVars] = sys_4vTrg()
    % -- define systems --
    % -- variables --
    syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11
    % -- parameters --
    syms c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 c16 c17 c18 c19 c20 c21 c22 c23 c24 c25 c26 c27 c28 c29 c30 c31 c32 c33 c34 c35
    
    numOfVars = 11;

    f(1) = 2*x1 + x3*x9*c1 + x4*x9*c2 + x9*c3 - 2*c28;
    f(2) = 2*x2 + x3*x9*c4 + x4*x9*c5 + x9*c6 - 2*c29;
    f(3) = x1*x9*c1 + x2*x9*c4 + 2*x3 + x5*x10*c10 + x6*x10*c11 + x9*c7 + x10*c12 - 2*c30;
    f(4) = x1*x9*c2 + x2*x9*c5 + 2*x4 + x5*x10*c13 + x6*x10*c14 + x9*c8 + x10*c15 - 2*c31;
    f(5) = x3*x10*c10 + x4*x10*c13 + 2*x5 + x7*x11*c19 + x8*x11*c20 + x10*c16 + x11*c21 - 2*c32;
    f(6) = x3*x10*c11 + x4*x10*c14 + 2*x6 + x7*x11*c22 + x8*x11*c23 + x10*c17 + x11*c24 - 2*c33;
    f(7) = x5*x11*c19 + x6*x11*c22 + 2*x7 + x11*c25 - 2*c34;
    f(8) = x5*x11*c20 + x6*x11*c23 + 2*x8 + x11*c26 - 2*c35;
    f(9) = x1*x3*c1 + x1*x4*c2 + x1*c3 + x2*x3*c4 + x2*x4*c5 + x2*c6 + x3*c7 + x4*c8 +c9;
    f(10) = x3*x5*c10 + x3*x6*c11 + x3*c12 + x4*x5*c13 + x4*x6*c14 + x4*c15 + x5*c16 + x6*c17 +c18;
    f(11) = x5*x7*c19 + x5*x8*c20 + x5*c21 + x6*x7*c22 + x6*x8*c23 + x6*c24 + x7*c25 + x8*c26 +c27;

end
