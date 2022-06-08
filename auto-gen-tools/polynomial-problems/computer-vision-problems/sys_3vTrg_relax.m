function [f, numOfVars] = sys_3vTrg_relax()
    % -- define systems --
    % -- variables --
    syms x1 x2 x3 x4 x5 x6 x7 x8
    % -- parameters --
    syms c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 c16 c17 c18 c19 c20 c21 c22 c23 c24
    
    numOfVars = 8;

    f(1) = 2*x1 + x3*x7*c1 + x4*x7*c4 + x7*c7 - 2*c19;
    f(2) = 2*x2 + x3*x7*c2 + x4*x7*c5 + x7*c8 - 2*c20;
    f(3) = x1*x7*c1 + x2*x7*c2 + 2*x3 + x5*x8*c10 + x6*x8*c13 + x7*c3 + x8*c16 - 2*c21;
    f(4) = x1*x7*c4 + x2*x7*c5 + 2*x4 + x5*x8*c11 + x6*x8*c14 + x7*c6 + x8*c17 - 2*c22;
    f(5) = x3*x8*c10 + x4*x8*c11 + 2*x5 + x8*c12 - 2*c23;
    f(6) = x3*x8*c13 + x4*x8*c14 + 2*x6 + x8*c15 - 2*c24;
    f(7) = x1*x3*c1 + x1*x4*c4 + x1*c7 + x2*x3*c2 + x2*x4*c5 + x2*c8 + x3*c3 + x4*c6 +c9;
    f(8) = x3*x5*c10 + x3*x6*c13 + x3*c16 + x4*x5*c11 + x4*x6*c14 + x4*c17 + x5*c12 + x6*c15 +c18;

end
