function [f, numOfVars] = sys_cyclic8()
    % -- define systems --
    syms x1 x2 x3 x4 x5 x6 x7 x8
    numOfVars = 8;

    f(1) = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8;
    f(2) = x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x6 + x6*x7 + x7*x8 + x8*x1;
    f(3) = x1*x2*x3 + x2*x3*x4 + x3*x4*x5 + x4*x5*x6 + x5*x6*x7 + x6*x7*x8 + x7*x8*x1 + x8*x1*x2;
    f(4) = x1*x2*x3*x4 + x2*x3*x4*x5 + x3*x4*x5*x6 + x4*x5*x6*x7 + x5*x6*x7*x8 + x6*x7*x8*x1 + x7*x8*x1*x2 + x8*x1*x2*x3;
    f(5) = x1*x2*x3*x4*x5 + x2*x3*x4*x5*x6 + x3*x4*x5*x6*x7 + x4*x5*x6*x7*x8 + x5*x6*x7*x8*x1 + x6*x7*x8*x1*x2 + x7*x8*x1*x2*x3 + x8*x1*x2*x3*x4;
    f(6) = x1*x2*x3*x4*x5*x6 + x2*x3*x4*x5*x6*x7 + x3*x4*x5*x6*x7*x8 + x4*x5*x6*x7*x8*x1 + x5*x6*x7*x8*x1*x2 + x6*x7*x8*x1*x2*x3 + x7*x8*x1*x2*x3*x4 + x8*x1*x2*x3*x4*x5;
    f(7) = x1*x2*x3*x4*x5*x6*x7 + x2*x3*x4*x5*x6*x7*x8 + x3*x4*x5*x6*x7*x8*x1 + x4*x5*x6*x7*x8*x1*x2 + x5*x6*x7*x8*x1*x2*x3 + x6*x7*x8*x1*x2*x3*x4 + x7*x8*x1*x2*x3*x4*x5 + x8*x1*x2*x3*x4*x5*x6;
    f(8) = x1*x2*x3*x4*x5*x6*x7*x8 - 1;
end
