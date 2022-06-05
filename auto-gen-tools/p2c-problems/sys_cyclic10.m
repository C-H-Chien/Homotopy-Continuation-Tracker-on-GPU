function [f, numOfVars] = sys_cyclic10()
    % -- define systems --
    syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10
    numOfVars = 10;

    f(1) = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10;
    f(2) = x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x6 + x6*x7 + x7*x8 + x8*x9 + x9*x10 + x10*x1;
    f(3) = x1*x2*x3 + x2*x3*x4 + x3*x4*x5 + x4*x5*x6 + x5*x6*x7 + x6*x7*x8 + x7*x8*x9 + x8*x9*x10 + x9*x10*x1 + x10*x1*x2;
    f(4) = x1*x2*x3*x4 + x2*x3*x4*x5 + x3*x4*x5*x6 + x4*x5*x6*x7 + x5*x6*x7*x8 + x6*x7*x8*x9 + x7*x8*x9*x10 + x8*x9*x10*x1 + x9*x10*x1*x2 + x10*x1*x2*x3;
    f(5) = x1*x2*x3*x4*x5 + x2*x3*x4*x5*x6 + x3*x4*x5*x6*x7 + x4*x5*x6*x7*x8 + x5*x6*x7*x8*x9 + x6*x7*x8*x9*x10 + x7*x8*x9*x10*x1 + x8*x9*x10*x1*x2 + x9*x10*x1*x2*x3 + x10*x1*x2*x3*x4;
    f(6) = x1*x2*x3*x4*x5*x6 + x2*x3*x4*x5*x6*x7 + x3*x4*x5*x6*x7*x8 + x4*x5*x6*x7*x8*x9 + x5*x6*x7*x8*x9*x10 + x6*x7*x8*x9*x10*x1 + x7*x8*x9*x10*x1*x2 + x8*x9*x10*x1*x2*x3 + x9*x10*x1*x2*x3*x4 + x10*x1*x2*x3*x4*x5;
    f(7) = x1*x2*x3*x4*x5*x6*x7 + x2*x3*x4*x5*x6*x7*x8 + x3*x4*x5*x6*x7*x8*x9 + x4*x5*x6*x7*x8*x9*x10 + x5*x6*x7*x8*x9*x10*x1 + x6*x7*x8*x9*x10*x1*x2 + x7*x8*x9*x10*x1*x2*x3 + x8*x9*x10*x1*x2*x3*x4 + x9*x10*x1*x2*x3*x4*x5 + x10*x1*x2*x3*x4*x5*x6;
    f(8) = x1*x2*x3*x4*x5*x6*x7*x8 + x2*x3*x4*x5*x6*x7*x8*x9 + x3*x4*x5*x6*x7*x8*x9*x10 + x4*x5*x6*x7*x8*x9*x10*x1 + x5*x6*x7*x8*x9*x10*x1*x2 + x6*x7*x8*x9*x10*x1*x2*x3 + x7*x8*x9*x10*x1*x2*x3*x4 + x8*x9*x10*x1*x2*x3*x4*x5 + x9*x10*x1*x2*x3*x4*x5*x6 + x10*x1*x2*x3*x4*x5*x6*x7;
    f(9) = x1*x2*x3*x4*x5*x6*x7*x8*x9 + x2*x3*x4*x5*x6*x7*x8*x9*x10 + x3*x4*x5*x6*x7*x8*x9*x10*x1 + x4*x5*x6*x7*x8*x9*x10*x1*x2 + x5*x6*x7*x8*x9*x10*x1*x2*x3 + x6*x7*x8*x9*x10*x1*x2*x3*x4 + x7*x8*x9*x10*x1*x2*x3*x4*x5 + x8*x9*x10*x1*x2*x3*x4*x5*x6 + x9*x10*x1*x2*x3*x4*x5*x6*x7 + x10*x1*x2*x3*x4*x5*x6*x7*x8;
    f(10) = x1*x2*x3*x4*x5*x6*x7*x8*x9*x10 - 1;
end
