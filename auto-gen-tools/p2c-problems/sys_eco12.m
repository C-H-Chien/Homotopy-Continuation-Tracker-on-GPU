function [f, numOfVars] = sys_eco12()
    % -- define systems --
    syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12
    numOfVars = 12;

    f(1) = x1*x12*x2 + x1*x12 + x10*x11*x12 + x10*x12*x9 + x12*x2*x3 + x12*x3*x4 + x12* x4*x5 + x12*x5*x6 + x12*x6*x7 + x12*x7*x8 + x12*x8*x9 - 1;
    f(2) = x1*x12*x3 + x10*x12*x8 + x11*x12*x9 + x12*x2*x4 + x12*x2 + x12*x3*x5 + x12* x4*x6 + x12*x5*x7 + x12*x6*x8 + x12*x7*x9 - 2;
    f(3) = x1*x12*x4 + x10*x12*x7 + x11*x12*x8 + x12*x2*x5 + x12*x3*x6 + x12*x3 + x12* x4*x7 + x12*x5*x8 + x12*x6*x9 - 3;
    f(4) = x1*x12*x5 + x10*x12*x6 + x11*x12*x7 + x12*x2*x6 + x12*x3*x7 + x12*x4*x8 + x12*x4 + x12*x5*x9 - 4;
    f(5) = x1*x12*x6 + x10*x12*x5 + x11*x12*x6 + x12*x2*x7 + x12*x3*x8 + x12*x4*x9 + x12*x5 - 5;
    f(6) = x1*x12*x7 + x10*x12*x4 + x11*x12*x5 + x12*x2*x8 + x12*x3*x9 + x12*x6 - 6;
    f(7) = x1*x12*x8 + x10*x12*x3 + x11*x12*x4 + x12*x2*x9 + x12*x7 - 7;
    f(8) = x1*x12*x9 + x10*x12*x2 + x11*x12*x3 + x12*x8 - 8;
    f(9) = x1*x10*x12 + x11*x12*x2 + x12*x9 - 9;
    f(10) = x1*x11*x12 + x10*x12 - 10;
    f(11) = x11*x12 - 11;
    f(12) = x1 + x10 + x11 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + 1;
end
