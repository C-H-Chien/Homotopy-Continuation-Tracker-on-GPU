function [f, numOfVars] = sys_katsura13()
    % -- define systems --
    syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14
    numOfVars = 14;

    f(1) = x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2 + 2*x5^2 + 2*x6^2 + 2*x7^2 + 2*x8^2 + 2*x9^2 + 2*x10^2 + 2*x11^2 + 2*x12^2 + 2*x13^2 + 2*x14^2 - x1;
    f(2) = 2*x1*x2 + 2*x2*x3 + 2*x3*x4 + 2*x4*x5 + 2*x5*x6 + 2*x6*x7 + 2*x7*x8 + 2*x8*x9 + 2*x9*x10 + 2*x10*x11 + 2*x11*x12 + 2*x12*x13 + 2*x13*x14 - x2;
    f(3) = 2*x1*x3 + x2^2 + 2*x2*x4 + 2*x3*x5 + 2*x4*x6 + 2*x5*x7 + 2*x6*x8 + 2*x7*x9 + 2*x8*x10 + 2*x9*x11 + 2*x10*x12 + 2*x11*x13 + 2*x12*x14 - x3;
    f(4) = 2*x1*x4 + 2*x2*x3 + 2*x2*x5 + 2*x3*x6 + 2*x4*x7 + 2*x5*x8 + 2*x6*x9 + 2*x7*x10 + 2*x8*x11 + 2*x9*x12 + 2*x10*x13 + 2*x11*x14 - x4;
    f(5) = 2*x1*x5 + 2*x2*x4 + 2*x2*x6 + x3^2 + 2*x3*x7 + 2*x4*x8 + 2*x5*x9 + 2*x6*x10 + 2*x7*x11 + 2*x8*x12 + 2*x9*x13 + 2*x10*x14 - x5;
    f(6) = 2*x1*x6 + 2*x2*x5 + 2*x2*x7 + 2*x3*x4 + 2*x3*x8 + 2*x4*x9 + 2*x5*x10 + 2*x6*x11 + 2*x7*x12 + 2*x8*x13 + 2*x9*x14 - x6;
    f(7) = 2*x1*x7 + 2*x2*x6 + 2*x2*x8 + 2*x3*x5 + 2*x3*x9 + x4^2 + 2*x4*x10 + 2*x5*x11 + 2*x6*x12 + 2*x7*x13 + 2*x8*x14 - x7;
    f(8) = 2*x1*x8 + 2*x2*x7 + 2*x2*x9 + 2*x3*x6 + 2*x3*x10 + 2*x4*x5 + 2*x4*x11 + 2*x5*x12 + 2*x6*x13 + 2*x7*x14 - x8;
    f(9) = 2*x1*x9 + 2*x2*x8 + 2*x2*x10 + 2*x3*x7 + 2*x3*x11 + 2*x4*x6 + 2*x4*x12 + x5^2 + 2*x5*x13 + 2*x6*x14 - x9;
    f(10) = 2*x1*x10 + 2*x2*x9 + 2*x2*x11 + 2*x3*x8 + 2*x3*x12 + 2*x4*x7 + 2*x4*x13 + 2*x5*x6 + 2*x5*x14 - x10;
    f(11) = 2*x1*x11 + 2*x2*x10 + 2*x2*x12 + 2*x3*x9 + 2*x3*x13 + 2*x4*x8 + 2*x4*x14 + 2*x5*x7 + x6^2 - x11;
    f(12) = 2*x1*x12 + 2*x2*x11 + 2*x2*x13 + 2*x3*x10 + 2*x3*x14 + 2*x4*x9 + 2*x5*x8 + 2*x6*x7 - x12;
    f(13) = 2*x1*x13 + 2*x2*x12 + 2*x2*x14 + 2*x3*x11 + 2*x4*x10 + 2*x5*x9 + 2*x6*x8 + x7^2 - x13;
    f(14) = x1 + 2*x2 + 2*x3 + 2*x4 + 2*x5 + 2*x6 + 2*x7 + 2*x8 + 2*x9 + 2*x10 + 2*x11 + 2*x12 + 2*x13 + 2*x14 - 1;
    
end
