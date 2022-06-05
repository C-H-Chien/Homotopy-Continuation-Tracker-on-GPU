function [f, numOfVars] = sys_d1()
    % -- define systems --
    syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12
    numOfVars = 12;

    f(1) = x1^2  + x2^2 - 1;
    f(2) = x3^2  + x4^2 - 1;
    f(3) = x5^2  + x6^2 - 1;
    f(4) = x7^2  + x8^2 - 1;
    f(5) = x9^2  + x10^2 - 1;
    f(6) = x11^2 + x12^2 - 1;
    f(7) = 3*x3 + 2*x5 + x7 - 3.9701;
    f(8) = 3*x1*x4 + 2*x1*x6 + x1*x8 - 1.7172;
    f(9) = 3*x2*x4 + 2*x2*x6 + x2*x8 - 4.0616;
    f(10) = x3*x9 + x5*x9 + x7*x9 - 1.9791;
    f(11) = x2*x4*x9 + x2*x6*x9 + x2*x8*x9 + x1*x10 - 1.9115;
    f(12) = -x3*x10*x11 - x5*x10*x11 - x7*x10*x11 + x4*x12 + x6*x12 + x8*x12 - 0.4077;
end
