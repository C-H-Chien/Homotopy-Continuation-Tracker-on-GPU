function [f, numOfVars] = sys_katsura6()
    % -- define systems --
    syms x1 x2 x3 x4 x5 x6 x7
    numOfVars = 7;

    f(1) = 1*x1+2*x2+2*x3+2*x4+2*x5+2*x6+2*x7-1;
    f(2) = 2*x4*x3+2*x5*x2+2*x6*x1+2*x7*x2-1*x6;
    f(3) = 1*x3^2+2*x4*x2+2*x5*x1+2*x6*x2+2*x7*x3-1*x5;
    f(4) = 2*x3*x2+2*x4*x1+2*x5*x2+2*x6*x3+2*x7*x4-1*x4;
    f(5) = 1*x2^2+2*x3*x1+2*x4*x2+2*x5*x3+2*x6*x4+2*x7*x5-1*x3;
    f(6) = 2*x2*x1+2*x3*x2+2*x4*x3+2*x5*x4+2*x6*x5+2*x7*x6-1*x2;
    f(7) = 1*x1^2+2*x2^2+2*x3^2+2*x4^2+2*x5^2+2*x6^2+2*x7^2-1*x1;
end
