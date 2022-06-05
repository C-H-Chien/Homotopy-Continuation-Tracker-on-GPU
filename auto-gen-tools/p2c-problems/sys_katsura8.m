function [f, numOfVars] = sys_katsura8()
    % -- define systems --
    syms x1 x2 x3 x4 x5 x6 x7 x8 x9
    numOfVars = 9;

    f(1) = -x1+2*x9^2+2*x8^2+2*x7^2+2*x6^2+2*x5^2+2*x4^2+2*x3^2+2*x2^2+x1^2;
    f(2) = -x2+2*x9*x8+2*x8*x7+2*x7*x6+2*x6*x5+2*x5*x4+2*x4*x3+2*x3*x2+2*x2*x1;
    f(3) = -x3+2*x9*x7+2*x8*x6+2*x7*x5+2*x6*x4+2*x5*x3+2*x4*x2+2*x3*x1+x2^2;
    f(4) = -x4+2*x9*x6+2*x8*x5+2*x7*x4+2*x6*x3+2*x5*x2+2*x4*x1+2*x3*x2;
    f(5) = -x5+2*x9*x5+2*x8*x4+2*x7*x3+2*x6*x2+2*x5*x1+2*x4*x2+x3^2;
    f(6) = -x6+2*x9*x4+2*x8*x3+2*x7*x2+2*x6*x1+2*x5*x2+2*x4*x3;
    f(7) = -x7+2*x9*x3+2*x8*x2+2*x7*x1+2*x6*x2+2*x5*x3+x4^2;
    f(8) = -x8+2*x9*x2+2*x8*x1+2*x7*x2+2*x6*x3+2*x5*x4;
    f(9) = -1+2*x9+2*x8+2*x7+2*x6+2*x5+2*x4+2*x3+2*x2+x1;
end
