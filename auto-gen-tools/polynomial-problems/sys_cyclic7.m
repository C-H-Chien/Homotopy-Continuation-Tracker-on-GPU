function [f, numOfVars] = sys_cyclic7()
    % -- define systems --
    syms x1 x2 x3 x4 x5 x6 x7
    numOfVars = 7;

    f(1) = x1+x2+x3+x4+x5+x6+x7;
    f(2) = x1*x2+x1*x7+x2*x3+x3*x4+x4*x5+x5*x6+x6*x7;
    f(3) = x1*x2*x3+x1*x2*x7+x1*x6*x7+x2*x3*x4+x3*x4*x5+x4*x5*x6+x5*x6*x7;
    f(4) = x1*x2*x3*x4+x1*x2*x3*x7+x1*x2*x6*x7+x1*x5*x6*x7+x2*x3*x4*x5+x3*x4*x5*x6+x4*x5*x6*x7;
    f(5) = x1*x2*x3*x4*x5+x1*x2*x3*x4*x7+x1*x2*x3*x6*x7+x1*x2*x5*x6*x7+x1*x4*x5*x6*x7+x2*x3*x4*x5*x6+x3*x4*x5*x6*x7;
    f(6) = x1*x2*x3*x4*x5*x6+x1*x2*x3*x4*x5*x7+x1*x2*x3*x4*x6*x7+x1*x2*x3*x5*x6*x7+x1*x2*x4*x5*x6*x7+x1*x3*x4*x5*x6*x7+x2*x3*x4*x5*x6*x7;
    f(7) = x1*x2*x3*x4*x5*x6*x7-1;
end
