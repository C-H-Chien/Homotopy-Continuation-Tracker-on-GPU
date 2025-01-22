function [f, numOfVars, num_of_params] = sys_P3P_10x10()
    
    %> define systems
    %> variables
    syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10
    %> parameters
    syms p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15
    
    numOfVars = 10;
    num_of_params = 15;
    
    %> polynomial equations
    f(1) = x1 - p7 + x8*p1 - 2*x8*x5^2*p1 - 2*x8*x6^2*p1 + 2*x8*x4*x6 + 2*x8*x5*x7 + 2*x8*x4*x5*p2 - 2*x8*x6*x7*p2;
    f(2) = x2 - p8 + x8*p2 - 2*x8*x4^2*p2 - 2*x8*x6^2*p2 - 2*x8*x4*x7 + 2*x8*x5*x6 + 2*x8*x4*x5*p1 + 2*x8*x6*x7*p1;
    f(3) = x8 - p9 + x3 - 2*x8*x4^2 - 2*x8*x5^2 + 2*x8*x4*x6*p1 - 2*x8*x5*x7*p1 + 2*x8*x4*x7*p2 + 2*x8*x5*x6*p2;
    f(4) = x1 - p10 + x9*p3 - 2*x9*x5^2*p3 - 2*x9*x6^2*p3 + 2*x9*x4*x6 + 2*x9*x5*x7 + 2*x9*x4*x5*p4 - 2*x9*x6*x7*p4;
    f(5) = x2 - p11 + x9*p4 - 2*x9*x4^2*p4 - 2*x9*x6^2*p4 - 2*x9*x4*x7 + 2*x9*x5*x6 + 2*x9*x4*x5*p3 + 2*x9*x6*x7*p3;
    f(6) = x9 - p12 + x3 - 2*x9*x4^2 - 2*x9*x5^2 + 2*x9*x4*x6*p3 - 2*x9*x5*x7*p3 + 2*x9*x4*x7*p4 + 2*x9*x5*x6*p4;
    f(7) = x1 - p13 + x10*p5 - 2*x10*x5^2*p5 - 2*x10*x6^2*p5 + 2*x10*x4*x6 + 2*x10*x5*x7 + 2*x10*x4*x5*p6 - 2*x10*x6*x7*p6;
    f(8) = x2 - p14 + x10*p6 - 2*x10*x4^2*p6 - 2*x10*x6^2*p6 - 2*x10*x4*x7 + 2*x10*x5*x6 + 2*x10*x4*x5*p5 + 2*x10*x6*x7*p5;
    f(9) = x10 - p15 + x3 - 2*x10*x4^2 - 2*x10*x5^2 + 2*x10*x4*x6*p5 - 2*x10*x5*x7*p5 + 2*x10*x4*x7*p6 + 2*x10*x5*x6*p6;
    f(10) = x4^2 + x5^2 + x6^2 + x7^2 - 1;
end