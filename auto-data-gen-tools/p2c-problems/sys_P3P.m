function [f, numOfVars, num_of_params] = sys_P3P()
    
    %> define systems
    %> variables
    syms x1 x2 x3
    %> parameters
    syms p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18
    
    numOfVars = 3;
    num_of_params = 18;
    
    %> polynomial equations
    f(1) = x1^2*p1^2 + x1^2*p2^2 + x1^2*p3^2 - 2*x1*x2*p1*p4 - 2*x1*x2*p2*p5 - 2*x1*x2*p3*p6 + x2^2*p4^2 + x2^2*p5^2 + x2^2*p6^2 - p10^2 + 2*p10*p13 - p11^2 + 2*p11*p14 - p12^2 + 2*p12*p15 - p13^2 - p14^2 - p15^2;
    f(2) = x1^2*p1^2 + x1^2*p2^2 + x1^2*p3^2 - 2*x1*x3*p1*p7 - 2*x1*x3*p2*p8 - 2*x1*x3*p3*p9 + x3^2*p7^2 + x3^2*p8^2 + x3^2*p9^2 - p10^2 + 2*p10*p16 - p11^2 + 2*p11*p17 - p12^2 + 2*p12*p18 - p16^2 - p17^2 - p18^2;
    f(3) = x2^2*p4^2 + x2^2*p5^2 + x2^2*p6^2 - 2*x2*x3*p4*p7 - 2*x2*x3*p5*p8 - 2*x2*x3*p6*p9 + x3^2*p7^2 + x3^2*p8^2 + x3^2*p9^2 - p13^2 + 2*p13*p16 - p14^2 + 2*p14*p17 - p15^2 + 2*p15*p18 - p16^2 - p17^2 - p18^2;
end