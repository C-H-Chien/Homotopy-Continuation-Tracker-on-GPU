function [f, numOfVars, num_of_params] = sys_distorted_2view_triangulation()
    
    %> define systems
    %> variables
    syms x1 x2 x3 x4 x5
    %> parameters
    syms p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15

    numOfVars = 5;
    num_of_params = 15;
    
    %> polynomial equations
    f(1) = p13 + x1*p7 + x2*p10 + x3*p11 + x4*p12 + x1^2*p13*p14 + x2^2*p13*p14 + x3^2*p13*p15 + x4^2*p13*p15 + x1*x3*p5 + x1*x4*p6 + x2*x3*p8 + x2*x4*p9 + x1*x3^2*p7*p15 + x1*x4^2*p7*p15 + x1^2*x3*p11*p14 + x2*x3^2*p10*p15 + x2^2*x3*p11*p14 + x2*x4^2*p10*p15 + x1^2*x4*p12*p14 + x2^2*x4*p12*p14 + x1^2*x3^2*p13*p14*p15 + x1^2*x4^2*p13*p14*p15 + x2^2*x3^2*p13*p14*p15 + x2^2*x4^2*p13*p14*p15;
    f(2) = 2*x1 - 2*p1 + 2*x5*p7 + 2*x3*x5*p5 + 2*x4*x5*p6 + 4*x1*x5*p13*p14 + 2*x3^2*x5*p7*p15 + 2*x4^2*x5*p7*p15 + 4*x1*x3*x5*p11*p14 + 4*x1*x4*x5*p12*p14 + 4*x1*x3^2*x5*p13*p14*p15 + 4*x1*x4^2*x5*p13*p14*p15;
    f(3) = 2*x2 - 2*p2 + 2*x5*p10 + 2*x3*x5*p8 + 2*x4*x5*p9 + 4*x2*x5*p13*p14 + 2*x3^2*x5*p10*p15 + 2*x4^2*x5*p10*p15 + 4*x2*x3*x5*p11*p14 + 4*x2*x4*x5*p12*p14 + 4*x2*x3^2*x5*p13*p14*p15 + 4*x2*x4^2*x5*p13*p14*p15;
    f(4) = 2*x3 - 2*p3 + 2*x5*p7 + 2*x1*x5*p5 + 2*x2*x5*p6 + 4*x3*x5*p13*p15 + 2*x1^2*x5*p7*p14 + 2*x2^2*x5*p7*p14 + 4*x1*x3*x5*p11*p15 + 4*x2*x3*x5*p12*p15 + 4*x1^2*x3*x5*p13*p14*p15 + 4*x2^2*x3*x5*p13*p14*p15;
    f(5) = 2*x4 - 2*p4 + 2*x5*p10 + 2*x1*x5*p8 + 2*x2*x5*p9 + 4*x4*x5*p13*p15 + 2*x1^2*x5*p10*p14 + 2*x2^2*x5*p10*p14 + 4*x1*x4*x5*p11*p15 + 4*x2*x4*x5*p12*p15 + 4*x1^2*x4*x5*p13*p14*p15 + 4*x2^2*x4*x5*p13*p14*p15;
end
