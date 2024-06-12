function [f, numOfVars, num_of_params] = sys_rel_pos_quiver()
    
    %> define systems
    %> variables
    syms x1 x2 x3 x4
    %> parameters
    syms p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20 p21 p22 p23 p24 p25 p26 p27 p28 p29 p30 p31 p32 p33 p34 p35 p36
    
    numOfVars = 4;
    num_of_params = 36;
    
    %> polynomial equations
    f(1) = -x1^2*x4*p33 + x1^2*p1 - x1^2*p17 + 2*x1*x2*p5 + 2*x1*x2*p13 + 2*x1*x3*x4*p9 + 2*x1*x3*p25 + 2*x1*x4*p21 - 2*x1*p29 - x2^2*x4*p33 - x2^2*p1 + x2^2*p17 + 2*x2*x3*x4*p21 + 2*x2*x3*p29 - 2*x2*x4*p9 + 2*x2*p25 + x3^2*x4*p33 - x3^2*p1 - x3^2*p17 + 2*x3*p5 - 2*x3*p13 + x4*p33 + p1 + p17;
    f(2) = -x1^2*x4*p34 + x1^2*p2 - x1^2*p18 + 2*x1*x2*p6 + 2*x1*x2*p14 + 2*x1*x3*x4*p10 + 2*x1*x3*p26 + 2*x1*x4*p22 - 2*x1*p30 - x2^2*x4*p34 - x2^2*p2 + x2^2*p18 + 2*x2*x3*x4*p22 + 2*x2*x3*p30 - 2*x2*x4*p10 + 2*x2*p26 + x3^2*x4*p34 - x3^2*p2 - x3^2*p18 + 2*x3*p6 - 2*x3*p14 + x4*p34 + p2 + p18;
    f(3) = -x1^2*x4*p35 + x1^2*p3 - x1^2*p19 + 2*x1*x2*p7 + 2*x1*x2*p15 + 2*x1*x3*x4*p11 + 2*x1*x3*p27 + 2*x1*x4*p23 - 2*x1*p31 - x2^2*x4*p35 - x2^2*p3 + x2^2*p19 + 2*x2*x3*x4*p23 + 2*x2*x3*p31 - 2*x2*x4*p11 + 2*x2*p27 + x3^2*x4*p35 - x3^2*p3 - x3^2*p19 + 2*x3*p7 - 2*x3*p15 + x4*p35 + p3 + p19;
    f(4) = -x1^2*x4*p36 + x1^2*p4 - x1^2*p20 + 2*x1*x2*p8 + 2*x1*x2*p16 + 2*x1*x3*x4*p12 + 2*x1*x3*p28 + 2*x1*x4*p24 - 2*x1*p32 - x2^2*x4*p36 - x2^2*p4 + x2^2*p20 + 2*x2*x3*x4*p24 + 2*x2*x3*p32 - 2*x2*x4*p12 + 2*x2*p28 + x3^2*x4*p36 - x3^2*p4 - x3^2*p20 + 2*x3*p8 - 2*x3*p16 + x4*p36 + p4 + p20;
end