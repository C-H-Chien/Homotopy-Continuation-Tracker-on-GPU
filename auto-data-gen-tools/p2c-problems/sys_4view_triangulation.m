function [f, numOfVars, num_of_params] = sys_4view_triangulation()
    % -- define systems --
    % -- variables --
    syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11
    % -- parameters --
    syms p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20 p21 p22 p23 p24 p25 p26 p27 p28 p29 p30 p31 p32 p33 p34 p35
    
    numOfVars = 11;
    num_of_params = 35;

    f(1)  = 2*x1 + p1*x3*x9 + p2*x4*x9 + p3*x9 - 2*p28*1;
    f(2)  = 2*x2 + p4*x3*x9 + p5*x4*x9 + p6*x9 - 2*p29*1;
    f(3)  = p1*x1*x9 + p4*x2*x9 + 2*x3 + p10*x5*x10 + p11*x6*x10 + p7*x9 + p12*x10 - 2*p30*1;
    f(4)  = p2*x1*x9 + p5*x2*x9 + 2*x4 + p13*x5*x10 + p14*x6*x10 + p8*x9 + p15*x10 - 2*p31*1;
    f(5)  = p10*x3*x10 + p13*x4*x10 + 2*x5 + p19*x7*x11 + p20*x8*x11 + p16*x10 + p21*x11 - 2*p32*1;
    f(6)  = p11*x3*x10 + p14*x4*x10 + 2*x6 + p22*x7*x11 + p23*x8*x11 + p17*x10 + p24*x11 - 2*p33*1;
    f(7)  = p19*x5*x11 + p22*x6*x11 + 2*x7 + p25*x11 - 2*p34*1;
    f(8)  = p20*x5*x11 + p23*x6*x11 + 2*x8 + p26*x11 - 2*p35*1;
    f(9)  = p1*x1*x3 + p2*x1*x4 + p3*x1 + p4*x2*x3 + p5*x2*x4 + p6*x2 + p7*x3 + p8*x4 + p9*1;
    f(10) = p10*x3*x5 + p11*x3*x6 + p12*x3 + p13*x4*x5 + p14*x4*x6 + p15*x4 + p16*x5 + p17*x6 + p18*1;
    f(11) = p19*x5*x7 + p20*x5*x8 + p21*x5 + p22*x6*x7 + p23*x6*x8 + p24*x6 + p25*x7 + p26*x8 + p27*1;

end
