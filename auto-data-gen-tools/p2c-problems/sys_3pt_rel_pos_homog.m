function [f, numOfVars] = sys_3pt_rel_pos_homog()
    % -- define systems --
    % -- variables --
    syms x1 x2 x3 x4 x5 x6 x7 x8
    % -- parameters --
    syms p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18
    
    numOfVars = 8;

    f(1) = x1*x6*x8*p1*p16 - x1*x6*p1*p13 + x2*x6*x8*p4*p16 - x2*x6*p4*p13 + x3*x6*x8*p7*p16 - x3*x6*p7*p13 - x4*p4*p16 - x5*p1*p16 + p7*p13;
    f(2) = -x1*x6*x7*p1*p16 + x1*x6*p1*p10 - x2*x6*x7*p4*p16 + x2*x6*p4*p10 - x3*x6*x7*p7*p16 + x3*x6*p7*p10 + x4*p1*p16 - x5*p4*p16 - p7*p10;
    f(3) = x1*x6*x8*p2*p17 - x1*x6*p2*p14 + x2*x6*x8*p5*p17 - x2*x6*p5*p14 + x3*x6*x8*p8*p17 - x3*x6*p8*p14 - x4*p5*p17 - x5*p2*p17 + p8*p14;
    f(4) = -x1*x6*x7*p2*p17 + x1*x6*p2*p11 - x2*x6*x7*p5*p17 + x2*x6*p5*p11 - x3*x6*x7*p8*p17 + x3*x6*p8*p11 + x4*p2*p17 - x5*p5*p17 - p8*p11;
    f(5) = x1*x6*x8*p3*p18 - x1*x6*p3*p15 + x2*x6*x8*p6*p18 - x2*x6*p6*p15 + x3*x6*x8*p9*p18 - x3*x6*p9*p15 - x4*p6*p18 - x5*p3*p18 + p9*p15;
    f(6) = -x1*x6*x7*p3*p18 + x1*x6*p3*p12 - x2*x6*x7*p6*p18 + x2*x6*p6*p12 - x3*x6*x7*p9*p18 + x3*x6*p9*p12 + x4*p3*p18 - x5*p6*p18 - p9*p12;
    f(7) = x1^2 + x2^2 + x3^2 - 1;
    f(8) = x4^2 + x5^2 - 1;

end
