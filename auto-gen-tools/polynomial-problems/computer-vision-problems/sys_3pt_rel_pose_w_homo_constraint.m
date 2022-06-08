function [f, numOfVars] = sys_3pt_rel_pose_w_homo_constraint()
    % -- define systems --
    % -- variables --
    syms x1 x2 x3 x4 x5 x6 x7 x8
    % -- parameters --
    syms c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 c16 c17 c18
    
    numOfVars = 8;

    f(1) = x1*x6*x8*c1*c16 - x1*x6*c1*c13 + x2*x6*x8*c4*c16 - x2*x6*c4*c13 + x3*x6*x8*c7*c16 - x3*x6*c7*c13 - x4*c4*c16 - x5*c1*c16 +c7*c13;
    f(2) = -x1*x6*x7*c1*c16 + x1*x6*c1*c10 - x2*x6*x7*c4*c16 + x2*x6*c4*c10 - x3*x6*x7*c7*c16 + x3*x6*c7*c10 + x4*c1*c16 - x5*c4*c16 -c7*c10;
    f(3) = x1*x6*x8*c2*c17 - x1*x6*c2*c14 + x2*x6*x8*c5*c17 - x2*x6*c5*c14 + x3*x6*x8*c8*c17 - x3*x6*c8*c14 - x4*c5*c17 - x5*c2*c17 +c8*c14;
    f(4) = -x1*x6*x7*c2*c17 + x1*x6*c2*c11 - x2*x6*x7*c5*c17 + x2*x6*c5*c11 - x3*x6*x7*c8*c17 + x3*x6*c8*c11 + x4*c2*c17 - x5*c5*c17 -c8*c11;
    f(5) = x1*x6*x8*c3*c18 - x1*x6*c3*c15 + x2*x6*x8*c6*c18 - x2*x6*c6*c15 + x3*x6*x8*c9*c18 - x3*x6*c9*c15 - x4*c6*c18 - x5*c3*c18 +c9*c15;
    f(6) = -x1*x6*x7*c3*c18 + x1*x6*c3*c12 - x2*x6*x7*c6*c18 + x2*x6*c6*c12 - x3*x6*x7*c9*c18 + x3*x6*c9*c12 + x4*c3*c18 - x5*c6*c18 -c9*c12;
    f(7) = x1*x1 + x2*x2 + x3*x3 - 1;
    f(8) = x4*x4 + x5*x5 - 1;

end
