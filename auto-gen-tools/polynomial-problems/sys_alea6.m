function [f, numOfVars] = sys_alea6()
    % -- define systems --
    syms x1 x2 x3 x4 x5 x6
    numOfVars = 6;

    f(1) = 5*x1^2*x4+37*x2*x4*x5+32*x2*x4*x6+21*x4*x6+55*x5*x6;
    f(2) = 39*x1*x2*x6+23*x2^2*x5+57*x2*x3*x5+56*x2*x5^2+10*x3^2+52*x4*x5*x6;
    f(3) = 33*x1^2*x4+51*x1^2+42*x1*x4*x6+51*x2^2*x5+32*x2*x4^2+x6^3;
    f(4) = 44*x1*x4^2+42*x2*x4+47*x2*x5^2+12*x3*x4+2*x3*x5*x6+43*x4*x5^2;
    f(5) = 49*x1^2*x3+11*x1*x2*x3+39*x1*x4*x5+44*x1*x4*x5+54*x1*x4+45*x2^2*x5;
    f(6) = 48*x1*x3*x4+2*x3^2*x4+59*x3^2*x6+17*x3+36*x4^3+45*x5;
end