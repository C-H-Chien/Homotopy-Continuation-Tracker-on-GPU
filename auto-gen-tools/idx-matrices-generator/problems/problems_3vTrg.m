function [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_3vTrg()
    % -- define the system -
    numOfVars = 9;
    numOfCoeff = 34;
    syms x_0 x_1 x_2 x_3 x_4 x_5 x_6 x_7 x_8
    X = [x_0 x_1 x_2 x_3 x_4 x_5 x_6 x_7 x_8];
    syms t;
    
    % -- coefficients --
    for i = 0:numOfCoeff-1
        % -- 1) target system --
        str_tgt_c = 'c_';
        cat_str_tgt_c = strcat(str_tgt_c, num2str(i));
        C(i+1) = str2sym(cat_str_tgt_c);

        % -- 2) start system --
        str_stt_c = 'd_';
        cat_str_stt_c = strcat(str_stt_c, num2str(i));
        D(i+1) = str2sym(cat_str_stt_c);
    end
    
    start_sys = [
        D(1)*x_0 + D(2)*x_2*x_6 + D(3)*x_3*x_6 + D(4)*x_4*x_8 + D(5)*x_5*x_8 + D(6)*x_6 + D(7)*x_8 + D(8)
        D(1)*x_1 + D(9)*x_2*x_6 + D(10)*x_3*x_6 + D(11)*x_4*x_8 + D(12)*x_5*x_8 + D(13)*x_6 + D(14)*x_8 + D(15)
        D(2)*x_0*x_6 + D(9)*x_1*x_6 + D(1)*x_2 + D(16)*x_4*x_7 + D(17)*x_5*x_7 + D(18)*x_6 + D(19)*x_7 + D(20)
        D(3)*x_0*x_6 + D(10)*x_1*x_6 + D(1)*x_3 + D(21)*x_4*x_7 + D(22)*x_5*x_7 + D(23)*x_6 + D(24)*x_7 + D(25)
        D(4)*x_0*x_8 + D(11)*x_1*x_8 + D(16)*x_2*x_7 + D(21)*x_3*x_7 + D(1)*x_4 + D(26)*x_7 + D(27)*x_8 + D(28)
        D(5)*x_0*x_8 + D(12)*x_1*x_8 + D(17)*x_2*x_7 + D(22)*x_3*x_7 + D(1)*x_5 + D(29)*x_7 + D(30)*x_8 + D(31)
        D(2)*x_0*x_2 + D(3)*x_0*x_3 + D(6)*x_0 + D(9)*x_1*x_2 + D(10)*x_1*x_3 + D(13)*x_1 + D(18)*x_2 + D(23)*x_3 + D(32)
        D(16)*x_2*x_4 + D(17)*x_2*x_5 + D(19)*x_2 + D(21)*x_3*x_4 + D(22)*x_3*x_5 + D(24)*x_3 + D(26)*x_4 + D(29)*x_5 + D(33)
        D(4)*x_0*x_4 + D(5)*x_0*x_5 + D(7)*x_0 + D(11)*x_1*x_4 + D(12)*x_1*x_5 + D(14)*x_1 + D(27)*x_4 + D(30)*x_5 + D(34)
    ];

    target_sys = [
        C(1)*x_0 + C(2)*x_2*x_6 + C(3)*x_3*x_6 + C(4)*x_4*x_8 + C(5)*x_5*x_8 + C(6)*x_6 + C(7)*x_8 + C(8)
        C(1)*x_1 + C(9)*x_2*x_6 + C(10)*x_3*x_6 + C(11)*x_4*x_8 + C(12)*x_5*x_8 + C(13)*x_6 + C(14)*x_8 + C(15)
        C(2)*x_0*x_6 + C(9)*x_1*x_6 + C(1)*x_2 + C(16)*x_4*x_7 + C(17)*x_5*x_7 + C(18)*x_6 + C(19)*x_7 + C(20)
        C(3)*x_0*x_6 + C(10)*x_1*x_6 + C(1)*x_3 + C(21)*x_4*x_7 + C(22)*x_5*x_7 + C(23)*x_6 + C(24)*x_7 + C(25)
        C(4)*x_0*x_8 + C(11)*x_1*x_8 + C(16)*x_2*x_7 + C(21)*x_3*x_7 + C(1)*x_4 + C(26)*x_7 + C(27)*x_8 + C(28)
        C(5)*x_0*x_8 + C(12)*x_1*x_8 + C(17)*x_2*x_7 + C(22)*x_3*x_7 + C(1)*x_5 + C(29)*x_7 + C(30)*x_8 + C(31)
        C(2)*x_0*x_2 + C(3)*x_0*x_3 + C(6)*x_0 + C(9)*x_1*x_2 + C(10)*x_1*x_3 + C(13)*x_1 + C(18)*x_2 + C(23)*x_3 + C(32)
        C(16)*x_2*x_4 + C(17)*x_2*x_5 + C(19)*x_2 + C(21)*x_3*x_4 + C(22)*x_3*x_5 + C(24)*x_3 + C(26)*x_4 + C(29)*x_5 + C(33)
        C(4)*x_0*x_4 + C(5)*x_0*x_5 + C(7)*x_0 + C(11)*x_1*x_4 + C(12)*x_1*x_5 + C(14)*x_1 + C(27)*x_4 + C(30)*x_5 + C(34)
    ];

    % -- construct the hpmotopy via the scalar variable t --
    Homotopy = start_sys * (1-t) + target_sys * t;

    % -- compute the Jacobian matrix --
    J = jacobian(Homotopy, [x_0 x_1 x_2 x_3 x_4 x_5 x_6 x_7 x_8 t]);
end