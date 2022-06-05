function [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_PnP_wo_principal_point()
    % -- define the system -
    numOfVars = 10;
    numOfCoeff = 34;
    syms x_0 x_1 x_2 x_3 x_4 x_5 x_6 x_7 x_8 x_9
    X = [x_0 x_1 x_2 x_3 x_4 x_5 x_6 x_7 x_8 x_9];
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
        D(1)*x_0 + D(2)*x_4*x_6 + D(3)*x_4*x_7 + D(4)*x_4*x_8 + D(5)*x_4*x_9
        D(1)*x_1 + D(6)*x_4*x_6 + D(7)*x_4*x_7 + D(8)*x_4*x_8 + D(9)*x_4*x_9
        D(1)*x_2 + D(10)*x_5*x_6 + D(11)*x_5*x_7 + D(12)*x_5*x_8 + D(13)*x_5*x_9
        D(1)*x_3 + D(14)*x_5*x_6 + D(15)*x_5*x_7 + D(16)*x_5*x_8 + D(17)*x_5*x_9
        D(2)*x_0*x_6 + D(3)*x_0*x_7 + D(4)*x_0*x_8 + D(5)*x_0*x_9 + D(6)*x_1*x_6 + D(7)*x_1*x_7 + D(8)*x_1*x_8 + D(9)*x_1*x_9 + D(18)*x_6 + D(19)*x_7 + D(20)*x_8 + D(21)*x_9
        D(10)*x_2*x_6 + D(11)*x_2*x_7 + D(12)*x_2*x_8 + D(13)*x_2*x_9 + D(14)*x_3*x_6 + D(15)*x_3*x_7 + D(16)*x_3*x_8 + D(17)*x_3*x_9 + D(22)*x_6 + D(23)*x_7 + D(24)*x_8 + D(25)*x_9
        D(2)*x_0*x_4 + D(6)*x_1*x_4 + D(10)*x_2*x_5 + D(14)*x_3*x_5 + D(18)*x_4 + D(22)*x_5 + D(26)
        D(3)*x_0*x_4 + D(7)*x_1*x_4 + D(11)*x_2*x_5 + D(15)*x_3*x_5 + D(19)*x_4 + D(23)*x_5 + D(27)
        D(4)*x_0*x_4 + D(8)*x_1*x_4 + D(12)*x_2*x_5 + D(16)*x_3*x_5 + D(20)*x_4 + D(24)*x_5 + D(28)
        D(5)*x_0*x_4 + D(9)*x_1*x_4 + D(13)*x_2*x_5 + D(17)*x_3*x_5 + D(21)*x_4 + D(25)*x_5 + D(29)
    ];

    target_sys = [
        C(1)*x_0 + C(2)*x_4*x_6 + C(3)*x_4*x_7 + C(4)*x_4*x_8 + C(5)*x_4*x_9
        C(1)*x_1 + C(6)*x_4*x_6 + C(7)*x_4*x_7 + C(8)*x_4*x_8 + C(9)*x_4*x_9
        C(1)*x_2 + C(10)*x_5*x_6 + C(11)*x_5*x_7 + C(12)*x_5*x_8 + C(13)*x_5*x_9
        C(1)*x_3 + C(14)*x_5*x_6 + C(15)*x_5*x_7 + C(16)*x_5*x_8 + C(17)*x_5*x_9
        C(2)*x_0*x_6 + C(3)*x_0*x_7 + C(4)*x_0*x_8 + C(5)*x_0*x_9 + C(6)*x_1*x_6 + C(7)*x_1*x_7 + C(8)*x_1*x_8 + C(9)*x_1*x_9 + C(18)*x_6 + C(19)*x_7 + C(20)*x_8 + C(21)*x_9
        C(10)*x_2*x_6 + C(11)*x_2*x_7 + C(12)*x_2*x_8 + C(13)*x_2*x_9 + C(14)*x_3*x_6 + C(15)*x_3*x_7 + C(16)*x_3*x_8 + C(17)*x_3*x_9 + C(22)*x_6 + C(23)*x_7 + C(24)*x_8 + C(25)*x_9
        C(2)*x_0*x_4 + C(6)*x_1*x_4 + C(10)*x_2*x_5 + C(14)*x_3*x_5 + C(18)*x_4 + C(22)*x_5 + C(26)
        C(3)*x_0*x_4 + C(7)*x_1*x_4 + C(11)*x_2*x_5 + C(15)*x_3*x_5 + C(19)*x_4 + C(23)*x_5 + C(27)
        C(4)*x_0*x_4 + C(8)*x_1*x_4 + C(12)*x_2*x_5 + C(16)*x_3*x_5 + C(20)*x_4 + C(24)*x_5 + C(28)
        C(5)*x_0*x_4 + C(9)*x_1*x_4 + C(13)*x_2*x_5 + C(17)*x_3*x_5 + C(21)*x_4 + C(25)*x_5 + C(29)
    ];

    % -- construct the hpmotopy via the scalar variable t --
    Homotopy = start_sys * (1-t) + target_sys * t;

    % -- compute the Jacobian matrix --
    J = jacobian(Homotopy, [x_0 x_1 x_2 x_3 x_4 x_5 x_6 x_7 x_8 x_9 t]);
end