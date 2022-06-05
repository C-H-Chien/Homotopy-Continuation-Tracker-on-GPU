function [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_3pt_rel_pose_w_homo_constraint()
    % -- define the system -
    numOfVars = 8;
    numOfCoeff = 44;
    syms x_0 x_1 x_2 x_3 x_4 x_5 x_6 x_7
    X = [x_0 x_1 x_2 x_3 x_4 x_5 x_6 x_7];
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
        D(1)*x_0*x_5*x_7 + D(2)*x_0*x_5 + D(3)*x_1*x_5*x_7 + D(4)*x_1*x_5 + D(5)*x_2*x_5*x_7 + D(6)*x_2*x_5 + D(7)*x_3 + D(8)*x_4 + D(9)
        D(8)*x_0*x_5*x_6 + D(10)*x_0*x_5 + D(7)*x_1*x_5*x_6 + D(11)*x_1*x_5 + D(12)*x_2*x_5*x_6 + D(13)*x_2*x_5 + D(1)*x_3 + D(7)*x_4 + D(14)
        D(15)*x_0*x_5*x_7 + D(16)*x_0*x_5 + D(17)*x_1*x_5*x_7 + D(18)*x_1*x_5 + D(19)*x_2*x_5*x_7 + D(20)*x_2*x_5 + D(21)*x_3 + D(22)*x_4 + D(23)
        D(22)*x_0*x_5*x_6 + D(24)*x_0*x_5 + D(21)*x_1*x_5*x_6 + D(25)*x_1*x_5 + D(26)*x_2*x_5*x_6 + D(27)*x_2*x_5 + D(15)*x_3 + D(21)*x_4 + D(28)
        D(29)*x_0*x_5*x_7 + D(30)*x_0*x_5 + D(31)*x_1*x_5*x_7 + D(32)*x_1*x_5 + D(33)*x_2*x_5*x_7 + D(34)*x_2*x_5 + D(35)*x_3 + D(36)*x_4 + D(37)
        D(36)*x_0*x_5*x_6 + D(38)*x_0*x_5 + D(35)*x_1*x_5*x_6 + D(39)*x_1*x_5 + D(40)*x_2*x_5*x_6 + D(41)*x_2*x_5 + D(29)*x_3 + D(35)*x_4 + D(42)
        D(43)*x_0^2 + D(43)*x_1^2 + D(43)*x_2^2 + D(44)
        D(43)*x_3^2 + D(43)*x_4^2 + D(44)
    ];

    target_sys = [
        C(1)*x_0*x_5*x_7 + C(2)*x_0*x_5 + C(3)*x_1*x_5*x_7 + C(4)*x_1*x_5 + C(5)*x_2*x_5*x_7 + C(6)*x_2*x_5 + C(7)*x_3 + C(8)*x_4 + C(9)
        C(8)*x_0*x_5*x_6 + C(10)*x_0*x_5 + C(7)*x_1*x_5*x_6 + C(11)*x_1*x_5 + C(12)*x_2*x_5*x_6 + C(13)*x_2*x_5 + C(1)*x_3 + C(7)*x_4 + C(14)
        C(15)*x_0*x_5*x_7 + C(16)*x_0*x_5 + C(17)*x_1*x_5*x_7 + C(18)*x_1*x_5 + C(19)*x_2*x_5*x_7 + C(20)*x_2*x_5 + C(21)*x_3 + C(22)*x_4 + C(23)
        C(22)*x_0*x_5*x_6 + C(24)*x_0*x_5 + C(21)*x_1*x_5*x_6 + C(25)*x_1*x_5 + C(26)*x_2*x_5*x_6 + C(27)*x_2*x_5 + C(15)*x_3 + C(21)*x_4 + C(28)
        C(29)*x_0*x_5*x_7 + C(30)*x_0*x_5 + C(31)*x_1*x_5*x_7 + C(32)*x_1*x_5 + C(33)*x_2*x_5*x_7 + C(34)*x_2*x_5 + C(35)*x_3 + C(36)*x_4 + C(37)
        C(36)*x_0*x_5*x_6 + C(38)*x_0*x_5 + C(35)*x_1*x_5*x_6 + C(39)*x_1*x_5 + C(40)*x_2*x_5*x_6 + C(41)*x_2*x_5 + C(29)*x_3 + C(35)*x_4 + C(42)
        C(43)*x_0^2 + C(43)*x_1^2 + C(43)*x_2^2 + C(44)
        C(43)*x_3^2 + C(43)*x_4^2 + C(44)
    ];

    % -- construct the hpmotopy via the scalar variable t --
    Homotopy = start_sys * (1-t) + target_sys * t;

    % -- compute the Jacobian matrix --
    J = jacobian(Homotopy, [x_0 x_1 x_2 x_3 x_4 x_5 x_6 x_7 t]);
end