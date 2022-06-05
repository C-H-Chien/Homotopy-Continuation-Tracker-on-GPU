function [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_3vTrg_relax()
    % -- define the system -
    numOfVars = 8;
    numOfCoeff = 25;
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
        D(1)*x_0 + D(2)*x_2*x_6 + D(3)*x_3*x_6 + D(4)*x_6 + D(5)
        D(1)*x_1 + D(6)*x_2*x_6 + D(7)*x_3*x_6 + D(8)*x_6 + D(9)
        D(2)*x_0*x_6 + D(6)*x_1*x_6 + D(1)*x_2 + D(10)*x_4*x_7 + D(11)*x_5*x_7 + D(12)*x_6 + D(13)*x_7 + D(14)
        D(3)*x_0*x_6 + D(7)*x_1*x_6 + D(1)*x_3 + D(15)*x_4*x_7 + D(16)*x_5*x_7 + D(17)*x_6 + D(18)*x_7 + D(19)
        D(10)*x_2*x_7 + D(15)*x_3*x_7 + D(1)*x_4 + D(20)*x_7 + D(21)
        D(11)*x_2*x_7 + D(16)*x_3*x_7 + D(1)*x_5 + D(22)*x_7 + D(23)
        D(2)*x_0*x_2 + D(3)*x_0*x_3 + D(4)*x_0 + D(6)*x_1*x_2 + D(7)*x_1*x_3 + D(8)*x_1 + D(12)*x_2 + D(17)*x_3 + D(24)
        D(10)*x_2*x_4 + D(11)*x_2*x_5 + D(13)*x_2 + D(15)*x_3*x_4 + D(16)*x_3*x_5 + D(18)*x_3 + D(20)*x_4 + D(22)*x_5 + D(25)
    ];

    target_sys = [
        C(1)*x_0 + C(2)*x_2*x_6 + C(3)*x_3*x_6 + C(4)*x_6 + C(5)
        C(1)*x_1 + C(6)*x_2*x_6 + C(7)*x_3*x_6 + C(8)*x_6 + C(9)
        C(2)*x_0*x_6 + C(6)*x_1*x_6 + C(1)*x_2 + C(10)*x_4*x_7 + C(11)*x_5*x_7 + C(12)*x_6 + C(13)*x_7 + C(14)
        C(3)*x_0*x_6 + C(7)*x_1*x_6 + C(1)*x_3 + C(15)*x_4*x_7 + C(16)*x_5*x_7 + C(17)*x_6 + C(18)*x_7 + C(19)
        C(10)*x_2*x_7 + C(15)*x_3*x_7 + C(1)*x_4 + C(20)*x_7 + C(21)
        C(11)*x_2*x_7 + C(16)*x_3*x_7 + C(1)*x_5 + C(22)*x_7 + C(23)
        C(2)*x_0*x_2 + C(3)*x_0*x_3 + C(4)*x_0 + C(6)*x_1*x_2 + C(7)*x_1*x_3 + C(8)*x_1 + C(12)*x_2 + C(17)*x_3 + C(24)
        C(10)*x_2*x_4 + C(11)*x_2*x_5 + C(13)*x_2 + C(15)*x_3*x_4 + C(16)*x_3*x_5 + C(18)*x_3 + C(20)*x_4 + C(22)*x_5 + C(25)
    ];

    % -- construct the hpmotopy via the scalar variable t --
    Homotopy = start_sys * (1-t) + target_sys * t;

    % -- compute the Jacobian matrix --
    J = jacobian(Homotopy, [x_0 x_1 x_2 x_3 x_4 x_5 x_6 x_7 t]);
end