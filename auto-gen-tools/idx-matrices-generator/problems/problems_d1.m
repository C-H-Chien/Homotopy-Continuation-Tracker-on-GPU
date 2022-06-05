function [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_d1()
    % -- define the system -
    numOfVars = 12;
    numOfCoeff = 10;
    syms x_0 x_1 x_2 x_3 x_4 x_5 x_6 x_7 x_8 x_9 x_10 x_11
    X = [x_0 x_1 x_2 x_3 x_4 x_5 x_6 x_7 x_8 x_9 x_10 x_11];
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
        D(1)*x_0^2 + D(1)*x_1^2 + D(2)*1
        D(1)*x_2^2 + D(1)*x_3^2 + D(2)*1
        D(1)*x_4^2 + D(1)*x_5^2 + D(2)*1
        D(1)*x_6^2 + D(1)*x_7^2 + D(2)*1
        D(1)*x_8^2 + D(1)*x_9^2 + D(2)*1
        D(1)*x_10^2 + D(1)*x_11^2 + D(2)*1
        D(3)*x_2 + D(4)*x_4 + D(1)*x_6 + D(5)*1
        D(3)*x_0*x_3 + D(4)*x_0*x_5 + D(1)*x_0*x_7 + D(6)*1
        D(3)*x_1*x_3 + D(4)*x_1*x_5 + D(1)*x_1*x_7 + D(7)*1
        D(1)*x_2*x_8 + D(1)*x_4*x_8 + D(1)*x_6*x_8 + D(8)*1
        D(1)*x_0*x_9 + D(1)*x_1*x_3*x_8 + D(1)*x_1*x_5*x_8 + D(1)*x_1*x_7*x_8 + D(9)*1
        D(2)*x_2*x_9*x_10 + D(1)*x_3*x_11 + D(2)*x_4*x_9*x_10 + D(1)*x_5*x_11 + D(2)*x_6*x_9*x_10 + D(1)*x_7*x_11 + D(10)
    ];

    target_sys = [
        C(1)*x_0^2 + C(1)*x_1^2 + C(2)*1
        C(1)*x_2^2 + C(1)*x_3^2 + C(2)*1
        C(1)*x_4^2 + C(1)*x_5^2 + C(2)*1
        C(1)*x_6^2 + C(1)*x_7^2 + C(2)*1
        C(1)*x_8^2 + C(1)*x_9^2 + C(2)*1
        C(1)*x_10^2 + C(1)*x_11^2 + C(2)*1
        C(3)*x_2 + C(4)*x_4 + C(1)*x_6 + C(5)*1
        C(3)*x_0*x_3 + C(4)*x_0*x_5 + C(1)*x_0*x_7 + C(6)*1
        C(3)*x_1*x_3 + C(4)*x_1*x_5 + C(1)*x_1*x_7 + C(7)*1
        C(1)*x_2*x_8 + C(1)*x_4*x_8 + C(1)*x_6*x_8 + C(8)*1
        C(1)*x_0*x_9 + C(1)*x_1*x_3*x_8 + C(1)*x_1*x_5*x_8 + C(1)*x_1*x_7*x_8 + C(9)*1
        C(2)*x_2*x_9*x_10 + C(1)*x_3*x_11 + C(2)*x_4*x_9*x_10 + C(1)*x_5*x_11 + C(2)*x_6*x_9*x_10 + C(1)*x_7*x_11 + C(10)
    ];

    % -- construct the hpmotopy via the scalar variable t --
    Homotopy = start_sys * (1-t) + target_sys * t;

    % -- compute the Jacobian matrix --
    J = jacobian(Homotopy, [x_0 x_1 x_2 x_3 x_4 x_5 x_6 x_7 x_8 x_9 x_10 x_11 t]);
end