function [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_cyclic8()
    % -- define the system -
    numOfVars = 8;
    numOfCoeff = 2;
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
        D(1)*x_0 + D(1)*x_1 + D(1)*x_2 + D(1)*x_3 + D(1)*x_4 + D(1)*x_5 + D(1)*x_6 + D(1)*x_7
        D(1)*x_0*x_1 + D(1)*x_0*x_7 + D(1)*x_1*x_2 + D(1)*x_2*x_3 + D(1)*x_3*x_4 + D(1)*x_4*x_5 + D(1)*x_5*x_6 + D(1)*x_6*x_7
        D(1)*x_0*x_1*x_2 + D(1)*x_0*x_1*x_7 + D(1)*x_0*x_6*x_7 + D(1)*x_1*x_2*x_3 + D(1)*x_2*x_3*x_4 + D(1)*x_3*x_4*x_5 + D(1)*x_4*x_5*x_6 + D(1)*x_5*x_6*x_7
        D(1)*x_0*x_1*x_2*x_3 + D(1)*x_0*x_1*x_2*x_7 + D(1)*x_0*x_1*x_6*x_7 + D(1)*x_0*x_5*x_6*x_7 + D(1)*x_1*x_2*x_3*x_4 + D(1)*x_2*x_3*x_4*x_5 + D(1)*x_3*x_4*x_5*x_6 + D(1)*x_4*x_5*x_6*x_7
        D(1)*x_0*x_1*x_2*x_3*x_4 + D(1)*x_0*x_1*x_2*x_3*x_7 + D(1)*x_0*x_1*x_2*x_6*x_7 + D(1)*x_0*x_1*x_5*x_6*x_7 + D(1)*x_0*x_4*x_5*x_6*x_7 + D(1)*x_1*x_2*x_3*x_4*x_5 + D(1)*x_2*x_3*x_4*x_5*x_6 + D(1)*x_3*x_4*x_5*x_6*x_7
        D(1)*x_0*x_1*x_2*x_3*x_4*x_5 + D(1)*x_0*x_1*x_2*x_3*x_4*x_7 + D(1)*x_0*x_1*x_2*x_3*x_6*x_7 + D(1)*x_0*x_1*x_2*x_5*x_6*x_7 + D(1)*x_0*x_1*x_4*x_5*x_6*x_7 + D(1)*x_0*x_3*x_4*x_5*x_6*x_7 + D(1)*x_1*x_2*x_3*x_4*x_5*x_6 + D(1)*x_2*x_3*x_4*x_5*x_6*x_7
        D(1)*x_0*x_1*x_2*x_3*x_4*x_5*x_6 + D(1)*x_0*x_1*x_2*x_3*x_4*x_5*x_7 + D(1)*x_0*x_1*x_2*x_3*x_4*x_6*x_7 + D(1)*x_0*x_1*x_2*x_3*x_5*x_6*x_7 + D(1)*x_0*x_1*x_2*x_4*x_5*x_6*x_7 + D(1)*x_0*x_1*x_3*x_4*x_5*x_6*x_7 + D(1)*x_0*x_2*x_3*x_4*x_5*x_6*x_7 + D(1)*x_1*x_2*x_3*x_4*x_5*x_6*x_7
        D(1)*x_0*x_1*x_2*x_3*x_4*x_5*x_6*x_7 + D(2)*1
    ];

    target_sys = [
        C(1)*x_0 + C(1)*x_1 + C(1)*x_2 + C(1)*x_3 + C(1)*x_4 + C(1)*x_5 + C(1)*x_6 + C(1)*x_7
        C(1)*x_0*x_1 + C(1)*x_0*x_7 + C(1)*x_1*x_2 + C(1)*x_2*x_3 + C(1)*x_3*x_4 + C(1)*x_4*x_5 + C(1)*x_5*x_6 + C(1)*x_6*x_7
        C(1)*x_0*x_1*x_2 + C(1)*x_0*x_1*x_7 + C(1)*x_0*x_6*x_7 + C(1)*x_1*x_2*x_3 + C(1)*x_2*x_3*x_4 + C(1)*x_3*x_4*x_5 + C(1)*x_4*x_5*x_6 + C(1)*x_5*x_6*x_7
        C(1)*x_0*x_1*x_2*x_3 + C(1)*x_0*x_1*x_2*x_7 + C(1)*x_0*x_1*x_6*x_7 + C(1)*x_0*x_5*x_6*x_7 + C(1)*x_1*x_2*x_3*x_4 + C(1)*x_2*x_3*x_4*x_5 + C(1)*x_3*x_4*x_5*x_6 + C(1)*x_4*x_5*x_6*x_7
        C(1)*x_0*x_1*x_2*x_3*x_4 + C(1)*x_0*x_1*x_2*x_3*x_7 + C(1)*x_0*x_1*x_2*x_6*x_7 + C(1)*x_0*x_1*x_5*x_6*x_7 + C(1)*x_0*x_4*x_5*x_6*x_7 + C(1)*x_1*x_2*x_3*x_4*x_5 + C(1)*x_2*x_3*x_4*x_5*x_6 + C(1)*x_3*x_4*x_5*x_6*x_7
        C(1)*x_0*x_1*x_2*x_3*x_4*x_5 + C(1)*x_0*x_1*x_2*x_3*x_4*x_7 + C(1)*x_0*x_1*x_2*x_3*x_6*x_7 + C(1)*x_0*x_1*x_2*x_5*x_6*x_7 + C(1)*x_0*x_1*x_4*x_5*x_6*x_7 + C(1)*x_0*x_3*x_4*x_5*x_6*x_7 + C(1)*x_1*x_2*x_3*x_4*x_5*x_6 + C(1)*x_2*x_3*x_4*x_5*x_6*x_7
        C(1)*x_0*x_1*x_2*x_3*x_4*x_5*x_6 + C(1)*x_0*x_1*x_2*x_3*x_4*x_5*x_7 + C(1)*x_0*x_1*x_2*x_3*x_4*x_6*x_7 + C(1)*x_0*x_1*x_2*x_3*x_5*x_6*x_7 + C(1)*x_0*x_1*x_2*x_4*x_5*x_6*x_7 + C(1)*x_0*x_1*x_3*x_4*x_5*x_6*x_7 + C(1)*x_0*x_2*x_3*x_4*x_5*x_6*x_7 + C(1)*x_1*x_2*x_3*x_4*x_5*x_6*x_7
        C(1)*x_0*x_1*x_2*x_3*x_4*x_5*x_6*x_7 + C(2)*1
    ];

    % -- construct the hpmotopy via the scalar variable t --
    Homotopy = start_sys * (1-t) + target_sys * t;

    % -- compute the Jacobian matrix --
    J = jacobian(Homotopy, [x_0 x_1 x_2 x_3 x_4 x_5 x_6 x_7 t]);
end