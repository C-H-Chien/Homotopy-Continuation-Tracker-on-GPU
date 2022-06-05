function [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_alea6()
    % -- define the system -
    numOfVars = 6;
    numOfCoeff = 29;
    syms x_0 x_1 x_2 x_3 x_4 x_5
    X = [x_0 x_1 x_2 x_3 x_4 x_5];
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
        D(1)*x_0^2*x_3 + D(2)*x_1*x_3*x_4 + D(3)*x_1*x_3*x_5 + D(4)*x_3*x_5 + D(5)*x_4*x_5
        D(6)*x_0*x_1*x_5 + D(7)*x_1^2*x_4 + D(8)*x_1*x_2*x_4 + D(9)*x_1*x_4^2 + D(10)*x_2^2 + D(11)*x_3*x_4*x_5
        D(12)*x_0^2*x_3 + D(13)*x_0^2 + D(14)*x_0*x_3*x_5 + D(13)*x_1^2*x_4 + D(3)*x_1*x_3^2 + D(15)*x_5^3
        D(16)*x_0*x_3^2 + D(14)*x_1*x_3 + D(17)*x_1*x_4^2 + D(18)*x_2*x_3 + D(19)*x_2*x_4*x_5 + D(20)*x_3*x_4^2
        D(21)*x_0^2*x_2 + D(22)*x_0*x_1*x_2 + D(23)*x_0*x_3*x_4 + D(24)*x_0*x_3 + D(25)*x_1^2*x_4
        D(26)*x_0*x_2*x_3 + D(19)*x_2^2*x_3 + D(27)*x_2^2*x_5 + D(28)*x_2 + D(29)*x_3^3 + D(25)*x_4
    ];

    target_sys = [
        C(1)*x_0^2*x_3 + C(2)*x_1*x_3*x_4 + C(3)*x_1*x_3*x_5 + C(4)*x_3*x_5 + C(5)*x_4*x_5
        C(6)*x_0*x_1*x_5 + C(7)*x_1^2*x_4 + C(8)*x_1*x_2*x_4 + C(9)*x_1*x_4^2 + C(10)*x_2^2 + C(11)*x_3*x_4*x_5
        C(12)*x_0^2*x_3 + C(13)*x_0^2 + C(14)*x_0*x_3*x_5 + C(13)*x_1^2*x_4 + C(3)*x_1*x_3^2 + C(15)*x_5^3
        C(16)*x_0*x_3^2 + C(14)*x_1*x_3 + C(17)*x_1*x_4^2 + C(18)*x_2*x_3 + C(19)*x_2*x_4*x_5 + C(20)*x_3*x_4^2
        C(21)*x_0^2*x_2 + C(22)*x_0*x_1*x_2 + C(23)*x_0*x_3*x_4 + C(24)*x_0*x_3 + C(25)*x_1^2*x_4
        C(26)*x_0*x_2*x_3 + C(19)*x_2^2*x_3 + C(27)*x_2^2*x_5 + C(28)*x_2 + C(29)*x_3^3 + C(25)*x_4
    ];

    % -- construct the hpmotopy via the scalar variable t --
    Homotopy = start_sys * (1-t) + target_sys * t;

    % -- compute the Jacobian matrix --
    J = jacobian(Homotopy, [x_0 x_1 x_2 x_3 x_4 x_5 t]);
end