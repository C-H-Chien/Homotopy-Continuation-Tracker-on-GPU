function [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_r6p()
    % -- define the system -
    numOfVars = 6;
    numOfCoeff = 96;
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
        D(1)*x_0*x_3 + D(2)*x_0*x_4 + D(3)*x_0*x_5 + D(4)*x_0 + D(5)*x_1*x_3 + D(6)*x_1*x_4 + D(7)*x_1*x_5 + D(8)*x_1 + D(9)*x_2*x_3 + D(10)*x_2*x_4 + D(11)*x_2*x_5 + D(12)*x_2 + D(13)*x_3 + D(14)*x_4 + D(15)*x_5 + D(16)
        D(17)*x_0*x_3 + D(18)*x_0*x_4 + D(19)*x_0*x_5 + D(20)*x_0 + D(21)*x_1*x_3 + D(22)*x_1*x_4 + D(23)*x_1*x_5 + D(24)*x_1 + D(25)*x_2*x_3 + D(26)*x_2*x_4 + D(27)*x_2*x_5 + D(28)*x_2 + D(29)*x_3 + D(30)*x_4 + D(31)*x_5 + D(32)
        D(33)*x_0*x_3 + D(34)*x_0*x_4 + D(35)*x_0*x_5 + D(36)*x_0 + D(37)*x_1*x_3 + D(38)*x_1*x_4 + D(39)*x_1*x_5 + D(40)*x_1 + D(41)*x_2*x_3 + D(42)*x_2*x_4 + D(43)*x_2*x_5 + D(44)*x_2 + D(45)*x_3 + D(46)*x_4 + D(47)*x_5 + D(48)
        D(49)*x_0*x_3 + D(50)*x_0*x_4 + D(51)*x_0*x_5 + D(52)*x_0 + D(53)*x_1*x_3 + D(54)*x_1*x_4 + D(55)*x_1*x_5 + D(56)*x_1 + D(57)*x_2*x_3 + D(58)*x_2*x_4 + D(59)*x_2*x_5 + D(60)*x_2 + D(61)*x_3 + D(62)*x_4 + D(63)*x_5 + D(64)
        D(65)*x_0*x_3 + D(66)*x_0*x_4 + D(67)*x_0*x_5 + D(68)*x_0 + D(69)*x_1*x_3 + D(70)*x_1*x_4 + D(71)*x_1*x_5 + D(72)*x_1 + D(73)*x_2*x_3 + D(74)*x_2*x_4 + D(75)*x_2*x_5 + D(76)*x_2 + D(77)*x_3 + D(78)*x_4 + D(79)*x_5 + D(80)
        D(81)*x_0*x_3 + D(82)*x_0*x_4 + D(83)*x_0*x_5 + D(84)*x_0 + D(85)*x_1*x_3 + D(86)*x_1*x_4 + D(87)*x_1*x_5 + D(88)*x_1 + D(89)*x_2*x_3 + D(90)*x_2*x_4 + D(91)*x_2*x_5 + D(92)*x_2 + D(93)*x_3 + D(94)*x_4 + D(95)*x_5 + D(96)
    ];

    target_sys = [
        C(1)*x_0*x_3 + C(2)*x_0*x_4 + C(3)*x_0*x_5 + C(4)*x_0 + C(5)*x_1*x_3 + C(6)*x_1*x_4 + C(7)*x_1*x_5 + C(8)*x_1 + C(9)*x_2*x_3 + C(10)*x_2*x_4 + C(11)*x_2*x_5 + C(12)*x_2 + C(13)*x_3 + C(14)*x_4 + C(15)*x_5 + C(16)
        C(17)*x_0*x_3 + C(18)*x_0*x_4 + C(19)*x_0*x_5 + C(20)*x_0 + C(21)*x_1*x_3 + C(22)*x_1*x_4 + C(23)*x_1*x_5 + C(24)*x_1 + C(25)*x_2*x_3 + C(26)*x_2*x_4 + C(27)*x_2*x_5 + C(28)*x_2 + C(29)*x_3 + C(30)*x_4 + C(31)*x_5 + C(32)
        C(33)*x_0*x_3 + C(34)*x_0*x_4 + C(35)*x_0*x_5 + C(36)*x_0 + C(37)*x_1*x_3 + C(38)*x_1*x_4 + C(39)*x_1*x_5 + C(40)*x_1 + C(41)*x_2*x_3 + C(42)*x_2*x_4 + C(43)*x_2*x_5 + C(44)*x_2 + C(45)*x_3 + C(46)*x_4 + C(47)*x_5 + C(48)
        C(49)*x_0*x_3 + C(50)*x_0*x_4 + C(51)*x_0*x_5 + C(52)*x_0 + C(53)*x_1*x_3 + C(54)*x_1*x_4 + C(55)*x_1*x_5 + C(56)*x_1 + C(57)*x_2*x_3 + C(58)*x_2*x_4 + C(59)*x_2*x_5 + C(60)*x_2 + C(61)*x_3 + C(62)*x_4 + C(63)*x_5 + C(64)
        C(65)*x_0*x_3 + C(66)*x_0*x_4 + C(67)*x_0*x_5 + C(68)*x_0 + C(69)*x_1*x_3 + C(70)*x_1*x_4 + C(71)*x_1*x_5 + C(72)*x_1 + C(73)*x_2*x_3 + C(74)*x_2*x_4 + C(75)*x_2*x_5 + C(76)*x_2 + C(77)*x_3 + C(78)*x_4 + C(79)*x_5 + C(80)
        C(81)*x_0*x_3 + C(82)*x_0*x_4 + C(83)*x_0*x_5 + C(84)*x_0 + C(85)*x_1*x_3 + C(86)*x_1*x_4 + C(87)*x_1*x_5 + C(88)*x_1 + C(89)*x_2*x_3 + C(90)*x_2*x_4 + C(91)*x_2*x_5 + C(92)*x_2 + C(93)*x_3 + C(94)*x_4 + C(95)*x_5 + C(96)
    ];

    % -- construct the hpmotopy via the scalar variable t --
    Homotopy = start_sys * (1-t) + target_sys * t;

    % -- compute the Jacobian matrix --
    J = jacobian(Homotopy, [x_0 x_1 x_2 x_3 x_4 x_5 t]);
end