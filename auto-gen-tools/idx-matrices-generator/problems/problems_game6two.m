function [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_game6two()
    % -- define the system -
    numOfVars = 6;
    numOfCoeff = 191;
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
        D(1)*x_1*x_2*x_3*x_4*x_5 + D(2)*x_1*x_2*x_3*x_4 + D(3)*x_1*x_2*x_3*x_5 + D(4)*x_1*x_2*x_3 + D(5)*x_1*x_2*x_4*x_5 + D(6)*x_1*x_2*x_4 + D(7)*x_1*x_2*x_5 + D(8)*x_1*x_2 + D(9)*x_1*x_3*x_4*x_5 + D(10)*x_1*x_3*x_4 + D(11)*x_1*x_3*x_5 + D(12)*x_1*x_3 + D(13)*x_1*x_4*x_5 + D(14)*x_1*x_4 + D(15)*x_1*x_5 + D(16)*x_1 + D(17)*x_2*x_3*x_4*x_5 + D(18)*x_2*x_3*x_4 + D(19)*x_2*x_3*x_5 + D(20)*x_2*x_3 + D(21)*x_2*x_4*x_5 + D(22)*x_2*x_4 + D(23)*x_2*x_5 + D(24)*x_2 + D(25)*x_3*x_4*x_5 + D(26)*x_3*x_4 + D(27)*x_3*x_5 + D(28)*x_3 + D(29)*x_4*x_5 + D(30)*x_4 + D(31)*x_5 + D(32)*1
        D(33)*x_0*x_2*x_3*x_4*x_5 + D(34)*x_0*x_2*x_3*x_4 + D(35)*x_0*x_2*x_3*x_5 + D(36)*x_0*x_2*x_3 + D(37)*x_0*x_2*x_4*x_5 + D(38)*x_0*x_2*x_4 + D(39)*x_0*x_2*x_5 + D(40)*x_0*x_2 + D(41)*x_0*x_3*x_4*x_5 + D(42)*x_0*x_3*x_4 + D(43)*x_0*x_3*x_5 + D(44)*x_0*x_3 + D(45)*x_0*x_4*x_5 + D(46)*x_0*x_4 + D(47)*x_0*x_5 + D(48)*x_0 + D(49)*x_2*x_3*x_4*x_5 + D(50)*x_2*x_3*x_4 + D(51)*x_2*x_3*x_5 + D(52)*x_2*x_3 + D(53)*x_2*x_4*x_5 + D(54)*x_2*x_4 + D(55)*x_2*x_5 + D(56)*x_2 + D(57)*x_3*x_4*x_5 + D(58)*x_3*x_4 + D(59)*x_3*x_5 + D(60)*x_3 + D(61)*x_4*x_5 + D(62)*x_4 + D(63)*x_5 + D(64)*1
        D(65)*x_0*x_1*x_3*x_4*x_5 + D(66)*x_0*x_1*x_3*x_4 + D(67)*x_0*x_1*x_3*x_5 + D(68)*x_0*x_1*x_3 + D(69)*x_0*x_1*x_4*x_5 + D(70)*x_0*x_1*x_4 + D(71)*x_0*x_1*x_5 + D(72)*x_0*x_1 + D(73)*x_0*x_3*x_4*x_5 + D(74)*x_0*x_3*x_4 + D(75)*x_0*x_3*x_5 + D(76)*x_0*x_3 + D(77)*x_0*x_4*x_5 + D(78)*x_0*x_4 + D(79)*x_0*x_5 + D(80)*x_0 + D(81)*x_1*x_3*x_4*x_5 + D(82)*x_1*x_3*x_4 + D(83)*x_1*x_3*x_5 + D(84)*x_1*x_3 + D(85)*x_1*x_4*x_5 + D(86)*x_1*x_4 + D(87)*x_1*x_5 + D(88)*x_1 + D(89)*x_3*x_4*x_5 + D(90)*x_3*x_4 + D(91)*x_3*x_5 + D(92)*x_3 + D(93)*x_4*x_5 + D(94)*x_4 + D(95)*x_5 + D(96)*1
        D(97)*x_0*x_1*x_2*x_4*x_5 + D(98)*x_0*x_1*x_2*x_4 + D(99)*x_0*x_1*x_2*x_5 + D(100)*x_0*x_1*x_2 + D(101)*x_0*x_1*x_4*x_5 + D(102)*x_0*x_1*x_4 + D(103)*x_0*x_1*x_5 + D(104)*x_0*x_1 + D(105)*x_0*x_2*x_4*x_5 + D(106)*x_0*x_2*x_4 + D(107)*x_0*x_2*x_5 + D(108)*x_0*x_2 + D(109)*x_0*x_4*x_5 + D(110)*x_0*x_4 + D(111)*x_0*x_5 + D(112)*x_0 + D(113)*x_1*x_2*x_4*x_5 + D(114)*x_1*x_2*x_4 + D(115)*x_1*x_2*x_5 + D(116)*x_1*x_2 + D(117)*x_1*x_4*x_5 + D(118)*x_1*x_4 + D(119)*x_1*x_5 + D(120)*x_1 + D(121)*x_2*x_4*x_5 + D(122)*x_2*x_4 + D(123)*x_2*x_5 + D(124)*x_2 + D(125)*x_4*x_5 + D(126)*x_4 + D(127)*x_5 + D(128)*1
        D(129)*x_0*x_1*x_2*x_3*x_5 + D(130)*x_0*x_1*x_2*x_3 + D(131)*x_0*x_1*x_2*x_5 + D(132)*x_0*x_1*x_2 + D(133)*x_0*x_1*x_3*x_5 + D(134)*x_0*x_1*x_3 + D(135)*x_0*x_1*x_5 + D(136)*x_0*x_1 + D(137)*x_0*x_2*x_3*x_5 + D(138)*x_0*x_2*x_3 + D(139)*x_0*x_2*x_5 + D(140)*x_0*x_2 + D(141)*x_0*x_3*x_5 + D(142)*x_0*x_3 + D(143)*x_0*x_5 + D(144)*x_0 + D(145)*x_1*x_2*x_3*x_5 + D(146)*x_1*x_2*x_3 + D(147)*x_1*x_2*x_5 + D(148)*x_1*x_2 + D(149)*x_1*x_3*x_5 + D(150)*x_1*x_3 + D(151)*x_1*x_5 + D(152)*x_1 + D(153)*x_2*x_3*x_5 + D(154)*x_2*x_3 + D(155)*x_2*x_5 + D(156)*x_2 + D(157)*x_3*x_5 + D(158)*x_3 + D(159)*x_5 + D(160)*1
        D(161)*x_0*x_1*x_2*x_3*x_4 + D(162)*x_0*x_1*x_2*x_3 + D(163)*x_0*x_1*x_2*x_4 + D(164)*x_0*x_1*x_2 + D(165)*x_0*x_1*x_3*x_4 + D(166)*x_0*x_1*x_3 + D(167)*x_0*x_1*x_4 + D(168)*x_0*x_1 + D(169)*x_0*x_2*x_3*x_4 + D(170)*x_0*x_2*x_3 + D(171)*x_0*x_2*x_4 + D(172)*x_0*x_2 + D(173)*x_0*x_3*x_4 + D(174)*x_0*x_3 + D(175)*x_0*x_4 + D(176)*x_0 + D(177)*x_1*x_2*x_3*x_4 + D(178)*x_1*x_2*x_3 + D(179)*x_1*x_2*x_4 + D(180)*x_1*x_2 + D(79)*x_1*x_3*x_4 + D(181)*x_1*x_3 + D(182)*x_1*x_4 + D(183)*x_1 + D(184)*x_2*x_3*x_4 + D(185)*x_2*x_3 + D(186)*x_2*x_4 + D(187)*x_2 + D(188)*x_3*x_4 + D(189)*x_3 + D(190)*x_4 + D(191)*1
    ];

    target_sys = [
        C(1)*x_1*x_2*x_3*x_4*x_5 + C(2)*x_1*x_2*x_3*x_4 + C(3)*x_1*x_2*x_3*x_5 + C(4)*x_1*x_2*x_3 + C(5)*x_1*x_2*x_4*x_5 + C(6)*x_1*x_2*x_4 + C(7)*x_1*x_2*x_5 + C(8)*x_1*x_2 + C(9)*x_1*x_3*x_4*x_5 + C(10)*x_1*x_3*x_4 + C(11)*x_1*x_3*x_5 + C(12)*x_1*x_3 + C(13)*x_1*x_4*x_5 + C(14)*x_1*x_4 + C(15)*x_1*x_5 + C(16)*x_1 + C(17)*x_2*x_3*x_4*x_5 + C(18)*x_2*x_3*x_4 + C(19)*x_2*x_3*x_5 + C(20)*x_2*x_3 + C(21)*x_2*x_4*x_5 + C(22)*x_2*x_4 + C(23)*x_2*x_5 + C(24)*x_2 + C(25)*x_3*x_4*x_5 + C(26)*x_3*x_4 + C(27)*x_3*x_5 + C(28)*x_3 + C(29)*x_4*x_5 + C(30)*x_4 + C(31)*x_5 + C(32)*1
        C(33)*x_0*x_2*x_3*x_4*x_5 + C(34)*x_0*x_2*x_3*x_4 + C(35)*x_0*x_2*x_3*x_5 + C(36)*x_0*x_2*x_3 + C(37)*x_0*x_2*x_4*x_5 + C(38)*x_0*x_2*x_4 + C(39)*x_0*x_2*x_5 + C(40)*x_0*x_2 + C(41)*x_0*x_3*x_4*x_5 + C(42)*x_0*x_3*x_4 + C(43)*x_0*x_3*x_5 + C(44)*x_0*x_3 + C(45)*x_0*x_4*x_5 + C(46)*x_0*x_4 + C(47)*x_0*x_5 + C(48)*x_0 + C(49)*x_2*x_3*x_4*x_5 + C(50)*x_2*x_3*x_4 + C(51)*x_2*x_3*x_5 + C(52)*x_2*x_3 + C(53)*x_2*x_4*x_5 + C(54)*x_2*x_4 + C(55)*x_2*x_5 + C(56)*x_2 + C(57)*x_3*x_4*x_5 + C(58)*x_3*x_4 + C(59)*x_3*x_5 + C(60)*x_3 + C(61)*x_4*x_5 + C(62)*x_4 + C(63)*x_5 + C(64)*1
        C(65)*x_0*x_1*x_3*x_4*x_5 + C(66)*x_0*x_1*x_3*x_4 + C(67)*x_0*x_1*x_3*x_5 + C(68)*x_0*x_1*x_3 + C(69)*x_0*x_1*x_4*x_5 + C(70)*x_0*x_1*x_4 + C(71)*x_0*x_1*x_5 + C(72)*x_0*x_1 + C(73)*x_0*x_3*x_4*x_5 + C(74)*x_0*x_3*x_4 + C(75)*x_0*x_3*x_5 + C(76)*x_0*x_3 + C(77)*x_0*x_4*x_5 + C(78)*x_0*x_4 + C(79)*x_0*x_5 + C(80)*x_0 + C(81)*x_1*x_3*x_4*x_5 + C(82)*x_1*x_3*x_4 + C(83)*x_1*x_3*x_5 + C(84)*x_1*x_3 + C(85)*x_1*x_4*x_5 + C(86)*x_1*x_4 + C(87)*x_1*x_5 + C(88)*x_1 + C(89)*x_3*x_4*x_5 + C(90)*x_3*x_4 + C(91)*x_3*x_5 + C(92)*x_3 + C(93)*x_4*x_5 + C(94)*x_4 + C(95)*x_5 + C(96)*1
        C(97)*x_0*x_1*x_2*x_4*x_5 + C(98)*x_0*x_1*x_2*x_4 + C(99)*x_0*x_1*x_2*x_5 + C(100)*x_0*x_1*x_2 + C(101)*x_0*x_1*x_4*x_5 + C(102)*x_0*x_1*x_4 + C(103)*x_0*x_1*x_5 + C(104)*x_0*x_1 + C(105)*x_0*x_2*x_4*x_5 + C(106)*x_0*x_2*x_4 + C(107)*x_0*x_2*x_5 + C(108)*x_0*x_2 + C(109)*x_0*x_4*x_5 + C(110)*x_0*x_4 + C(111)*x_0*x_5 + C(112)*x_0 + C(113)*x_1*x_2*x_4*x_5 + C(114)*x_1*x_2*x_4 + C(115)*x_1*x_2*x_5 + C(116)*x_1*x_2 + C(117)*x_1*x_4*x_5 + C(118)*x_1*x_4 + C(119)*x_1*x_5 + C(120)*x_1 + C(121)*x_2*x_4*x_5 + C(122)*x_2*x_4 + C(123)*x_2*x_5 + C(124)*x_2 + C(125)*x_4*x_5 + C(126)*x_4 + C(127)*x_5 + C(128)*1
        C(129)*x_0*x_1*x_2*x_3*x_5 + C(130)*x_0*x_1*x_2*x_3 + C(131)*x_0*x_1*x_2*x_5 + C(132)*x_0*x_1*x_2 + C(133)*x_0*x_1*x_3*x_5 + C(134)*x_0*x_1*x_3 + C(135)*x_0*x_1*x_5 + C(136)*x_0*x_1 + C(137)*x_0*x_2*x_3*x_5 + C(138)*x_0*x_2*x_3 + C(139)*x_0*x_2*x_5 + C(140)*x_0*x_2 + C(141)*x_0*x_3*x_5 + C(142)*x_0*x_3 + C(143)*x_0*x_5 + C(144)*x_0 + C(145)*x_1*x_2*x_3*x_5 + C(146)*x_1*x_2*x_3 + C(147)*x_1*x_2*x_5 + C(148)*x_1*x_2 + C(149)*x_1*x_3*x_5 + C(150)*x_1*x_3 + C(151)*x_1*x_5 + C(152)*x_1 + C(153)*x_2*x_3*x_5 + C(154)*x_2*x_3 + C(155)*x_2*x_5 + C(156)*x_2 + C(157)*x_3*x_5 + C(158)*x_3 + C(159)*x_5 + C(160)*1
        C(161)*x_0*x_1*x_2*x_3*x_4 + C(162)*x_0*x_1*x_2*x_3 + C(163)*x_0*x_1*x_2*x_4 + C(164)*x_0*x_1*x_2 + C(165)*x_0*x_1*x_3*x_4 + C(166)*x_0*x_1*x_3 + C(167)*x_0*x_1*x_4 + C(168)*x_0*x_1 + C(169)*x_0*x_2*x_3*x_4 + C(170)*x_0*x_2*x_3 + C(171)*x_0*x_2*x_4 + C(172)*x_0*x_2 + C(173)*x_0*x_3*x_4 + C(174)*x_0*x_3 + C(175)*x_0*x_4 + C(176)*x_0 + C(177)*x_1*x_2*x_3*x_4 + C(178)*x_1*x_2*x_3 + C(179)*x_1*x_2*x_4 + C(180)*x_1*x_2 + C(79)*x_1*x_3*x_4 + C(181)*x_1*x_3 + C(182)*x_1*x_4 + C(183)*x_1 + C(184)*x_2*x_3*x_4 + C(185)*x_2*x_3 + C(186)*x_2*x_4 + C(187)*x_2 + C(188)*x_3*x_4 + C(189)*x_3 + C(190)*x_4 + C(191)*1
    ];

    % -- construct the hpmotopy via the scalar variable t --
    Homotopy = start_sys * (1-t) + target_sys * t;

    % -- compute the Jacobian matrix --
    J = jacobian(Homotopy, [x_0 x_1 x_2 x_3 x_4 x_5 t]);
end