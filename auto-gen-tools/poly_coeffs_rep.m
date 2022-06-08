% -- transform parametric polynomials into coefficient polynomials --
clear;
format long;

% -- define output file-write directories --
problem = 'alea6/';
fileFolder = '/home/chchien/BrownU/research/22-issac/';
category = 'benchmark-problems/';

% -- output file names --
outputFileName_targetCoeffs = 'target_coeffs.txt';
outputFileName_repCoeffs = 'rep_coeffs.txt';
outputFileName_repProblem = 'rep_problem.txt';
fullOutputFileName_targetCoeffs = fullfile(fileFolder, category, problem, outputFileName_targetCoeffs);
fullOutputFileName_repCoeffs = fullfile(fileFolder, category, problem, outputFileName_repCoeffs);
fullOutputFileName_repProblem = fullfile(fileFolder, category, problem, outputFileName_repProblem);
outputFileWr_targetCoeffs = fopen(fullOutputFileName_targetCoeffs, 'w');
outputFileWr_repCoeffs = fopen(fullOutputFileName_repCoeffs, 'w');
outputFileWr_repProblem = fopen(fullOutputFileName_repProblem, 'w');

% -- define systems --
if strcmp(problem, 'katsura6/')
    [f, numOfVars] = sys_katsura6();
elseif strcmp(problem, 'katsura7/')
    [f, numOfVars] = sys_katsura7();
elseif strcmp(problem, 'katsura8/')
    [f, numOfVars] = sys_katsura8();
elseif strcmp(problem, 'katsura9/')
    [f, numOfVars] = sys_katsura9();
elseif strcmp(problem, 'katsura10/')
    [f, numOfVars] = sys_katsura10();
elseif strcmp(problem, 'katsura11/')
    [f, numOfVars] = sys_katsura11();
elseif strcmp(problem, 'katsura12/')
    [f, numOfVars] = sys_katsura12();
elseif strcmp(problem, 'katsura13/')
    [f, numOfVars] = sys_katsura13();
elseif strcmp(problem, 'katsura14/')
    [f, numOfVars] = sys_katsura14();
elseif strcmp(problem, 'katsura15/')
    [f, numOfVars] = sys_katsura15();
elseif strcmp(problem, 'katsura20/')
    [f, numOfVars] = sys_katsura20();
elseif strcmp(problem, 'katsura23/')
    [f, numOfVars] = sys_katsura23();
elseif strcmp(problem, 'cyclic7/')
    [f, numOfVars] = sys_cyclic7();
elseif strcmp(problem, 'cyclic8/')
    [f, numOfVars] = sys_cyclic8();
elseif strcmp(problem, 'cyclic9/')
    [f, numOfVars] = sys_cyclic9();
elseif strcmp(problem, 'cyclic10/')
    [f, numOfVars] = sys_cyclic10();
elseif strcmp(problem, 'cyclic11/')
    [f, numOfVars] = sys_cyclic11();
elseif strcmp(problem, 'alea6/')
    [f, numOfVars] = sys_alea6();
elseif strcmp(problem, 'eco12/')
    [f, numOfVars] = sys_eco12();
elseif strcmp(problem, 'd1/')
    [f, numOfVars] = sys_d1();
elseif strcmp(problem, 'game6two/')
    [f, numOfVars] = sys_game6two();
elseif strcmp(problem, 'game7two/')
    [f, numOfVars] = sys_game7two();
elseif strcmp(problem, 'pole28sys/')
    [f, numOfVars] = sys_pole28sys();
elseif strcmp(problem, 'alea6-extend/')
    [f, numOfVars] = sys_alea6_extend();
end

% -- create symbolic symbols of variables --
for i = 1:numOfVars
    str_vars = 'x';
    cat_str_vars = strcat(str_vars, num2str(i));
    X(i) = str2sym(cat_str_vars);
end

% -- find the variables and coefficients of the system --
coefficient_stack = '';
collect_sys = strings(numOfVars, 1);
for p = 1:numOfVars
    [coefficient, variable] = coeffs(f(p), X);
    
    for i = 1:size(variable, 2)
        consistent_digits_str = sprintf('%.6f',coefficient(i));
        collect_sys(p,1) = strcat(collect_sys(p,1), '(', consistent_digits_str, ')*', char(variable(i)));
        if i < size(variable, 2)
            collect_sys(p,1) = strcat(collect_sys(p,1), {' '}, '+', {' '});
        end
        coefficient_stack = [coefficient_stack string(consistent_digits_str)];
    end
end

coefficient_stack(1) = [];
% -- construct unique coefficients in the coefficient_stack array --
[uniq_coefficients, uniqueIdx] = unique(coefficient_stack);
uniqueIdx = sort(uniqueIdx, 'ascend');
p2c_stack = coefficient_stack(uniqueIdx);

% -- write independent parameter-coefficient correspondences to a file --
as_stack = strings(size(p2c_stack, 2), 1);
p2c_stack_str = strings(size(p2c_stack, 2), 1);
for i = 1:size(p2c_stack, 2)
    as_str = strcat('a(', num2str(i), ')');
    as_stack(i, 1) = as_str;
    p2c_lhs = strcat(as_str, {' '}, '=', {' '});
    
    consistent_digits_str = sprintf('%.6f',p2c_stack(i));
    p2c_rhs = strcat('(', consistent_digits_str, ')');
    cplx_coeffs = strcat(string(double(p2c_stack(i))), '\t', '0.0', '\n');
    fprintf(outputFileWr_targetCoeffs, cplx_coeffs);
    
    p2c_stack_str(i, 1) = p2c_rhs;
    full_p2c = strcat(p2c_lhs, p2c_rhs);
    fprintf(outputFileWr_repCoeffs, string(full_p2c));
    fprintf(outputFileWr_repCoeffs, '\n');
end

% -- replace the coefficients of the collect_sys by p2c correspondences --
fprintf('replacing coefficients expressions ');
for p = 1:numOfVars
    for i = 1:size(p2c_stack_str, 1)
        %if ~strcmp(string(p2c_stack(i)), '1')
            if contains(collect_sys(p,1), p2c_stack_str(i, 1))
                collect_sys(p,1) = strrep(collect_sys(p,1), p2c_stack_str(i, 1), as_stack(i, 1));
            end
        %end
    end
    fprintf('. ')
    fprintf(outputFileWr_repProblem, collect_sys(p,1));
    fprintf(outputFileWr_repProblem, '\n');
end
fprintf('\n');

% -- extract all possible monomials --
collect_sys_monomial = [];
for p = 1:numOfVars
    split_f = strsplit(collect_sys(p,1), '+');
    for i = 1:size(split_f, 2)
        a_indx = extractBetween(split_f(1,i), 'a', '*');
        str_as = strcat('a', a_indx, '*');
        split_f_w_monomial = strrep(split_f(1,i), str_as, '');
        collect_sys_monomial = [collect_sys_monomial; split_f_w_monomial];
    end
end
Z = cellstr(collect_sys_monomial);
ux = unique(Z);

% -- generate system matrix for running on phcpack --
% -- initialize one row of a system matrix for phcpack --
sys_mat = strings(numOfVars+1, 1);

% -- close the files --
fclose(outputFileWr_repCoeffs);
fclose(outputFileWr_repProblem);
fclose(outputFileWr_targetCoeffs);
