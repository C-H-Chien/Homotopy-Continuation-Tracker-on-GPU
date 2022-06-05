% -- transform parametric polynomials into coefficient polynomials --
clear;
format long;

% -- define directories --
% -- top priority --
problem = 'alea6-extend/';
category = 'benchmark-problems/';
fileFolder = '/home/chchien/BrownU/research/22-issac/';
outputFileName_coeffs = 'target_coeffs.txt';
fullOutputFileName_coeffs = fullfile(fileFolder, category, problem, outputFileName_coeffs);
outputFileWr_coeffs = fopen(fullOutputFileName_coeffs, 'w');

% -- generate phcpack polynomial format --
target_folder = 'phcpack/';
outputFileName_phcpack_target_sys = 'phcpack-target-sys-formulations.txt';
outputFileName_phcpack_start_sys = 'phcpack-start-sys-formulations.txt';
fullOutputFileName_phcpck_target_sys = fullfile(fileFolder, category, problem, target_folder, outputFileName_phcpack_target_sys);
fullOutputFileName_phcpck_start_sys = fullfile(fileFolder, category, problem, target_folder, outputFileName_phcpack_start_sys);
outputFileWr_phcpack_target_sys = fopen(fullOutputFileName_phcpck_target_sys, 'w');
outputFileWr_phcpack_start_sys = fopen(fullOutputFileName_phcpck_start_sys, 'w');

outputFileName_p2c = 'params2coeffs.txt';
outputFileName_problem = 'collected_problem.txt';
outputFileName_coefficient_problem = 'rep_problem.txt';
fullOutputFileName_p2c = fullfile(fileFolder, category, problem, outputFileName_p2c);
fullOutputFileName_problem = fullfile(fileFolder, category, problem, outputFileName_problem);
fullOutputFileName_rep_problem = fullfile(fileFolder, category, problem, outputFileName_coefficient_problem);
outputFileWr_p2c = fopen(fullOutputFileName_p2c, 'w');
outputFileWr_problem = fopen(fullOutputFileName_problem, 'w');
outputFileWr_rep_problem = fopen(fullOutputFileName_rep_problem, 'w');

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
    
    str_collect_sys = string(collect_sys(p,1));
    fprintf(outputFileWr_problem, str_collect_sys);
    fprintf(outputFileWr_problem, '\n');
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
    fprintf(outputFileWr_coeffs, cplx_coeffs);
    
    p2c_stack_str(i, 1) = p2c_rhs;
    full_p2c = strcat(p2c_lhs, p2c_rhs);
    fprintf(outputFileWr_p2c, string(full_p2c));
    fprintf(outputFileWr_p2c, '\n');
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
    fprintf(outputFileWr_rep_problem, collect_sys(p,1));
    fprintf(outputFileWr_rep_problem, '\n');
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

% -- phcpack system matrix for target system --
for p = 1:numOfVars
    % -- 1) find the variables and coefficients of the target system --
    [coefficient, variable] = coeffs(f(p), X);
    
    % -- split terms of replaced system --
    split_f = strsplit(collect_sys(p,1), '+');
    
    % -- 2) loop over all terms --
    for i = 1:size(variable, 2)
        % -- 2-1) initialize all elements in a row to 0 --
        for k = 1:numOfVars+1
            sys_mat(k,1) = '0';
        end
        
        % -- 2-2) store the coefficient first --
        str_coeff = string(coefficient(i));
        % -- check if the coefficient in the string format is a fraction --
        if contains(str_coeff, '/')
            numer = double(extractBefore(str_coeff, '/'));
            denom = double(extractAfter(str_coeff, '/'));
            converted_number = numer / denom;
            sys_mat(1,1) = num2str(converted_number, '%.6f');
        else
            sys_mat(1,1) = str_coeff;
        end
        
%         frac = sscanf(str_coeff,'%d/%d');
%         if size(frac, 1) > 1
%             % -- if it is a fractional number, convert it to a floating
%             % point number in the string format --
%             converted_number = frac(1) / frac(2);
%             sys_mat(1,1) = num2str(converted_number, '%.6f');
%         else
%             sys_mat(1,1) = str_coeff;
%         end
        
        % -- 2-3) split the unknowns by multiplier * -- 
        unknowns = strsplit(string(variable(i)), '*');
        
        % -- 2-4) loop over all unknowns --
        for j = 1:size(unknowns, 2)
            % -- check the power of each unknown --
            if contains(unknowns(j), '^')
                degree = extractAfter(unknowns(j), '^');
                x_idx = extractBetween(unknowns(j), 'x', '^');
                sys_mat(str2double(x_idx)+1, 1) = degree;
            elseif contains(unknowns(j), 'x')
                x_idx = extractAfter(unknowns(j), 'x');
                sys_mat(str2double(x_idx)+1, 1) = '1';
            else
                continue;
            end
        end
        
        % -- 2-5) write the target results to the file --
        for k = 1:numOfVars+1
            fprintf(outputFileWr_phcpack_target_sys, sys_mat(k,1));
            fprintf(outputFileWr_phcpack_target_sys, " ");
        end
        fprintf(outputFileWr_phcpack_target_sys, ';\n');
        
        % -- 2-6) for start coefficients, aka a's --
        a_indx = extractBetween(split_f(1,i), 'a', '*');
        a_indx = strrep(a_indx, '(', '');
        a_indx = strrep(a_indx, ')', '');
        str_as = strcat('a', a_indx);
        sys_mat(1,1) = str_as;
        
        % -- 2-7) write the start results to the file --
        for k = 1:numOfVars+1
            fprintf(outputFileWr_phcpack_start_sys, sys_mat(k,1));
            fprintf(outputFileWr_phcpack_start_sys, " ");
        end
        fprintf(outputFileWr_phcpack_start_sys, ';\n');
    end
    
    % -- 3) write the sperator to the file --
    for k = 1:numOfVars+1
        fprintf(outputFileWr_phcpack_target_sys, '0');
        fprintf(outputFileWr_phcpack_target_sys, " ");
        fprintf(outputFileWr_phcpack_start_sys, '0');
        fprintf(outputFileWr_phcpack_start_sys, " ");
    end
    fprintf(outputFileWr_phcpack_target_sys, ';\n');
    fprintf(outputFileWr_phcpack_start_sys, ';\n');
end

% -- close the files --
fclose(outputFileWr_p2c);
fclose(outputFileWr_problem);
fclose(outputFileWr_rep_problem);
fclose(outputFileWr_coeffs);
fclose(outputFileWr_phcpack_target_sys);
fclose(outputFileWr_phcpack_start_sys);
