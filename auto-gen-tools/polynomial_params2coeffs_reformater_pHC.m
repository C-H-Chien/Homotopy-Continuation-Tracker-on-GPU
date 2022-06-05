% -- transform parametric polynomials into coefficient polynomials --
clear;
format long;

% -- define directories --
% -- top priority --
problem = '3view_unknownf_pHC/';
category = 'computer-vision-problems/';
fileFolder = '/home/chchien/BrownU/research/22-issac/';
outputFileName_coeffs = 'target_coeffs.txt';
fullOutputFileName_coeffs = fullfile(fileFolder, category, problem, outputFileName_coeffs);
outputFileWr_coeffs = fopen(fullOutputFileName_coeffs, 'w');

make_coefficientHC = 0;

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
fprintf("Building the polynomials ...\n");
if strcmp(problem, 'chicago/')
    [f, numOfVars] = sys_chicago();
elseif strcmp(problem, '4vTrg/')
    [f, numOfVars] = sys_4vTrg();
elseif strcmp(problem, '3vTrg/')
    [f, numOfVars] = sys_3vTrg();
elseif strcmp(problem, '3vTrg_relax/')
    [f, numOfVars] = sys_3vTrg_relax();
elseif strcmp(problem, '5pt_rel_pose_w_depth_recon/')
    [f, numOfVars] = sys_5pt_rel_pose_w_depth_recon();
elseif strcmp(problem, 'PnP_wo_principal_point/')
    [f, numOfVars] = sys_PnP_wo_principal_point();
elseif strcmp(problem, 'optimalPnP_w_quaternion/')
    [f, numOfVars] = sys_optimalPnP_w_quaternion();
elseif strcmp(problem, '3pt_rel_pose_w_homo_constraint/')
    [f, numOfVars] = sys_3pt_rel_pose_w_homo_constraint();
elseif strcmp(problem, 'r6p/')
    [f, numOfVars] = sys_r6p();
elseif strcmp(problem, 'refractive_p5p/')
    [f, numOfVars] = sys_refractive_p5p();
elseif strcmp(problem, 'refractive_p6p/')
    [f, numOfVars] = sys_refractive_p6p();
elseif strcmp(problem, '3view_unknownf_pHC/')
    [f, numOfVars] = sys_3view_unknownf_pHC();
end


for i = 1:numOfVars
    str_vars = 'x';
    cat_str_vars = strcat(str_vars, num2str(i));
    X(i) = str2sym(cat_str_vars);
end

% -- find the variables and coefficients of the system --
fprintf("Staking all possible coefficients ");
coefficient_stack = '';
collect_sys = strings(numOfVars, 1);
for p = 1:numOfVars
    [coefficient, variable] = coeffs(f(p), X);
    
    for i = 1:size(variable, 2)
        
        collect_sys(p,1) = strcat(collect_sys(p,1), '(', string(coefficient(i)), ')*', string(variable(i)));
        if i < size(variable, 2)
            collect_sys(p,1) = strcat(collect_sys(p,1), {' '}, '+', {' '});
        end
        coefficient_stack = [coefficient_stack string(coefficient(i))];
    end
    
    str_collect_sys = string(collect_sys(p,1));
    fprintf(outputFileWr_problem, str_collect_sys);
    fprintf(outputFileWr_problem, '\n');
    
    fprintf('. ');
end
fprintf('\n');

% -- TODO: here we need to replace duplicate coefficients with opposit signs --


coefficient_stack(1) = [];
% -- construct unique coefficients in the coefficient_stack array --
[uniq_coefficients, uniqueIdx] = unique(coefficient_stack);
uniqueIdx = sort(uniqueIdx, 'ascend');
p2c_stack = coefficient_stack(uniqueIdx);

% -- write independent parameter-coefficient correspondences to a file --
fprintf("Complete writing individual coefficients to a file ");
as_stack = strings(size(p2c_stack, 2), 1);
p2c_stack_str = strings(size(p2c_stack, 2), 1);
for i = 1:size(p2c_stack, 2)
    %as_str = strcat('a(', num2str(i), ')');
    %as_str = strcat('a', num2str(i));
    
    if make_coefficientHC
        as_str = strcat('h_coeffs[', num2str(i-1), ']');
    else
        as_str = strcat('a', num2str(i));
    end
    
    as_stack(i, 1) = as_str;
    p2c_lhs = strcat(as_str, {' '}, '=', {' '});
    
    consistent_digits_str = p2c_stack(i);
    p2c_rhs = strcat('(', consistent_digits_str, ')');
    
    if make_coefficientHC
        p2c_rhs = strcat(p2c_rhs, ';');
    end
    
    p2c_stack_str(i, 1) = p2c_rhs;
    full_p2c = strcat(p2c_lhs, p2c_rhs);
    fprintf(outputFileWr_p2c, string(full_p2c));
    fprintf(outputFileWr_p2c, '\n');
    
    fprintf('. ');
end
fprintf('\n');


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
    fprintf('. ');
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

% -- close the files --
fclose(outputFileWr_p2c);
fclose(outputFileWr_problem);
fclose(outputFileWr_rep_problem);
fclose(outputFileWr_coeffs);
