% -- generate  --
clear;

%sympref('FloatingPointOutput',true);

% -- define directories --
fileFolder = '/home/chchien/BrownU/research/22-issac/computer-vision-problems/';
problem = 'chicago/';
numOfParams = 56;

% -- input file --
% -- 1) parameters to coefficient conversion file --
inputName_p2c = 'params2coeffs.txt';
fullInputFileName_p2c = fullfile(fileFolder, problem, inputName_p2c);
inputFileWr_p2c = fopen(fullInputFileName_p2c, 'r');

% -- output files --
outputFileName_phomotopy = 'params_homotopy_coeffs.txt';
outputFileName_coeffs_of_func_t = 'coeffs_of_func_t.txt';
outputFileName_coeffs_of_func_t_Ht = 'coeffs_of_func_t_Ht.txt';
outputFileName_param_arrays = 'params_arrays.txt';
outputFileNmae_phomotopy_raw = 'params_homotopy_coeffs_raw.txt';
fullOutputFileName_phomotopy = fullfile(fileFolder, problem, outputFileName_phomotopy);
fullOutputFileName_coeffs_of_func_t = fullfile(fileFolder, problem, outputFileName_coeffs_of_func_t);
fullOutputFileName_coeffs_of_func_t_Ht = fullfile(fileFolder, problem, outputFileName_coeffs_of_func_t_Ht);
fullOutputFileName_param_arrays = fullfile(fileFolder, problem, outputFileName_param_arrays);
fullOutputFileName_phomotopy_raw = fullfile(fileFolder, problem, outputFileNmae_phomotopy_raw);
outputFileWr_phomotopy = fopen(fullOutputFileName_phomotopy, 'w');
outputFileWr_coeffs_of_func_t = fopen(fullOutputFileName_coeffs_of_func_t, 'w');
outputFileWr_coeffs_of_func_t_Ht = fopen(fullOutputFileName_coeffs_of_func_t_Ht, 'w');
outputFileWr_param_arrays = fopen(fullOutputFileName_param_arrays, 'w');
outputFileWr_phomotopy_raw = fopen(fullOutputFileName_phomotopy_raw, 'w');

% -- read start and target parameters --
% inputName_start_params = 'start_params.txt';
% inputName_target_params = 'target_params.txt';
% fullInputFileName_start_params = fullfile(fileFolder, problem, inputName_start_params);
% fullInputFileName_target_params = fullfile(fileFolder, problem, inputName_target_params);
% inputFileWr_start_p = fopen(fullInputFileName_start_params, 'r');
% inputFileWr_target_p = fopen(fullInputFileName_target_params, 'r');

ldata = textscan(inputFileWr_p2c, '%s', 'Delimiter', ',', 'whitespace', '{}', 'CollectOutput', true);
line_p2c = ldata{1};
% ldata = textscan(inputFileWr_start_p, '%s', 'Delimiter', ',', 'whitespace', '{}', 'CollectOutput', true  );
% line_sp = ldata{1};
% ldata = textscan(inputFileWr_target_p, '%s', 'Delimiter', ',', 'whitespace', '{}', 'CollectOutput', true  );
% line_tp = ldata{1};

% -- convert start and target string parameters into complex numbers --
% p = complex(size(line_sp, 1), 1); % -- target --
% q = complex(size(line_sp, 1), 1); % -- start --
% function_t = strings(size(line_p2c, 1), 1);
% for i = 1:size(line_sp, 1)
%     % -- start --
%     split_param_str = strsplit(string(line_sp{i,1}), {' '});
%     str_complex = strcat(split_param_str(1,1), '+', split_param_str(1,2), 'i');
%     q(i,1) = str2double(str_complex);
%     
%     % -- target --
%     split_param_str = strsplit(string(line_tp{i,1}), {' '});
%     str_complex = strcat(split_param_str(1,1), '+', split_param_str(1,2), 'i');
%     p(i,1) = str2double(str_complex);
% end

sym_p_homotopy_all = [];
max_t_size = 0;
acc_numOfCoeffs = 0;
fprintf('construct parameter homotopy for each coefficient ');
for i = 1:size(line_p2c, 1)
    % -- 1) extract the RHS and convert it to a string, and remove parantheses --
    bothSides = strsplit(line_p2c{i,1}, '=');
    RHS = bothSides{1,2};
    str_RHS = string(strrep(RHS, {' '}, ''));
    str_RHS = strrep(str_RHS, '(', '');
    str_RHS = strrep(str_RHS, ')', '');
    
    % -- 2) splite the RHS into terms and arithmetic symbols --
    [p_terms, arith_sym] = strsplit(str_RHS, {'+','-'});
    
    % -- 3) loop over all terms and replace it by parameter homotopy --
    full_p_homotopy_str = '';
    for j = 1:size(p_terms, 2)
        % -- 3-1) check whether there are any power in each term --
        if contains(p_terms(1,j), '^')
            c_power_check = extractAfter(p_terms(1,j), '^');
            if contains(c_power_check, '*')
                c_power = extractBetween(p_terms(1,j), '^', '*');
            else
                c_power = c_power_check;
            end
            c_idx_of_power = extractBetween(p_terms(1,j), 'c', '^');
            
            % -- generate the multiplication 
            c_exp_power = '';
            for np = 1:double(c_power)
                c_exp_power = strcat(c_exp_power, 'c', c_idx_of_power, '*');
            end
            
            % -- remove the last character which is '*' --
            c_exp_power_char = char(c_exp_power);
            c_exp_power = c_exp_power_char(1:end-1);
            
            % -- replace the c power by a series of multiplications --
            c_power_part = strcat('c', c_idx_of_power, '^', c_power);
            p_terms(1,j) = strrep(p_terms(1,j), c_power_part, c_exp_power);
        end
        
        % -- 3-2) within each term, sperate based on symbol c, then replace each part with parameter homotopy --
        c_parts = strsplit(p_terms(1,j), '*');
        rep_term = '';
        for k = 1:size(c_parts, 2)
            if contains(c_parts(1,k), 'c')
                c_idx = extractAfter(c_parts(1,k), 'c');
                p_homotopy = strcat('(p', c_idx, '*t+q', c_idx,'*(1-t))');
                rep_term = strcat(rep_term, p_homotopy, '*');
            else
                rep_term = strcat(rep_term, c_parts(1,k), '*');
            end
        end
        char_term = char(rep_term);
        char_term_exact = char_term(1:end-1);
        
        % -- 3-3) cacetenate all terms with original arithmatic symbols --
        if j == size(p_terms, 2)
            full_p_homotopy_str = strcat(full_p_homotopy_str, char_term_exact);
        else
            full_p_homotopy_str = strcat(full_p_homotopy_str, char_term_exact, arith_sym(1,j));
        end
    end
    
    wr_lhs_c = strcat('c', num2str(i), '=');
    fprintf(outputFileWr_phomotopy_raw, full_p_homotopy_str);
    fprintf(outputFileWr_phomotopy_raw, '\n');
    
    
    % -- 4) convert string to symbolic and factorize it based on t --
    syms t;
    sym_p_homotopy = str2sym(full_p_homotopy_str);
    exp_sym_p_homotopy = expand(sym_p_homotopy);
    collect_sym_p_homotopy = collect(exp_sym_p_homotopy, t);
    
    % -- 5) find the maximal power of t among all coefficient --
    function_t(i,1) = string(collect_sym_p_homotopy);
    [~, var] = coeffs(collect_sym_p_homotopy,t);
    if size(var, 2) > max_t_size
        max_t_size = size(var, 2);
    end
    
    % -- 6) TODO: check whether there is a square in either p or q of the coefficients of the function t --
    % -- TODO: change to a generic power --
    %[pq_terms, arith_sym_pq] = strsplit(string(coefficient(1, nc)), {'+','-'});
    %square_power = strfind(string(coefficient(1, nc)), '^2');
    
    % -- 7) write the coefficients of funciton t into a file --
    fprintf(outputFileWr_phomotopy, string(collect_sym_p_homotopy));
    fprintf(outputFileWr_phomotopy, '\n');

% ===============================================================================================
%     % -- 5) plug in complex numbers of parameters into the symbolic
%     % expressions --
%     % -- 5-1) extract all ps and qs --
%     pqs = symvar(collect_sym_p_homotopy);
%     % -- 5-2) loop over all vars --
%     for v = 1:size(pqs, 2)
%         str_var = string(pqs(v));
%         if contains(str_var, 'p')
%             var_idx = extractAfter(str_var, 'p');
%             collect_sym_p_homotopy = subs(collect_sym_p_homotopy, pqs(v), p(double(var_idx)));
%         elseif contains(str_var, 'q')
%             var_idx = extractAfter(str_var, 'q');
%             collect_sym_p_homotopy = subs(collect_sym_p_homotopy, pqs(v), q(double(var_idx)));
%         end
%     end
%     
%     % -- 6) store function t --
%     function_t(i,1) = string(collect_sym_p_homotopy);
%     
%     % -- 7) extract coefficients and variables of the function t and write into a file --
%     %[coefficient, var] = coeffs(collect_sym_p_homotopy,t);
%     sym_p_homotopy_all = [sym_p_homotopy_all; collect_sym_p_homotopy];
%     
%     % -- 8) find the maximal power of t among all coefficient --
%     [coefficient, var] = coeffs(collect_sym_p_homotopy,t);
%     if size(var, 2) > max_t_size
%         max_t_size = size(var, 2);
%     end
    
%     for nc = 1:size(coefficient, 2)
%         fprintf(outputFileWr_coeffs_of_func_t, string(double(coefficient(1,nc))));
%         fprintf(outputFileWr_coeffs_of_func_t, '\t');
%     end
%     fprintf(outputFileWr_coeffs_of_func_t, '\n');
    
    % -- show running progress --
    if mod(i,10) == 0
        fprintf('. ');
    end
end
fprintf('\n');

% -- write all ceifficient of the function t in each homotopy problem coefficient in the order of [t^0, t^1, t^2, ...] --
% coeffs_of_func_t = strings(size(sym_p_homotopy_all, 1), max_t_size);
% for i = 1:size(sym_p_homotopy_all, 1)
%     coeff_idx = 1;
%     [coefficient, var] = coeffs(str2sym(function_t(i,1)), t);
%     for j = size(coefficient, 2):-1:1
%         coeffs_of_func_t(i, coeff_idx) = double(coefficient(1, j));
%         coeff_idx = coeff_idx + 1;
%     end
% end

% -- write the coefficients of the function t to a file --
acc_numOfCoeffs = 0;
fprintf('writing the coefficients of function t for H into a file ');
for i = 1:size(line_p2c, 1)
    [coefficient, var] = coeffs(str2sym(function_t(i,1)), t);
    
    for nc = size(coefficient, 2):-1:1
        % -- find whether the degree of parameters is greater than 1 --
        [p_terms, arith_sym] = strsplit(string(coefficient(1, nc)), {'+','-'});
        full_coeff = '';
        
        % -- for each term --
        for j = 1:size(p_terms, 2)
            
            % -- split the term by multiplication operator --
            parts = strsplit(p_terms(1,j), '*');
            
            % -- for each part --
            p_terms(1,j) = '';
            for p = 1:size(parts, 2)
                if contains(parts(1,p), '^')
                    degree = extractAfter(parts(1,p), '^');
                    if contains(parts(1,p), 'p')
                        param_sym = 'p';
                        idx = extractBetween(parts(1,p), 'p', '^');
                        p_sym = extractBefore(parts(1,p), 'p');
                    else
                        param_sym = 'q';
                        idx = extractBetween(parts(1,p), 'q', '^');
                        p_sym = extractBefore(parts(1,p), 'q');
                    end
                    
                    parts(1,p) = '';
                    for d = 1:str2double(degree)
                        if d < str2double(degree)
                            parts(1,p) = strcat(parts(1,p), param_sym, idx, '*');
                        else
                            parts(1,p) = strcat(parts(1,p), param_sym, idx);
                        end
                    end
                end
                
                if p < size(parts, 2)
                    p_terms(1,j) = strcat(p_terms(1,j), parts(1,p), '*');
                else
                    p_terms(1,j) = strcat(p_terms(1,j), parts(1,p));
                end
            end
            
            % -- cacetenate all terms with original arithmatic symbols --
            if j == size(p_terms, 2)
                full_coeff = strcat(full_coeff, p_terms(1,j));
            else
                full_coeff = strcat(full_coeff, p_terms(1,j), arith_sym(1,j));
            end
        end

        % -- write to the file --
        %wrstr = strcat('h_phc_coeffs_H[', num2str(acc_numOfCoeffs), ']=', string(coefficient(1, nc)), ';');
        wrstr = strcat('h_phc_coeffs_H[', num2str(acc_numOfCoeffs), ']=', full_coeff, ';');
        fprintf(outputFileWr_coeffs_of_func_t, wrstr);
        fprintf(outputFileWr_coeffs_of_func_t, '\n');
        acc_numOfCoeffs = acc_numOfCoeffs + 1;
    end
    
    % -- pad with zeros --
    if size(coefficient, 2) < max_t_size
        for j = size(coefficient, 2):max_t_size-1
            wrstr = strcat('h_phc_coeffs_H[', num2str(acc_numOfCoeffs), ']= MAGMA_C_ZERO;');
            fprintf(outputFileWr_coeffs_of_func_t, wrstr);
            fprintf(outputFileWr_coeffs_of_func_t, '\n');
            acc_numOfCoeffs = acc_numOfCoeffs + 1;
        end
    end
    
    % -- show running progress --
    if mod(i,10) == 0
        fprintf('. ');
    end
end
fprintf('\n');


% -- differentiating the function t to get the coefficient of function t
% for Ht --
function_t_Jt = strings(size(line_p2c, 1), 1);
fprintf('writing the coefficients of function t for Ht into a file ');
for i = 1:size(function_t, 1)
    Jt = jacobian(str2sym(function_t(i,1)), t);
    exp_Jt = expand(Jt);
    collect_exp_Jt = collect(exp_Jt, t);
    [coefficient, var] = coeffs(collect_exp_Jt, t);
    
    for nc = size(coefficient, 2):-1:1
        
        % -- find whether the degree of parameters is greater than 1 --
        [p_terms, arith_sym] = strsplit(string(coefficient(1, nc)), {'+','-'});
        full_coeff = '';
        
        for j = 1:size(p_terms, 2)
            % -- split the term by multiplication operator --
            parts = strsplit(p_terms(1,j), '*');

            % -- for each part --
            p_terms(1,j) = '';
            for p = 1:size(parts, 2)
                if contains(parts(1,p), '^')
                    degree = extractAfter(parts(1,p), '^');
                    if contains(parts(1,p), 'p')
                        param_sym = 'p';
                        idx = extractBetween(parts(1,p), 'p', '^');
                        p_sym = extractBefore(parts(1,p), 'p');
                    else
                        param_sym = 'q';
                        idx = extractBetween(parts(1,p), 'q', '^');
                        p_sym = extractBefore(parts(1,p), 'q');
                    end

                    parts(1,p) = '';
                    for d = 1:str2double(degree)
                        if d < str2double(degree)
                            parts(1,p) = strcat(parts(1,p), param_sym, idx, '*');
                        else
                            parts(1,p) = strcat(parts(1,p), param_sym, idx);
                        end
                    end
                end

                if p < size(parts, 2)
                    p_terms(1,j) = strcat(p_terms(1,j), parts(1,p), '*');
                else
                    p_terms(1,j) = strcat(p_terms(1,j), parts(1,p));
                end
            end
            
            % -- cacetenate all terms with original arithmatic symbols --
            if j == size(p_terms, 2)
                full_coeff = strcat(full_coeff, p_terms(1,j));
            else
                full_coeff = strcat(full_coeff, p_terms(1,j), arith_sym(1,j));
            end
        end
        
        % -- write to the file --
        %wrstr = strcat('h_phc_coeffs_Ht[', num2str(acc_numOfCoeffs), ']=', string(coefficient(1, nc)), ';');
        wrstr = strcat('h_phc_coeffs_Ht[', num2str(acc_numOfCoeffs), ']=', full_coeff, ';');
        fprintf(outputFileWr_coeffs_of_func_t_Ht, wrstr);
        fprintf(outputFileWr_coeffs_of_func_t_Ht, '\n');
        acc_numOfCoeffs = acc_numOfCoeffs + 1;
    end
    
    % -- pad with zeros --
    if size(coefficient, 2) < max_t_size-1
        for j = size(coefficient, 2):max_t_size-2
            wrstr = strcat('h_phc_coeffs_Ht[', num2str(acc_numOfCoeffs), ']= MAGMA_C_ZERO;');
            fprintf(outputFileWr_coeffs_of_func_t_Ht, wrstr);
            fprintf(outputFileWr_coeffs_of_func_t_Ht, '\n');
            acc_numOfCoeffs = acc_numOfCoeffs + 1;
        end
    end
    
    % -- show running progress --
    if mod(i,10) == 0
        fprintf('. ');
    end
end
fprintf('\n');

fprintf("Maximal order of function t is:");
fprintf(string(max_t_size-1));
fprintf('\n');

% -- write the correspondences between start and target params arrays and
% ps as well as qs --
% -- 1) target and ps --
for i = 1:numOfParams
    wrstr = strcat('magmaFloatComplex', {' '}, 'p', num2str(i), {' '}, '=', {' '}, 'h_targetParams[', num2str(i-1), '];\n');
    fprintf(outputFileWr_param_arrays, string(wrstr));
end
% -- 2) start and qs --
for i = 1:numOfParams
    wrstr = strcat('magmaFloatComplex', {' '}, 'q', num2str(i), {' '}, '=', {' '}, 'h_startParams[', num2str(i-1), '];\n');
    fprintf(outputFileWr_param_arrays, string(wrstr));
end

fclose(outputFileWr_phomotopy);
fclose(outputFileWr_coeffs_of_func_t);
fclose(outputFileWr_coeffs_of_func_t_Ht);
fclose(outputFileWr_param_arrays);
fclose(outputFileWr_phomotopy_raw);
