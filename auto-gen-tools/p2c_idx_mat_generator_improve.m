% -- VERSION 2 --

% -- generate  --
clear;

% -- define directories --
fileFolder = '/home/chchien/BrownU/research/22-issac/computer-vision-problems/';
problem = 'idx-generator-test/';
numOfParams = 81;

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

ldata = textscan(inputFileWr_p2c, '%s', 'Delimiter', ',', 'whitespace', '{}', 'CollectOutput', true);
line_p2c = ldata{1};

sym_p_homotopy_all = [];
max_t_size = 0;
acc_numOfCoeffs = 0;
function_t = strings(size(line_p2c, 1), 1);
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

    %fprintf(outputFileWr_phomotopy_raw, full_p_homotopy_str);
    %fprintf(outputFileWr_phomotopy_raw, '\n');
    
    
    % -- 4) convert string to symbolic and factorize it based on t --
    syms t;
    sym_p_homotopy = str2sym(full_p_homotopy_str);
    exp_sym_p_homotopy = expand(sym_p_homotopy);
    collect_sym_p_homotopy = collect(exp_sym_p_homotopy, t);
    
    % -- 5) find the maximal power of t among all coefficient --
    [coefficient, var] = coeffs(collect_sym_p_homotopy,t);
    if size(var, 2) > max_t_size
        max_t_size = size(var, 2);
    end
    
    % -- represent the string properly --
    function_t(i,1) = '';
    for c = 1:size(coefficient, 2)
        if c < size(coefficient, 2)
            function_t(i,1) = strcat(function_t(i,1), '(', string(coefficient(1,c)), ')*', string(var(1,c)), " + ");
        else
            function_t(i,1) = strcat(function_t(i,1), '(', string(coefficient(1,c)), ')*', string(var(1,c)));
        end
    end
%     function_t(i,1) = string(collect_sym_p_homotopy);
    
    % -- 6) write the coefficients of funciton t into a file --
    fprintf(outputFileWr_phomotopy, string(collect_sym_p_homotopy));
    fprintf(outputFileWr_phomotopy, '\n');

    % -- show running progress --
    if mod(i,10) == 0
        fprintf('. ');
    end
end
fprintf('\n');

% -- write the coefficients of the function t to a file for Hx --
acc_numOfCoeffs = 0;
fprintf('writing the coefficients of function t for Hx into a file ');
for i = 1:size(line_p2c, 1)
    % -- extract the coefficients of the polynomial t --
    if contains(function_t(i,1), 'p') || contains(function_t(i,1), 'q')
        coefficient = strsplit(function_t(i,1), " + (");
        for tc = 1:size(coefficient, 2)
            coefficient(1, tc) = extractBefore(coefficient(1, tc), ")");
            coefficient(1, tc) = strrep(coefficient(1, tc), '(', '');
        end
    else
        coefficient = function_t(i,1);
    end
   
    %[coefficient, var] = coeffs(str2sym(function_t(i,1)), t);
    
    for nc = size(coefficient, 2):-1:1
        % -- find whether the degree of parameters is greater than 1 --
        [p_terms, arith_sym] = strsplit(coefficient(1, nc), {'+','-'});
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
acc_numOfCoeffs = 0;
fprintf('writing the coefficients of function t for Ht into a file ');
for i = 1:size(function_t, 1)
    % -- 1) extract the coefficients of the polynomial t --
    if contains(function_t(i,1), 'p') || contains(function_t(i,1), 'q')
        coefficient_str_raw = strsplit(function_t(i,1), " + (");
        coefficient = strings(1, size(coefficient_str_raw, 2)-1);

        for tc = 1:size(coefficient_str_raw, 2)-1
            % -- 2) extract the degree of polynomial t --
            if contains(coefficient_str_raw(1, tc), 't^')
                str_t_degree = extractAfter(coefficient_str_raw(1, tc), 't^');
            elseif contains(coefficient_str_raw(1, tc), 't')
                str_t_degree = '1';
            end

            % -- 3) polish the string --
            coefficient(1, tc) = extractBefore(coefficient_str_raw(1, tc), ")");
            coefficient(1, tc) = strrep(coefficient(1, tc), '(', '');

            % -- 4) multiply the coefficient with the degree of t --
            if ~strcmp(str_t_degree, '1')
                coefficient(1, tc) = strcat(str_t_degree, '*(', coefficient(1, tc), ')');
            end
        end
    else
        coefficient = '';
    end
    
%     Jt = jacobian(str2sym(function_t(i,1)), t);
%     exp_Jt = expand(Jt);
%     collect_exp_Jt = collect(exp_Jt, t);
%     [coefficient, var] = coeffs(collect_exp_Jt, t);
    
    for nc = size(coefficient, 2):-1:1
        
        % -- find whether the degree of parameters is greater than 1 --
        [p_terms, arith_sym] = strsplit(coefficient(1, nc), {'+','-'});
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
