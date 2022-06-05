% -- brute-force parametric homotopy approach: homogenize the Hx and Ht coefficients in terms of parameter homotopy --
clear;

%sympref('FloatingPointOutput',true);

% -- define directories --
fileFolder = '/home/chchien/BrownU/research/22-issac/parametric-HC/';
problem = 'pHC-3vRelPose/';
numOfParams = 36;
numOfCoeffs = 153;

% -- input file --
% -- 1) parameters to coefficient conversion file --
inputName_p2c = 'params2coeffs.txt';
fullInputFileName_p2c = fullfile(fileFolder, problem, inputName_p2c);
inputFileWr_p2c = fopen(fullInputFileName_p2c, 'r');

% -- output files --
outputFileName_phomotopy_exp_c = 'params_homotopy_coeffs_expanded_c.txt';
outputFileName_coeffs_idx_mat = 'coeffs_idx_mat.txt';
outpufFileName_diff_phomotopy_c = 'params_homotopy_coeffs_derivative_c.txt';
fullOutputFileName_phomotopy_exp_c = fullfile(fileFolder, problem, outputFileName_phomotopy_exp_c);
fullOutputFileName_coeffs_idx_mat = fullfile(fileFolder, problem, outputFileName_coeffs_idx_mat);
fullOutputFileName_phomotopy_diff_c = fullfile(fileFolder, problem, outpufFileName_diff_phomotopy_c);
outputFileWr_phomotopy_exp_c = fopen(fullOutputFileName_phomotopy_exp_c, 'w');
outputFileWr_coeffs_idx_mat = fopen(fullOutputFileName_coeffs_idx_mat, 'w');
outputFileWr_phomotopy_diff_c = fopen(fullOutputFileName_phomotopy_diff_c, 'w');

ldata = textscan(inputFileWr_p2c, '%s', 'Delimiter', ',', 'whitespace', '{}', 'CollectOutput', true);
line_p2c = ldata{1};

max_terms = 0;
max_parts = 0;

phomotopyOfCoeffs = strings(size(line_p2c, 1), 1);

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
    expanded_phomotopy_str = '';
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

        % -- 3-2) concatenate all terms with original arithmatic symbols --
        if j == size(p_terms, 2)
            expanded_phomotopy_str = strcat(expanded_phomotopy_str, p_terms(1,j));
        else
            expanded_phomotopy_str = strcat(expanded_phomotopy_str, p_terms(1,j), arith_sym(1,j));
        end
        
        % -- 3-3) find maximal parts --
        [elements, multiplication] = strsplit(p_terms(1,j), '*');
        if size(elements, 2) > max_parts
            max_parts = size(elements, 2);
        end
    end
    
    % -- 4) find maximal terms --
    if size(p_terms, 2) > max_terms
        max_terms = size(p_terms, 2);
    end
    
    % -- 5) store all expanded paramtric homotopy --
    phomotopyOfCoeffs(i,1) = expanded_phomotopy_str;
    
    % -- 7) write the coefficients of funciton t into a file --
    fprintf(outputFileWr_phomotopy_exp_c, expanded_phomotopy_str);
    fprintf(outputFileWr_phomotopy_exp_c, '\n');
    
    % -- show running progress --
    if mod(i,10) == 0
        fprintf('. ');
    end
end
fprintf('\n');

fprintf("max terms = ");
fprintf(string(max_terms));
fprintf('\n');
fprintf("max parts = ");
fprintf(string(max_parts));
fprintf('\n');

% -- declare index matrix for parametric homotopy --
bf_index_mat = zeros(numOfCoeffs, max_terms*max_parts);
for i = 1:size(phomotopyOfCoeffs, 1)
    % -- 1) splite terms based on arithmetic operations --
    [p_terms, arith_sym] = strsplit(phomotopyOfCoeffs(i,1), {'+','-'});
    
    % -- 2) for each term... --
    for j = 1:size(p_terms, 2)
        % -- 3) splite each term into parts --
        [elements, multiplication] = strsplit(p_terms(1,j), '*');
        %parts_idx = 1;
        store_idx = (j-1)*max_parts + 1;
        
        % -- 4) for each part scalar ... --
        if contains(elements(1,1), 'c')
            % -- no explicit scalar --
            cat_arith = "";
            if j > 1
                if strcmp(arith_sym(1,j-1), '-')
                    cat_arith = arith_sym(1,j-1);
                end
            end
            bf_index_mat(i,store_idx) = double(strcat(cat_arith, "1"));
            store_idx = store_idx + 1;
            
            % -- for the c parts --
            for k = 1:size(elements, 2)
                c_idx = extractAfter(elements(1,k), 'c');
                bf_index_mat(i,store_idx) = double(c_idx)-1;
                store_idx = store_idx + 1;
            end
        elseif strcmp(elements(1,1), "")
            % -- only a scalar --
            bf_index_mat(i,store_idx) = double(strcat(arith_sym, p_terms(1,2)));
            store_idx = store_idx + 1;
            for s = store_idx:max_parts
                bf_index_mat(i,store_idx) = numOfParams;
                store_idx = store_idx + 1;
            end
            break;
        else
            % -- explicit scalar --
            cat_arith = "";
            if j > 1
                if strcmp(arith_sym(1,j-1), '-')
                    cat_arith = arith_sym(1,j-1);
                end
            end
            bf_index_mat(i,store_idx) = double(strcat(cat_arith, elements(1,1)));
            store_idx = store_idx + 1;
            
            % -- for the c parts --
            for k = 2:size(elements, 2)
                c_idx = extractAfter(elements(1,k), 'c');
                bf_index_mat(i,store_idx) = double(c_idx)-1;
                store_idx = store_idx + 1;
            end
        end      
        
        % -- 6) pad the remaining parts of the bf_index_mat --
        if store_idx <= j*max_parts
            for s = store_idx:j*max_parts
                bf_index_mat(i,store_idx) = numOfParams;
                store_idx = store_idx + 1;
            end
        end
    end
end

% -- construct a derivative of parametric homotopy w.r.t. t --
for i = 1:size(phomotopyOfCoeffs, 1)
    % -- 0) define array name in the gpu code --
    lin_interp_array_name = 'c';
    derivative_array_name = 'd';
    %lin_interp_array_name = 's_lin_interp_params[';
    %derivative_array_name = 'd_const_vec_cd[';
    
    % -- 1) splite terms based on arithmetic operations --
    [p_terms, arith_sym] = strsplit(phomotopyOfCoeffs(i,1), {'+','-'});
    
    % -- 2) for each term... --
    derivative_phomotopy_str = '';
    for j = 1:size(p_terms, 2)
        % -- 3) splite each term into parts --
        [elements, multiplication] = strsplit(p_terms(1,j), '*');
        
        % -- 4) check how many multiplied cs are there in each term --
        numOfC = 0;
        c_indices = strings(2,1);
        for k = 1:size(elements, 2)
            if contains(elements(1,k), 'c')
                numOfC = numOfC + 1;
                c_indices(numOfC, 1) = extractAfter(elements(1,k), 'c');
            end
        end
        
        % -- 5) replace the term with the derivative version --
        if numOfC == 2
            target_str = strcat('c', c_indices(1,1), '*c', c_indices(2,1));
            c_indices(1,1) = string(double(c_indices(1,1)) - 1);
            c_indices(2,1) = string(double(c_indices(2,1)) - 1);
            %derivative_representation = strcat('(', lin_interp_array_name, c_indices(1,1), ']*', derivative_array_name, c_indices(2,1), ']+', lin_interp_array_name, c_indices(2,1), ']*', derivative_array_name, c_indices(1,1), '])');
            derivative_representation = strcat('(', lin_interp_array_name, c_indices(1,1), '*', derivative_array_name, c_indices(2,1), '+', lin_interp_array_name, c_indices(2,1), '*', derivative_array_name, c_indices(1,1), ')');
            derivative_p_terms = strrep(p_terms(1,j), target_str, derivative_representation);
        elseif numOfC == 0
            % -- only scalars, set it as zero directly --
            derivative_p_terms = "0";
        else
            target_str = strcat('c', c_indices(1,1));
            c_indices(1,1) = string(double(c_indices(1,1)) - 1);
            %derivative_representation = strcat(derivative_array_name, c_indices(1,1), ']');
            derivative_representation = strcat(derivative_array_name, c_indices(1,1));
            derivative_p_terms = strrep(p_terms(1,j), target_str, derivative_representation);
        end

        % -- 6) concatenate all terms with original arithmatic symbols --
        if j == size(p_terms, 2)
            derivative_phomotopy_str = strcat(derivative_phomotopy_str, derivative_p_terms);
        else
            derivative_phomotopy_str = strcat(derivative_phomotopy_str, derivative_p_terms, arith_sym(1,j));
        end
    end
    
    % -- 7) write to a file --
    lhs = strcat('s_phc_coeffs[', num2str(i-1), ']=');
    write_str = strcat(lhs, derivative_phomotopy_str, ';');
    fprintf(outputFileWr_phomotopy_diff_c, write_str);
    fprintf(outputFileWr_phomotopy_diff_c, '\n');
end

% -- write bf_indx_mat into a file --
for i = 1:size(bf_index_mat, 1)
    for j = 1:size(bf_index_mat, 2)
        fprintf(outputFileWr_coeffs_idx_mat, num2str(bf_index_mat(i,j)));
        fprintf(outputFileWr_coeffs_idx_mat, '\t');
    end
    fprintf(outputFileWr_coeffs_idx_mat, '\n');
end

fclose(outputFileWr_phomotopy_exp_c);
fclose(outputFileWr_coeffs_idx_mat);
fclose(outputFileWr_phomotopy_diff_c);
