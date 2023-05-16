function [max_order_of_func_t] = createUnivPoly4Coeffs(line_p2c, outputFileWr_coeffs_Hx, outputFileWr_coeffs_Ht)

%     % -- write the correspondences between start and target params arrays and
%     % ps as well as qs --
%     % -- 1) target and ps --
%     for i = 1:numOfParams
%         wrstr = strcat('magmaFloatComplex', {' '}, 'p', num2str(i), {' '}, '=', {' '}, 'h_targetParams[', num2str(i-1), '];\n');
%         fprintf(outputFileWr_coeffs_Hx, string(wrstr));
%     end
%     % -- 2) start and qs --
%     for i = 1:numOfParams
%         wrstr = strcat('magmaFloatComplex', {' '}, 'q', num2str(i), {' '}, '=', {' '}, 'h_startParams[', num2str(i-1), '];\n');
%         fprintf(outputFileWr_coeffs_Hx, string(wrstr));
%     end
%     fprintf(outputFileWr_coeffs_Hx, "\n");

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

        %> 2) splite the RHS into terms and arithmetic symbols --
        [p_terms, arith_sym] = strsplit(str_RHS, {'+','-'});

        %> 3) loop over all terms 
        full_p_homotopy_str = '';
        for j = 1:size(p_terms, 2)
            
            %> for each term, split into parts
            [p_parts, ~] = strsplit(p_terms(1,j), "*");
            
            %> loop over all parts do:
            %> (i) check whether it is a scalar,
            %> (ii) check whether there is a power,
            %> (iii) construct parameter homotopy, and
            %> (iv) concatenate all together
            target_term = "";
            for k = 1:size(p_parts, 2)
                %> if it is a parmeter
                if contains(p_parts(1,k), 'p')
                    %> check the power
                    if contains(p_parts(1,k), '^')
                        c_power = extractAfter(p_parts(1,k), '^');
                        c_idx_of_power = extractBetween(p_parts(1,k), 'p', '^');
                        c_exp_power = '';
                        %> generate the power expansion
                        for np = 1:double(c_power)
                            c_exp_power = strcat(c_exp_power, 'p', c_idx_of_power, '*');
                        end
                        
                        %> remove the last character which is '*'
                        c_exp_power_char = char(c_exp_power);
                        c_exp_power = c_exp_power_char(1:end-1);
                        cat_part = string(c_exp_power);
                    else
                        cat_part = p_parts(1,k);
                    end
                else
                    %> it is a scalar
                    cat_part = p_parts(1,k);
                end
                
                target_term = strcat(target_term, cat_part, '*');
            end
            
            %> remove the last character which is '*'
            target_term_char = char(target_term);
            target_term = target_term_char(1:end-1);
            target_term = string(target_term);

            % -- 3-2) within each term, replace each part with parameter homotopy --
            %c_parts = strsplit(p_terms(1,j), '*');
            c_parts = strsplit(target_term, '*');
            rep_term = '';
            for k = 1:size(c_parts, 2)
                if contains(c_parts(1,k), 'p')
                    c_idx = extractAfter(c_parts(1,k), 'p');
                    p_homotopy = strcat('(p', c_idx, '*t+q', c_idx,'*(1-t))');
                    rep_term = strcat(rep_term, p_homotopy, '*');
                else
                    rep_term = strcat(rep_term, c_parts(1,k), '*');
                end
            end
            
            %> remove the last character which is '*'
            char_term = char(rep_term);
            char_term_exact = char_term(1:end-1);

            % -- 3-3) cacetenate all terms with original arithmatic symbols --
            if j == size(p_terms, 2)
                full_p_homotopy_str = strcat(full_p_homotopy_str, char_term_exact);
            else
                full_p_homotopy_str = strcat(full_p_homotopy_str, char_term_exact, arith_sym(1,j));
            end
        end

        wr_lhs_c = strcat('p', num2str(i), '=');
        %fprintf(outputFileWr_phomotopy_raw, full_p_homotopy_str);
        %fprintf(outputFileWr_phomotopy_raw, '\n');


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
        %fprintf(outputFileWr_phomotopy, string(collect_sym_p_homotopy));
        %fprintf(outputFileWr_phomotopy, '\n');

        % -- show running progress --
        if mod(i,10) == 0
            fprintf('. ');
        end
    end
    fprintf('\n');

    % -- write the coefficients of the function t to a file --
    acc_numOfCoeffs = 0;
    fprintf('Creating and writing the coefficients of function t for the Jacobian Hx ');
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

            %> Issues found by Yaqing Ding: converting full_coeff from
            %  string to sym and do simpify() is unnecessary. Otherwise if
            %  the parameter has a power, simplify() will convert the
            %  expansion form back to the power form.
            %sym_full_coeff = str2sym(full_coeff);
            %simp_sym_full_coeff = simplify(sym_full_coeff);
            %simp_str_full_coeff = string(simp_sym_full_coeff);
            %wrstr = strcat('h_phc_coeffs_Hx[', num2str(acc_numOfCoeffs), ']=', simp_str_full_coeff, ';');
            wrstr = strcat('h_phc_coeffs_Hx[', num2str(acc_numOfCoeffs), ']=', full_coeff, ';');

            % -- write to the file --
            %wrstr = strcat('h_phc_coeffs_H[', num2str(acc_numOfCoeffs), ']=', full_coeff, ';');
            fprintf(outputFileWr_coeffs_Hx, wrstr);
            fprintf(outputFileWr_coeffs_Hx, '\n');
            acc_numOfCoeffs = acc_numOfCoeffs + 1;
        end

        % -- pad with zeros --
        if size(coefficient, 2) < max_t_size
            for j = size(coefficient, 2):max_t_size-1
                wrstr = strcat('h_phc_coeffs_Hx[', num2str(acc_numOfCoeffs), "] = MAGMA_C_ZERO;");
                fprintf(outputFileWr_coeffs_Hx, wrstr);
                fprintf(outputFileWr_coeffs_Hx, '\n');
                acc_numOfCoeffs = acc_numOfCoeffs + 1;
            end
        end

        % -- show running progress --
        if mod(i,10) == 0
            fprintf('. ');
        end
    end
    fprintf('\n');
    fprintf(outputFileWr_coeffs_Hx, '\n');

    %> Differentiating the function t to get the coefficients of function t
    % for Ht --
    acc_numOfCoeffs = 0;
    function_t_Jt = strings(size(line_p2c, 1), 1);
    fprintf('Creating and writing the coefficients of function t for the Jacobian Ht ');
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

            %> Issues found by Yaqing Ding: converting full_coeff from
            %  string to sym and do simpify() is unnecessary. Otherwise if
            %  the parameter has a power, simplify() will convert the
            %  expansion form back to the power form.
            %sym_full_coeff = str2sym(full_coeff);
            %simp_sym_full_coeff = simplify(sym_full_coeff);
            %simp_str_full_coeff = string(simp_sym_full_coeff);
            %wrstr = strcat('h_phc_coeffs_Hx[', num2str(acc_numOfCoeffs), ']=', simp_str_full_coeff, ';');
            wrstr = strcat('h_phc_coeffs_Ht[', num2str(acc_numOfCoeffs), ']=', full_coeff, ';');
            
            
%             sym_full_coeff = str2sym(full_coeff);
%             simp_sym_full_coeff = simplify(sym_full_coeff);
%             simp_str_full_coeff = string(simp_sym_full_coeff);
%             wrstr = strcat('h_phc_coeffs_Ht[', num2str(acc_numOfCoeffs), ']=', simp_str_full_coeff, ';');
            %wrstr = strcat('h_phc_coeffs_Ht[', num2str(acc_numOfCoeffs), ']=', full_coeff, ';');
            fprintf(outputFileWr_coeffs_Ht, wrstr);
            fprintf(outputFileWr_coeffs_Ht, '\n');
            acc_numOfCoeffs = acc_numOfCoeffs + 1;
        end

        % -- pad with zeros --
        if size(coefficient, 2) < max_t_size-1
            for j = size(coefficient, 2):max_t_size-2
                wrstr = strcat('h_phc_coeffs_Ht[', num2str(acc_numOfCoeffs), ']= MAGMA_C_ZERO;');
                fprintf(outputFileWr_coeffs_Ht, wrstr);
                fprintf(outputFileWr_coeffs_Ht, '\n');
                acc_numOfCoeffs = acc_numOfCoeffs + 1;
            end
        end

        % -- show running progress --
        if mod(i,10) == 0
            fprintf('. ');
        end
    end
    fprintf('\n');

    %> return output
    max_order_of_func_t = max_t_size-1;
    
%     fprintf("Maximal order of function t is:");
%     fprintf(string(max_t_size-1));
%     fprintf('\n');

   
end

