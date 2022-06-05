function cdt_terms = extract_cdt_term(sys)
    f_sys = expand(sys);
    %f_sys = sys;
    split_f_sys = split(char(f_sys), {'+', '-'});
    
%     if contains(char(f_sys), '[')
%         %split_f_sys = split(char(sys), 'z');
%         split_f_sys = split(char(sys), '-');
%     else
%         split_f_sys = split(char(f_sys), {'+', '-'});
%     end
    
    % -- loop over all split elements and construct an array of cdt terms --
    cdt_terms = [];
    for i = 1:size(split_f_sys, 1)
        if size(split_f_sys, 1) == 1
            cdt_terms = split_f_sys{i};
            break;
        else
            if contains(split_f_sys{i}, 'c_')
                % -- if the term contains monomial --
                if contains(string(split_f_sys{i}), 't*')
                    monomial = extractAfter(string(split_f_sys{i}), 't*');

                    c_indx = extractBetween(string(split_f_sys{i}), 'c_', '*t');

                    % -- extract the scalar --
                    if strfind(string(split_f_sys{i}), 'c_') ~= 2
                        cdt_scalar = extractBefore(string(split_f_sys{i}), 'c_');
                    else
                        cdt_scalar = '';
                    end

                    str_of_cdt_term = strcat(cdt_scalar, 'c_', c_indx, '*t - ', cdt_scalar, 'd_', c_indx, '*(t - 1)');
                    cdt_with_monomial = strcat(monomial, '*(', str_of_cdt_term, ')');
                    cdt_with_monomial = strrep(cdt_with_monomial, ' ', '');
                    cdt_terms = [cdt_terms; cdt_with_monomial];
                else % -- otherwise, when the term only has a cdt without any monomial --
                    c_indx = extractBetween(string(split_f_sys{i}), 'c_', '*t');

                    % -- extract the scalar --
                    if strfind(string(split_f_sys{i}), 'c_') ~= 2
                        cdt_scalar = extractBefore(string(split_f_sys{i}), 'c_');
                    else
                        cdt_scalar = '';
                    end

                    str_of_cdt_term = strcat(cdt_scalar, 'c_', c_indx, '*t - ', cdt_scalar, 'd_', c_indx, '*(t - 1)');
                    str_of_cdt_term = strrep(str_of_cdt_term, ' ', '');
                    cdt_terms = [cdt_terms; str_of_cdt_term];
                end
            else
                continue;
            end
        end
    end
end