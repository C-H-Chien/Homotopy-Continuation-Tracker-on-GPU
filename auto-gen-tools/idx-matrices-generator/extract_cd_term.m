function cd_terms = extract_cd_term(sys)
    %f_sys = factor(sys);
    f_sys = sys;
    split_f_sys = split(char(f_sys), {'+', '-'});
    
    % -- loop over all split elements and construct an array of cdt terms --
    cd_terms = [];
    for i = 1:size(split_f_sys, 1)
        split_f_sys{i} = strrep(split_f_sys{i}, ' ', '');
        if contains(split_f_sys{i}, 'c_')
            
            % -- if the term contains monomial --
            if contains(split_f_sys{i}, '*') 
                if strcmp(split_f_sys{i}(1), 'c')
                    monomial = extractAfter(string(split_f_sys{i}), '*');

                    c_indx = extractBetween(string(split_f_sys{i}), 'c_', '*');

                    str_of_cd_term = strcat('c_', c_indx, ' - ', 'd_', c_indx);
                    cd_with_monomial = strcat(monomial, '*(', str_of_cd_term, ')');
                    if ~isempty(cd_with_monomial)
                        cd_with_monomial_rep = strrep(cd_with_monomial, ' ', '');
                    end
                else
                    scalar = extractBefore(string(split_f_sys{i}), '*');
                    sstr = strcat(scalar, '*');
                    cd_monomial = extractAfter(string(split_f_sys{i}), sstr);
                    
                    monomial = extractAfter(cd_monomial, '*');

                    c_indx = extractBetween(string(split_f_sys{i}), 'c_', '*');

                    str_of_cd_term = strcat('c_', c_indx, ' - ', 'd_', c_indx);
                    full_str = strcat(sstr, monomial, '*(', str_of_cd_term, ')');
                    if ~isempty(full_str)
                        cd_with_monomial_rep = strrep(full_str, ' ', '');
                    end
                end
            else % -- otherwise, when the term only has cd --
                c_indx = extractAfter(string(split_f_sys{i}), 'c_');
                cd_with_monomial = strcat('c_', c_indx, ' - ', 'd_', c_indx);
                cd_with_monomial_rep = strrep(cd_with_monomial, ' ', '');
            end
            cd_terms = [cd_terms; cd_with_monomial_rep];
        else
            continue;
        end
    end
end