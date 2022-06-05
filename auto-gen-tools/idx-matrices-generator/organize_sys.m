function organized_sys = organize_sys(sys, partial_derivative)
    
    if strcmp(partial_derivative, 'Hx')
        % -- get a collection of cdt terms --
        sys_terms = extract_cdt_term(sys);
    else
        % -- get a collection of cdt terms --
        sys_terms = extract_cd_term(sys);
    end
    
    % -- concatenate to form a newly organized system --
    sys_cat_term = '';
    for i = 1:size(sys_terms, 1)
        if i < size(sys_terms, 1)
            sys_cat_term = strcat(sys_cat_term, sys_terms(i, 1), '+');
        else
            sys_cat_term = strcat(sys_cat_term, sys_terms(i, 1));
        end
    end

    if isempty(sys_cat_term)
        organized_sys = '0';
    else
        organized_sys = sys_cat_term;
    end

    
end