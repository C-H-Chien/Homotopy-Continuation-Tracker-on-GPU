function ret_terms = inter_terms_division(ch, cat_terms, flag)
    recursive_flag = 0;
    terms = [];
    for k = size(ch, 2):-1:1
        v = symvar(ch{k});
        if isempty(v)
            terms = ch{k};
            break;
        end
        
        if size(v,2) == 1
            solution_part = ch{k};
            terms = [cat_terms, solution_part];
            cat_terms = terms;
        elseif size(v,2) == 3
            f = factor(ch{k});
            if size(f,2) > 1
                scalar_part = f(1,1);
                cdt_part = f(1,2);
                terms = [cat_terms, scalar_part, cdt_part];
            else
                cdt_part = f;
                terms = [cat_terms, cdt_part];
            end
            cat_terms = terms;
        else
            ch2 = children(ch{k});
            recursive_ret = inter_terms_division(ch2, [], 1);
            recursive_flag = 1;
        end
    end
    
    if flag == 1
        ret_terms = terms;
    elseif flag == 0 && recursive_flag == 1
        ret_terms = [terms, recursive_ret];
    else
        ret_terms = terms;
    end
end