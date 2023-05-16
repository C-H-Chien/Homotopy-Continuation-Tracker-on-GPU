function [has_power, expanded_input_term] = expand_powers_in_poly_term(input_term)

%> Code Description: Given a term from the polynomial system, check whether 
%  the variables in that term has a power. If it does, expand the variable
%  so that no power symbol exists in the term. If not, simply return the
%  input term.
%
%> Inputs: 
%     input_term:          The term of the polynomial system in string format
%
%> Outputs:
%     has_power:           whether the input term has a power symbol
%     expanded_input_term: expanded term represetation without a power symbol      
%
%> (c) LEMS, Brown University
%> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
%> Completed in Oct. 1st, 2022

    %> check if there is a power
    if contains(input_term, "^")
        has_power = 1;
        
        %> first split the input_term ny multiplication operator
        parts = strsplit(input_term, "*");
        
        %> loop over all parts splitted by the input input_term
        for k = 1:size(parts, 2)
            
            %> find the part that has the power symbol
            if contains(parts(1,k), "^")
                %> find variable indices with powers
                power_vars = extractBetween(parts(1,k), "x", "^");
                
                %> find the power
                power = extractAfter(parts(1,k), "^");
                be_replaced = strcat("x", power_vars, "^", power);
                rep_part = "";
                
                %> concatenate the variable based on the power
                for p = 1:double(power)
                    if p > 1
                        rep_part = strcat(rep_part, "*x", power_vars);
                    else
                        rep_part = strcat(rep_part, "x", power_vars);
                    end
                end

                %> replace the original part with the power expanded representation
                parts(1,k) = strrep(parts(1,k), be_replaced, rep_part);
                
            else
                continue;
            end
        end
        
        %> concatenate all the parts in an expanded input_term
        expanded_input_term = "";
        for k = 1:size(parts, 2)
            if k < size(parts, 2)
                expanded_input_term = strcat(expanded_input_term, parts(1,k), "*");
            else
                expanded_input_term = strcat(expanded_input_term, parts(1,k));
            end
        end
    else
        %> if there exists no power in the input input_term
        has_power = 0;
        expanded_input_term = input_term;
    end

end