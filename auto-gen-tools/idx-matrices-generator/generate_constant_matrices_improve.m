function [const_matrix_Hx, const_matrix_Ht] = generate_constant_matrices_improve(numOfVars, numOfCoeff, Hx, Ht, Hx_maximal_terms, Hx_maximal_parts, Ht_maximal_terms, Ht_maximal_parts)

    Hx_numOfX = Hx_maximal_parts - 2;
    Ht_numOfX = Ht_maximal_parts - 2;
    
    const_matrix_Hx = numOfVars*ones(numOfVars, numOfVars*Hx_maximal_terms*Hx_maximal_parts);
    const_matrix_Ht = numOfVars*ones(numOfVars, Ht_maximal_terms*Ht_maximal_parts);
    
    % -- initialize default index of coefficients in Hx --
    for i = Hx_maximal_parts:Hx_maximal_parts:size(const_matrix_Hx, 2)
        const_matrix_Hx(:,i) = numOfCoeff;
    end
    for i = 1:Hx_maximal_parts:size(const_matrix_Hx, 2)-Hx_maximal_parts+1
        const_matrix_Hx(:,i) = 0;
    end
    
    % -- initialize default index of coefficients in Ht --
    for i = Ht_maximal_parts:Ht_maximal_parts:size(const_matrix_Ht, 2)
        const_matrix_Ht(:,i) = numOfCoeff;
    end
    for i = 1:Ht_maximal_parts:size(const_matrix_Ht, 2)-Ht_maximal_parts+1
        const_matrix_Ht(:,i) = 0;
    end
    
    
    for i = 1:size(Hx, 1)
    
        % -- Hx --
        for j = 1:size(Hx, 2)

            Hx_idx = 1 + (j-1)*Hx_maximal_parts*Hx_maximal_terms;
            for k = 1:size(Hx{i,j}, 2)

                % -- 1) scalar --
                const_matrix_Hx(i, Hx_idx) = Hx{i,j}{k}(1);
                
                % -- 2) monomial --
                v = symvar(Hx{i,j}{k});
                if ~isempty(v)
                    x_index_backward = find(has(Hx{i,j}{k}, v(4:end)));
                    for xind = 1:size(x_index_backward, 2)
                        str_v = char(Hx{i,j}{k}(x_index_backward(xind)));
                        const_matrix_Hx(i, Hx_idx+xind) = str2double(erase(str_v, 'x_'));
                    end
                end
                
                % -- 3) coefficient --
                if ~isempty(v)
                    const_matrix_Hx(i, Hx_idx+Hx_maximal_parts-1) = str2double(erase(char(v(1)), 'c_'));
                end
                
                Hx_idx = Hx_idx + Hx_maximal_parts;
            end
        end

        % -- Ht --
        Ht_idx = 1;
        for m = 1:size(Ht{i}, 2)
            % -- 1) scalar --
            const_matrix_Ht(i, Ht_idx) = Ht{i}{m}(1);
            
            % -- 2) monomial --
            v = symvar(Ht{i}{m});
            x_index_backward = find(has(Ht{i}{m}, v(3:end)));
            for xind = 1:size(x_index_backward, 2)
                str_v = char(Ht{i}{m}(x_index_backward(xind)));
                const_matrix_Ht(i, Ht_idx+xind) = str2double(erase(str_v, 'x_'));
            end
            
            % -- 3) coefficient --
            const_matrix_Ht(i, Ht_idx+Ht_maximal_parts-1) = str2double(erase(char(v(1)), 'c_'));
            
            Ht_idx = Ht_idx + Ht_maximal_parts;
        end

        fprintf(". ");
    end
    fprintf("\n");
end