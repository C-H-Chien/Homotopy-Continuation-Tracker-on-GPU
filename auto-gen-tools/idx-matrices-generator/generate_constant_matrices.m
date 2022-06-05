function [Ht_numOfX, const_matrix_Hx_scalar, const_matrix_Hx_x, const_matrix_Hx_Y, const_matrix_Ht_scalar, const_matrix_Ht_x, const_matrix_Ht_Y] = ...
    generate_constant_matrices(numOfVars, numOfCoeff, Hx, Ht, Hx_maximal_terms, Hx_maximal_parts, Ht_maximal_terms, Ht_maximal_parts)

    Hx_numOfX = Hx_maximal_parts - 2;
    Ht_numOfX = Ht_maximal_parts - 2;
    
    const_matrix_Hx_scalar = zeros(numOfVars*numOfVars, Hx_maximal_terms);
    const_matrix_Hx_x = numOfVars*ones(numOfVars*numOfVars, Hx_maximal_terms, Hx_numOfX);
    const_matrix_Hx_Y = numOfCoeff*ones(numOfVars*numOfVars, Hx_maximal_terms);
    const_matrix_Ht_scalar = zeros(numOfVars, Ht_maximal_terms);
    const_matrix_Ht_x = numOfVars*ones(numOfVars, Ht_maximal_terms, Ht_numOfX);
    const_matrix_Ht_Y = zeros(numOfVars, Ht_maximal_terms);

    for i = 1:size(Hx, 1)
    
        % -- Hx --
        for j = 1:size(Hx, 2)
            for k = 1:size(Hx{j,i}, 2)

                % -- 1) scalar matrix --
                const_matrix_Hx_scalar((i-1)*numOfVars+j, k) = Hx{j,i}{k}(1);

                v = symvar(Hx{j,i}{k});

                % -- 2) X matrices --
                % -- use forward and backward find to make sure the power of x
                % is involved --
                %x_index_forward = find(has(v, Hx{j,i}{k}));
                if ~isempty(v)
                    x_index_backward = find(has(Hx{j,i}{k}, v(4:end)));
                    for xind = 1:size(x_index_backward, 2)
                        str_v = char(Hx{j,i}{k}(x_index_backward(xind)));
                        const_matrix_Hx_x((i-1)*numOfVars+j, k, xind) = str2double(erase(str_v, 'x_'));
                    end
                end

                % -- 3) Y matrix --
                if ~isempty(v)
                    const_matrix_Hx_Y((i-1)*numOfVars+j, k) = str2double(erase(char(v(1)), 'c_'));
                end

            end
        end

        % -- Ht --
        for m = 1:size(Ht{i}, 2)
            % -- 1) scalar matrix --
            const_matrix_Ht_scalar(i, m) = Ht{i}{m}(1);

            v = symvar(Ht{i}{m});
            % -- 2) X matrices --
            % -- use forward and backward find to make sure the power of x
            % is involved --
            % -- c and d are two different terms... --
            x_index_backward = find(has(Ht{i}{m}, v(3:end)));
            for xind = 1:size(x_index_backward, 2)
                str_v = char(Ht{i}{m}(x_index_backward(xind)));
                const_matrix_Ht_x(i, m, xind) = str2double(erase(str_v, 'x_'));
            end

            % -- 3) Y matrix --
            
            const_matrix_Ht_Y(i, m) = str2double(erase(char(v(1)), 'c_'));
        end

%         % -- H --
%         for h_parts_indx = 1:size(H{i}, 2)
%             % -- 1) scalar matrix --
%             const_matrix_H_scalar(i, h_parts_indx) = H{i}{h_parts_indx}(1);
% 
%             v = symvar(H{i}{h_parts_indx});
%             % -- 2) X matrices --
%             x_index_backward = find(has(H{i}{h_parts_indx}, v(4:end)));
%             for xind = 1:size(x_index_backward, 2)
%                 str_v = char(H{i}{h_parts_indx}(x_index_backward(xind)));
%                 const_matrix_H_x(i, h_parts_indx, xind) = str2double(erase(str_v, 'x_'));
%             end
% 
%             % -- 3) Y matrix --
%             const_matrix_H_Y(i, h_parts_indx) = str2double(erase(char(v(1)), 'c_'));
%         end

        fprintf(". ");
    end
    fprintf("\n");
end