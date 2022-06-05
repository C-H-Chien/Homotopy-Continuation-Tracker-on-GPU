function [Hx, Ht, H] = split_parts_in_each_term(partial_H_partial_x, partial_H_partial_t, Ht_scalar_part)

    for i = 1:size(partial_H_partial_x, 1)
        
        % -- Hx --
        for j = 1:size(partial_H_partial_x, 2)

            split_res = strsplit(partial_H_partial_x(i,j),'+');

            % -- factorize again and feed to recursive inter terms division --
            Hx_tm_indx = 1;
            for tm = 1:size(split_res, 2)
                sys_term = str2sym(split_res{tm});

                % -- check whether sys_term is cdt-only --
                if contains(split_res{tm}, 'x_')
                    ch_g = children(sys_term);
                    Hx{i,j}{Hx_tm_indx} = inter_terms_division(ch_g, [], 0);
                else
                    if strcmp(split_res{tm}(1), 'c')
                        Hx{i,j}{Hx_tm_indx} = sys_term;
                    else
                        %scalar_star = strcat(split_res{tm}(1), '*');
                        %remove_scalar_sym = str2sym(erase(split_res{tm}, scalar_star));
                        %exp_remove_scalar_sym = expand(remove_scalar_sym);
                        %Hx{i,j}{Hx_tm_indx} = [str2sym(split_res{tm}(1)), exp_remove_scalar_sym];
                        Hx{i,j}{Hx_tm_indx} = factor(sys_term);
                        Hx_tm_indx = Hx_tm_indx + 1;
                        continue;
                    end
                end

                % -- make sure that the "power" of x is expanded to
                % multiple terms --
                for p = 1:size(Hx{i,j}{Hx_tm_indx}, 2)
                    sf = strfind(char(Hx{i,j}{Hx_tm_indx}(p)), '^');
                    if isempty(sf)
                        continue;
                    else
                        power = str2double(extractAfter(char(Hx{i,j}{Hx_tm_indx}(p)), '^'));
                        Hx{i,j}{Hx_tm_indx}(p) = extractBefore(char(Hx{i,j}{Hx_tm_indx}(p)), '^');
                        for ex_parts = 1:power-1
                            Hx{i,j}{Hx_tm_indx} = [Hx{i,j}{Hx_tm_indx}, Hx{i,j}{Hx_tm_indx}(p)];
                        end
                    end
                end

                % -- put a scalar in the foremost index if the term does not have one --
                res = strfind(char(Hx{i,j}{Hx_tm_indx}(1)), 'c_');
                if ~isempty(res)
                    Hx{i,j}{Hx_tm_indx} = [str2sym('1'), Hx{i,j}{Hx_tm_indx}];
                end

                Hx_tm_indx = Hx_tm_indx + 1;
            end
        end
        
        fprintf('hx. ');

        % -- Ht --
        % -- not sure why it is needed to convert back and forth between sym and string to get correct results --
        %tmp = str2sym(char(partial_H_partial_t(i,1)));
        split_res = strsplit(partial_H_partial_t(i,1),'+');

        % -- factorize again and feed to recursive inter terms division --
        Ht_tm_indx = 1;
        special_scalar = 0;
        for ht_tm = 1:size(split_res, 2)
            sys_term = str2sym(split_res{ht_tm});

            % -- check whether sys_term is cd-only --
            if contains(split_res{ht_tm}, 'x_')
                split_str = split(split_res{ht_tm}, '*');
                if strcmp(split_str(1), Ht_scalar_part)                    
                    scalar_star = strcat(Ht_scalar_part, '*');
                    remove_scalar_str = extractAfter(split_res{ht_tm}, scalar_star);
                    ch_g = children(str2sym(remove_scalar_str));
                    Ht{i}{Ht_tm_indx} = inter_terms_division(ch_g, [], 0);
                    special_scalar = 1;
                else
                    ch_g = children(sys_term);
                    Ht{i}{Ht_tm_indx} = inter_terms_division(ch_g, [], 0);
                end
            else
                %Ht{i}{Ht_tm_indx} = sys_term;
                Ht{i}{Ht_tm_indx} = ['1', sys_term];
                Ht_tm_indx = Ht_tm_indx + 1;
                continue;
            end

            % -- make sure that the "power" of x is expanded to
            % multiple terms --
            for ht_p = 1:size(Ht{i}{Ht_tm_indx}, 2)
                sf = strfind(char(Ht{i}{Ht_tm_indx}(ht_p)), '^');
                if isempty(sf)
                    continue;
                else
                    power = str2double(extractAfter(char(Ht{i}{Ht_tm_indx}(ht_p)), '^'));
                    Ht{i}{Ht_tm_indx}(ht_p) = extractBefore(char(Ht{i}{Ht_tm_indx}(ht_p)), '^');
                    for ex_parts = 1:power-1
                        Ht{i}{Ht_tm_indx} = [Ht{i}{Ht_tm_indx}, Ht{i}{Ht_tm_indx}(ht_p)];
                    end
                end
            end

            % -- the scalar part of Ht --
            if special_scalar
                Ht{i}{Ht_tm_indx} = [str2sym(Ht_scalar_part), Ht{i}{Ht_tm_indx}];
                special_scalar = 0;
            else
                Ht{i}{Ht_tm_indx} = ['1', Ht{i}{Ht_tm_indx}];
            end
            

            Ht_tm_indx = Ht_tm_indx + 1;
        end
        
        fprintf('ht. ')

%         % -- H --
%         % -- not sure why it is needed to convert back and forth between sym and string to get correct results --
%         %tmp = str2sym(char(org_H(i,1)));
%         %split_res = strsplit(char(tmp),'+');
%         %split_res = split_terms_of_cdt(org_H(i,1), C, D);
%         split_res = strsplit(org_H(i,1),'+');
% 
%         % -- factorize again and feed to recursive inter terms division --
%         H_tm_indx = 1;
%         for h_tm = 1:size(split_res, 2)
%             sys_term = str2sym(split_res{h_tm});
% 
%             % -- check whether sys_term is cd-only --
%             if contains(split_res{h_tm}, 'x_')
%                 ch_g = children(sys_term);
%                 H{i}{H_tm_indx} = inter_terms_division(ch_g, [], 0);
%             else
%                 H{i}{H_tm_indx} = sys_term;
%             end
% 
%             % -- make sure that the "power" of x is expanded to
%             % multiple terms --
%             for ht_p = 1:size(H{i}{H_tm_indx}, 2)
%                 sf = strfind(char(H{i}{H_tm_indx}(ht_p)), '^');
%                 if isempty(sf)
%                     continue;
%                 else
%                     power = str2double(extractAfter(char(H{i}{H_tm_indx}(ht_p)), '^'));
%                     H{i}{H_tm_indx}(ht_p) = extractBefore(char(H{i}{H_tm_indx}(ht_p)), '^');
%                     for ex_parts = 1:power-1
%                         H{i}{H_tm_indx} = [H{i}{H_tm_indx}, H{i}{H_tm_indx}(ht_p)];
%                     end
%                 end
%             end
% 
%             % -- the scalar part of Ht is always 1 --
%             H{i}{H_tm_indx} = [str2sym('1'), H{i}{H_tm_indx}];
% 
%             H_tm_indx = H_tm_indx + 1;
%         end

        fprintf(". ");
    end
end