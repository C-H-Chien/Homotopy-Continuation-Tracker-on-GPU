function [Hx_maximal_terms, Hx_maximal_parts, Ht_maximal_terms, Ht_maximal_parts] = find_maximal_terms_and_parts(Hx, Ht)
    
    % -- initialization --
    Hx_maximal_terms = 0;
    Hx_maximal_parts = 0;
    Ht_maximal_terms = 0;
    Ht_maximal_parts = 0;
%     H_maximal_terms = 0;
%     H_maximal_parts = 0;
    
    for i = 1:size(Hx, 1)
        for j = 1:size(Hx, 2)
            % -- obtain maximal terms of Hx --
            if size(Hx{i,j}, 2) > Hx_maximal_terms
                Hx_maximal_terms = size(Hx{i,j}, 2);
            end

            % -- obtain maximal number of parts of Hx --
            for k = 1:size(Hx{i,j}, 2)

                if size(Hx{i,j}{k}, 2) > Hx_maximal_parts 
                    Hx_maximal_parts = size(Hx{i,j}{k}, 2);
                end
            end
        end

        % -- obtain maximal terms of Ht --
        if size(Ht{i}, 2) > Ht_maximal_terms
            Ht_maximal_terms = size(Ht{i}, 2);
        end

        % -- obtain maximal number of parts of Ht --
        for m = 1:size(Ht{i}, 2)
            if size(Ht{i}{m}, 2) > Ht_maximal_parts
                tmp = Ht{i};
                Ht_maximal_parts = size(Ht{i}{m}, 2);
            end
        end

%         % -- obtain maximal terms of H --
%         if size(H{i}, 2) > H_maximal_terms
%             H_maximal_terms = size(H{i}, 2);
%         end
% 
%         % -- obtain maximal number of parts of H --
%         for m = 1:size(H{i}, 2)
%             if size(H{i}{m}, 2) > H_maximal_parts
%                 tmp = H{i};
%                 H_maximal_parts = size(H{i}{m}, 2);
%             end
%         end
    end

    % -- because cd term is considered two parts so the Ht_maximal_parts should
    % minus 1 --
    Ht_maximal_parts = Ht_maximal_parts - 1;

end