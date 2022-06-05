% -- params input initializer --
% for i = 0:47
%     pstr = strcat('magmaFloatComplex', {' '}, 'x', num2str(i+6), {' '}, '=', {' '}, 'h_target_params[', num2str(i), '];');
%     fprintf(string(pstr));
%     fprintf('\n');
% end

% -- params convert to coefficients --
fileFolder = '/home/chchien/BrownU/research/22-cvpr/self-minimal_problems/';
problem = 'trifocal_rel_pose_f/';

inputCoefsFileName = 'params2coeffs.txt';
fullInputFileName = fullfile(fileFolder, problem, inputCoefsFileName);
coefsFileRd = fopen(fullInputFileName, 'r');

ldata = textscan(coefsFileRd, '%s', 'Delimiter', ',', 'whitespace', '{}', 'CollectOutput', true  );
line = ldata{1};

coeffs_symbol = 'x';

for i = 1:size(line, 1)
    split_equality = strsplit(string(line{i,1}), '=');
    
    % -- 1) replace the lhs --
    lhs = split_equality(1,1);
    rep_lhs = strcat('h_target_coeffs[', num2str(i-1), ']', {' '});
    new_lhs_str = strrep(string(line{i,1}), lhs, rep_lhs);
    
    % -- 2) replace '^2' by double multiplication in the rhs --
    rhs = split_equality(1,2);
    rhs_terms = strsplit(rhs, {'+', '-', '*'});
    update_rhs = rhs;
    for j = 1:size(rhs_terms, 2)
        if contains(rhs_terms(1,j), '^2')
            x_indx = extractBetween(rhs_terms(1,j), coeffs_symbol, '^2');
            x_term_square_rep = strcat(coeffs_symbol, x_indx, '*', coeffs_symbol, x_indx);
            x_term_square_org = strcat(coeffs_symbol, x_indx, '^2');
            update_rhs = strrep(update_rhs, x_term_square_org, x_term_square_rep);
        elseif contains(rhs_terms(1,j), '^3')
            x_indx = extractBetween(rhs_terms(1,j), coeffs_symbol, '^3');
            x_term_square_rep = strcat(coeffs_symbol, x_indx, '*', coeffs_symbol, x_indx, '*', coeffs_symbol, x_indx);
            x_term_square_org = strcat(coeffs_symbol, x_indx, '^3');
            update_rhs = strrep(update_rhs, x_term_square_org, x_term_square_rep);
        end
    end
    
    new_full_str = strcat(string(rep_lhs), '=', update_rhs, ';\n');
    fprintf(new_full_str);
end