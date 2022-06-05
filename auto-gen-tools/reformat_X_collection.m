% -- params input initializer --
% for i = 0:47
%     pstr = strcat('magmaFloatComplex', {' '}, 'x', num2str(i+6), {' '}, '=', {' '}, 'h_target_params[', num2str(i), '];');
%     fprintf(string(pstr));
%     fprintf('\n');
% end

% -- params convert to coefficients --
fileFolder = '/home/chchien/BrownU/research/22-cvpr/self-minimal_problems/';
problem = '3view_unknownf/';

inputCoefsFileName = 'X4';
outputSolName = 'new_X4';
fullInputFileName = fullfile(fileFolder, problem, inputCoefsFileName);
fullOutputFileName = fullfile(fileFolder, outputSolName);
coefsFileRd = fopen(fullInputFileName, 'r');
newXFileWr = fopen(fullOutputFileName, 'w');

ldata = textscan(coefsFileRd, '%s', 'Delimiter', ',', 'whitespace', '{}', 'CollectOutput', true  );
line = ldata{1};

coeffs_symbol = 'x';

for i = 1:size(line, 1)
    split_equality = strsplit(string(line{i,1}), '=');
    
    % -- 1) replace the lhs --
    lhs = split_equality(1,1);
    rep_lhs = strcat('X[', num2str(i-1), ']', {' '});
    new_lhs_str = strrep(string(line{i,1}), lhs, rep_lhs);

    fprintf(newXFileWr, new_lhs_str);
    fprintf(newXFileWr, '\n');
end

fclose(newXFileWr);