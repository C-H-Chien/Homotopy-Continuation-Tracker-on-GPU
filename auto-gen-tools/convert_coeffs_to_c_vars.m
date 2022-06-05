% -- generate  --
clear;
close all;

% -- define directories --
%fileFolder = '/home/chchien/BrownU/research/22-issac/computer-vision-problems/';
fileFolder = '/home/chchien/BrownU/research/22-issac/benchmark-problems/';
problem = 'alea6-extend/';
numOfVars = 6;

% -- input file --
inputName_raw_poly = 'polynomial-formulations-raw';
fullFileName_raw_poly = fullfile(fileFolder, problem, inputName_raw_poly);
inputFileWr_raw_poly = fopen(fullFileName_raw_poly, 'r');

% -- output file --
outputFileName_rep_poly = 'polynomial-formulations-with-cs.txt';
fullOutputFileName_rep_poly = fullfile(fileFolder, problem, outputFileName_rep_poly);
outputFileWr_rep_poly = fopen(fullOutputFileName_rep_poly, 'w');

% -- 1) read data from file --
ldata = textscan(inputFileWr_raw_poly, '%s', 'Delimiter', ',', 'whitespace', '{}', 'CollectOutput', true);
line_raw_poly = ldata{1};
str_raw_poly = string(line_raw_poly);

% -- 2) loop over all polynomials
for i = 1:size(line_raw_poly, 1)
    % -- 2-1) splite string by arithmetic operators --
    [p_terms, arith_sym] = strsplit(str_raw_poly(i,1), {'+', '-'});
    
    % -- for each term --
    cat_terms = '';
    for j = 1:size(p_terms, 2)
        % -- 2-2) split each term by multiplier to get parts --
        parts = strsplit(p_terms(1,j), '*');
        
        % -- for each part --
        cat_parts = '';
        for k = 1:size(parts, 2)
            % -- 2-3) extract the index of x to see whether it is a
            % parameter --
            if contains(parts(1,k), 'x')
                % -- check whether the degree is greater than 1 --
                if contains(parts(1,k), '^')
                    x_idx = extractBetween(parts(1,k), 'x', '^');
                    degree = extractAfter(parts(1,k), '^');
                    if str2double(x_idx) > numOfVars
                        shifted_idx = str2double(x_idx) - numOfVars;
                        target_sym = 'c';
                    else
                        shifted_idx = str2double(x_idx);
                        target_sym = 'x';
                    end
                    % -- concatenate the varaibles up to its degree --
                    cat_sym_mul = '';
                    for d = 1:str2double(degree)
                        if d < str2double(degree)
                            cat_sym_mul = strcat(cat_sym_mul, target_sym, string(shifted_idx), '*');
                        else
                            cat_sym_mul = strcat(cat_sym_mul, target_sym, string(shifted_idx));
                        end
                    end
                    rep_part = strcat({' '}, cat_sym_mul);
                else
                    x_idx = extractAfter(parts(1,k), 'x');
                    if str2double(x_idx) > numOfVars
                        shifted_idx = str2double(x_idx) - numOfVars;
                        rep_part = strcat('c', string(shifted_idx));
                    else
                        rep_part = parts(1,k);
                    end
                end
            else
                rep_part = parts(1,k);
            end
            
            cat_parts = strcat(cat_parts, rep_part, '*');
        end
        cat_parts = char(cat_parts);
        cat_parts = cat_parts(1:end-1);
        
        % -- 3-3) cacetenate all parts with original arithmatic symbols --
        if j == size(p_terms, 2)
            cat_terms = strcat(cat_terms, cat_parts);
        else
            cat_terms = strcat(cat_terms, cat_parts, {' '}, arith_sym(1,j));
        end
    end
    
    fprintf(outputFileWr_rep_poly, cat_terms);
    fprintf(outputFileWr_rep_poly, '\n');
end

fclose(outputFileWr_rep_poly);
