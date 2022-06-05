% -- generate  --
clear;

%sympref('FloatingPointOutput',true);

% -- define directories --
fileFolder = '/home/chchien/BrownU/research/22-issac/parametric-HC/';
problem = 'chicago/';
numOfParams = 56;

% -- input file --
% -- M2 ccode evaluations --
inputName_eval = 'cpu-eval-HxHt';
fullInputFileName_eval = fullfile(fileFolder, problem, inputName_eval);
inputFileWr_eval = fopen(fullInputFileName_eval, 'r');
% -- output, convert to sym --
inputName_out = 'cpu-eval-Hx-out';
fullInputFileName_out = fullfile(fileFolder, problem, inputName_out);
inputFileWr_out = fopen(fullInputFileName_out, 'r');

% -- output files --
outputFileName_organized_eval = 'organized_eval_HxHt.txt';
fullOutputFileName_organized_eval = fullfile(fileFolder, problem, outputFileName_organized_eval);
outputFileWr_organized_eval = fopen(fullOutputFileName_organized_eval, 'w');
outputFileName_organized_out = 'organized_eval_sym_out.txt';
fullOutputFileName_organized_out = fullfile(fileFolder, problem, outputFileName_organized_out);
outputFileWr_organized_out = fopen(fullOutputFileName_organized_out, 'w');

% -- 1) extract data from file --
ldata = textscan(inputFileWr_eval, '%s', 'Delimiter', ',', 'whitespace', '{}', 'CollectOutput', true);
line_eval = ldata{1};
str_line_eval = string(line_eval);
ldata = textscan(inputFileWr_out, '%s', 'Delimiter', ',', 'whitespace', '{}', 'CollectOutput', true);
line_out = ldata{1};
str_line_out = string(line_out);

% -- 2) store left hand side and right hand side of each equation --
lhs = strings(size(str_line_eval, 1), 1);
rhs = strings(size(str_line_eval, 1), 1);
for i = 1:size(str_line_eval, 1)
    split_equality = strsplit(str_line_eval(i,1), '=');
    
    lhs(i,1) = split_equality(1,1);
    rhs(i,1) = split_equality(1,2);
    
    % -- 2-1) replace semicolon with a white space --
    rhs(i,1) = strrep(rhs(i,1), ';', {' '});
    
    % -- 2-2) padd a white space before each lhs --
    lhs(i,1) = strcat({' '}, lhs(i,1));
end

% -- 3) loop over all rhs and do replacements --
for i = 1:size(str_line_eval, 1)
    % -- 3-1) split the rhs by arithmetic operator --
    rhs_parts = strsplit(rhs(i,1), {'*', '+'});
    
    % -- loop over all parts of rhs --
    for j = 1:size(rhs_parts, 2)
        % -- if the target rhs part contains G which is required to be replaced --
        if contains(rhs_parts(1,j), 'G')
            % -- loop over all lhs --
            for k = 1:size(lhs, 1)
                if contains(rhs_parts(1,j), lhs(k,1))
                    % -- replacing --
                    target_rep = strcat('(', rhs(k,1), ')');
                    rhs(i,1) = strrep(rhs(i,1), rhs_parts(1,j), target_rep);
                    break;
                end
            end
        else
            continue;
        end
    end
    
%     % -- write the replaced string rhs to the file --
%     str_wr = strcat(lhs(i,1), {' '}, '=', {' '}, rhs(i,1));
%     fprintf(outputFileWr_organized_eval, str_wr);
%     fprintf(outputFileWr_organized_eval, '\n');
    
%     % -- write the symbolic rhs to the file --
%     sym_wr = str2sym(str_wr);
%     fprintf(outputFileWr_sym_org_eval, string(sym_wr));
%     fprintf(outputFileWr_sym_org_eval, '\n');
    
    % -- monitor the progress --
    if mod(i,300) == 0
        fprintf(' .');
    end
end
fprintf('\n');

% -- 4) replace C0 with (-1) and write the results to the file --
for i = 1:size(str_line_eval, 1)
    wr_eq = strrep(rhs(i,1), "C0", "(-1)");
    str_wr = strcat(lhs(i,1), {' '}, '=', {' '}, wr_eq);
    fprintf(outputFileWr_organized_eval, str_wr);
    fprintf(outputFileWr_organized_eval, '\n');
end

% -- 5) extract rhs of the outputs --
rhs_out = strings(size(str_line_out, 1), 1);
for i = 1:size(str_line_out, 1)
    split_equality = strsplit(str_line_out(i,1), '=');
    
    rhs_out(i,1) = split_equality(1,2);
    rhs_out(i,1) = strrep(rhs_out(i,1), " ", '');
end

% -- 6) convert to syms for only the outputs --
for j = 1:size(rhs_out, 1)
    if strcmp(rhs_out(j,1), '0')
        fprintf(outputFileWr_organized_out, 'MAGMA_C_ZERO');
        fprintf(outputFileWr_organized_out, '\n\n');
        continue;
    end
    
    for i = 1:size(str_line_eval, 1)
        cmp_lhs = strrep(lhs(i,1), " ", '');
        if strcmp(cmp_lhs, rhs_out(j,1))
            rhs_str = strrep(rhs(i,1), "C0", "(-1)");
            wr_sym = str2sym(rhs_str);
            fprintf(outputFileWr_organized_out, string(wr_sym));
            fprintf(outputFileWr_organized_out, '\n\n');
            break;
        end
    end
end

fclose(outputFileWr_organized_eval);
fclose(outputFileWr_organized_out);