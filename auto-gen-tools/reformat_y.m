
% -- 1) read only-y files --
% -- 1-1) define directories --
fileFolder = '/home/chchien/BrownU/research/22-issac/';
%category = 'computer-vision-problems/';
category = 'benchmark-problems/';
problem = 'katsura21/';
%problem = '3view_unknownf_pHC/coefficient-HxHtH/';
inputFileName = 'M2-HxH-y-raw';

% -- 1-2) define parameters --
numOfVars = 22;

% -- 1-3) read from directory --
fullInputFileName = fullfile(fileFolder, category, problem, inputFileName);
FileRd = fopen(fullInputFileName, 'r');
ldata = textscan(FileRd, '%s', 'Delimiter', ',', 'whitespace', '{}', 'CollectOutput', true );
udata = ldata;

% -- 2) write to the file --
% -- 2-1) define write directory and file name --
if contains(inputFileName, 'HxH-')
    outH = 'HxH-';
else
    outH = 'HxHt-';
end
outputFileName = strcat('output-', outH, 'y-for-cpu');
fullOutputFileName = fullfile(fileFolder, category, problem, outputFileName);
FileWr = fopen(fullOutputFileName, 'w');

% -- 3) loop over all data and construct the matrix --
% -- 3-1) r_cgesvA --
matrix_y_A = strings(numOfVars);
row_y = 1;
col_y = 1;
for i = 1:numOfVars*numOfVars
    Gvar = string(extractBetween(ldata{1,1}{i}, '= ', ';'));
    matrix_y_A(row_y, col_y) = Gvar;
    row_y = row_y + 1;
    if mod(i, numOfVars) == 0
        row_y = 1;
        col_y = col_y + 1;
    end
end

% -- 3-2) r_cgesvB --
matrix_y_b = strings(numOfVars, 1);
y_b_indx = 1;
for i = numOfVars*numOfVars+1:size(ldata{1,1}, 1)
    Gvar = string(extractBetween(ldata{1,1}{i}, '= ', ';'));
    matrix_y_b(y_b_indx, 1) = Gvar;
    y_b_indx = y_b_indx + 1;
end

% -- 4) write y-matrix to the file --
if contains(inputFileName, 'HxH-')
    rhs_sign = '] = ';
else
    rhs_sign = '] = -';
end
for i = 1:numOfVars
    case_indx = num2str(i-1);
    for j = 1:numOfVars
        str_A = strcat('r_cgesvA[', num2str(j-1 + (i-1)*numOfVars) ,'] =', {' '}, matrix_y_A(j,i), ';\n');
        fprintf(FileWr, str_A);
    end
    
    str_b = strcat('r_cgesvB[', num2str(i-1), rhs_sign, matrix_y_b(i,1), ';\n');
    fprintf(FileWr, str_b);
    fprintf(FileWr, '\n');
end

fclose(FileWr);