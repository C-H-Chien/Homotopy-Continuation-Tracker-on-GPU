sympref('FloatingPointOutput',false);
clear;

fileFolder = '/home/chchien/BrownU/research/';
%problem_name = '3vTrg/coefficient-HxHtH/v2-p2c/';
problem_name = 'alea6-extend/';
%category = 'computer-vision-problems';
category = 'benchmark-problems';

if strcmp(category, 'benchmark-problems')
    inputSolName = strcat('22-issac/benchmark-problems/', problem_name,'julia-start-sols-raw');
    outputSolName = strcat('22-issac/benchmark-problems/', problem_name, 'start-sols.txt');
else
    inputSolName = strcat('22-issac/computer-vision-problems/', problem_name,'julia-start-sols-raw');
    outputSolName = strcat('22-issac/computer-vision-problems/', problem_name, 'start-sols.txt');
end

fullInputFileName = fullfile(fileFolder, inputSolName);
fullOutputFileName = fullfile(fileFolder, outputSolName);

solFileRd = fopen(fullInputFileName, 'r');
solFileWr = fopen(fullOutputFileName, 'w');

result = [];
ldata = textscan(solFileRd, '%s', 'Delimiter', {'im', ','});
line = ldata{1,1};

real_parts = strings(size(line, 1), 1);
imag_parts = strings(size(line, 1), 1);

indx = 1;
for i = 1:size(line, 1)
    
    cvt_str = string(line{i,1});
    clean_str = strrep(cvt_str, '[', '');
    
    cplx_str = strcat(clean_str,'i');
    cplx_num = double(cplx_str);
    
    real_parts(indx,1) = num2str(real(cplx_num), 15);
    imag_parts(indx,1) = num2str(imag(cplx_num), 15);
    indx = indx + 1;
    wr_str = strcat(real_parts(i,1), '\t', imag_parts(i,1), '\n');
    fprintf(solFileWr, wr_str);
    
    if mod(i,1000) == 0
        fprintf('. ');
    end
end
fprintf('\n');

fclose(solFileWr);