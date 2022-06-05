sympref('FloatingPointOutput',false);

fileFolder = '/home/chchien/BrownU/research/22-issac/';
category = 'benchmark-problems/';
problem = 'katsura9/';
inputSolName = 'julia-start-sols-raw';
outputSol_bertini = 'bertini-start-sols.txt';

numOfSols = 512;
numOfvars = 10;

fullInputFileName_sols = fullfile(fileFolder, category, problem, inputSolName);
fullOutputFileName_phcjulia = fullfile(fileFolder, category, problem, outputSol_bertini);

solFileRd_sols = fopen(fullInputFileName_sols, 'r');
solFileWr_sols = fopen(fullOutputFileName_phcjulia, 'w');

result = [];
ldata_sols = textscan(solFileRd_sols, '%s', 'Delimiter', {'im', ','});
line = ldata_sols{1,1};

real_parts = strings(size(line, 1), 1);
imag_parts = strings(size(line, 1), 1);

fprintf(solFileWr_sols, num2str(numOfSols));
fprintf(solFileWr_sols, '\n');
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
    fprintf(solFileWr_sols, wr_str);
    if mod(i,numOfvars) == 0
        fprintf(solFileWr_sols, '\n');
    end
    
    if mod(i,1000) == 0
        fprintf('. ');
    end
end
fprintf('\n');

fclose(solFileWr_sols);