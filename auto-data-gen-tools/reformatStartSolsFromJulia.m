sympref('FloatingPointOutput',false);
clear;

fileFolder    = "/your/path/to/the/start/solutions/generated/by/Julia/";
inputSolName  = "julia-start-sols-raw";
outputSolName = "start-sols.txt";
fullInputFileName  = fullfile(fileFolder, inputSolName);
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
