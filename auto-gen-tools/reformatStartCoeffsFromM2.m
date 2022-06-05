% -- reformate coefficient files used to be read by the script --
% -- define directories --
fileFolder = '/home/chchien/BrownU/research/22-cvpr/self-minimal_problems/';
problem = 'R6P1lin/';
inputCoefsFileName = 'M2-start-coeffs-raw';
outputCoefsFileName = 'output-start-coeffs.txt';

% -- define parameters --
changeLine = 0;
numOfVar = 18;
varCount = 0;
solCount = 0;

fullInputFileName = fullfile(fileFolder, problem, inputCoefsFileName);
fullOutputFileName = fullfile(fileFolder, problem, outputCoefsFileName);

coefsFileRd = fopen(fullInputFileName, 'r');
coefsFileWr = fopen(fullOutputFileName, 'w');

result = [];
if strcmp(problem, 'chicago14/') || strcmp(problem, 'cleveland14/')
    ldata = textscan(coefsFileRd, '%f%f', 'Delimiter', ',', 'whitespace', '{}', 'CollectOutput', true  );
    line = ldata{1};
    soluMat = reshape(line, [], 2);
    result = [result soluMat];
    
    for row = 1:size(result, 1)
        % -- for the use of GPU --
        fprintf(coefsFileWr, '%f\t', result(row, 1));
        fprintf(coefsFileWr, '%f\n', result(row, 2));
        
        % -- change to the next lin in the write txt file --
        varCount = varCount + 1;
        if changeLine == 1
            if varCount == numOfVar
                fprintf(coefsFileWr, '\n');
                varCount = 0;
                solCount = solCount + 1;
            end
        end
    end
else
    ldata = textscan(coefsFileRd, '%f%f%s %f%f%s %f%f%s %f%f%s', 'Delimiter', ',', 'whitespace', '{}', 'CollectOutput', true );

    % -- read data from the solution file and reofrmate as a matrix --
    for i = 1:2:8
        line = ldata{i};
        soluMat = reshape(line, [], 2);
        result = [result soluMat];
    end

    % -- write data from soluMat to the output solution file --
    for row = 1:size(result, 1)
        for col = 1:2:8
            % -- format for the use of GPU --
            fprintf(coefsFileWr, '%20f\t', result(row, col));
            fprintf(coefsFileWr, '%20f\n', result(row, col+1));

            % -- format for the use of CPU --
            %fprintf(coefsFileWr, '{%10f,\t', result(row, col));
            %fprintf(coefsFileWr, '%10f},\n', result(row, col+1));

            varCount = varCount + 1;
            if changeLine == 1
                if varCount == numOfVar
                    fprintf(coefsFileWr, '\n');
                    varCount = 0;
                    solCount = solCount + 1;
                end
            end
        end
    end
end

fprintf('%d\n', solCount);
fclose(coefsFileWr);