
% -- define directories --
fileFolder = '/home/chchien/BrownU/research/22-issac/benchmark-problems/';
problem = 'katsura7/';
inputSolName = 'M2-start-sols-raw';
outputSolName = 'output-start-sols.txt';

fullInputFileName = fullfile(fileFolder, problem, inputSolName);
fullOutputFileName = fullfile(fileFolder, problem, outputSolName);

solFileRd = fopen(fullInputFileName, 'r');
solFileWr = fopen(fullOutputFileName, 'w');

% -- define parameters --
changeLine = 0;
numOfVar = 18;
varCount = 0;
solCount = 0;

% -- start reformating --
result = [];
if strcmp(problem, 'chicago14/') || strcmp(problem, 'cleveland14/')
    ldata = textscan(solFileRd, '%f%f', 'Delimiter', ',', 'whitespace', '{}', 'CollectOutput', true  );
    line = ldata{1};
    soluMat = reshape(line, [], 2);
    result = [result soluMat];
    
    for row = 1:size(result, 1)
        % -- for the use of GPU --
        fprintf(solFileWr, '%f\t', result(row, 1));
        fprintf(solFileWr, '%f\n', result(row, 2));
        
        % -- change to the next lin in the write txt file --
        varCount = varCount + 1;
        if changeLine == 1
            if varCount == numOfVar
                fprintf(solFileWr, '\n');
                varCount = 0;
                solCount = solCount + 1;
            end
        end
    end
else
    ldata = textscan(solFileRd, '%f%f%s %f%f%s %f%f%s %f%f%s', 'Delimiter', ',', 'whitespace', '{}', 'CollectOutput', true  );
    % -- read data from the solution file and reofrmate as a matrix --
    for i = 1:2:8
        line = ldata{i};
        soluMat = reshape(line, [], 2);
        result = [result soluMat];
    end

    % -- write data from soluMat to the output solution file --
    for row = 1:size(result, 1)
        for col = 1:2:8
            % -- for the use of GPU --
            fprintf(solFileWr, '%f\t', result(row, col));
            fprintf(solFileWr, '%f\n', result(row, col+1));

            % -- for the use of CPU --
            %fprintf(solFileWr, '{%f,\t', result(row, col));
            %fprintf(solFileWr, '%f},\n', result(row, col+1));

            % -- change to the next lin in the write txt file --
            varCount = varCount + 1;
            if changeLine == 1
                if varCount == numOfVar
                    fprintf(solFileWr, '\n');
                    varCount = 0;
                    solCount = solCount + 1;
                end
            end
        end
    end
end

% -- make sure that the # of written solutions are correct --
%fprintf('%d\n', solCount);
fclose(solFileWr);
