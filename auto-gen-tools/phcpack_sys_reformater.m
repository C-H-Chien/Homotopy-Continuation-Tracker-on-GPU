% -- transform polynomial systems into  --
clear;
format long;

% -- define directories --
% -- top priority --
problem = 'd1/';
category = 'benchmark-problems/';
fileFolder = '/home/chchien/BrownU/research/22-issac/';
target_folder = 'phcpack/';
outputFileName_phcpack_sys = 'phcpack-sys-formulations.txt';
fullOutputFileName_phcpck_sys = fullfile(fileFolder, category, problem, target_folder, outputFileName_phcpack_sys);
outputFileWr_phcpack_sys = fopen(fullOutputFileName_phcpck_sys, 'w');

% -- define systems --
if strcmp(problem, 'katsura6/')
    [f, numOfVars] = sys_katsura6();
elseif strcmp(problem, 'katsura7/')
    [f, numOfVars] = sys_katsura7();
elseif strcmp(problem, 'katsura8/')
    [f, numOfVars] = sys_katsura8();
elseif strcmp(problem, 'katsura9/')
    [f, numOfVars] = sys_katsura9();
elseif strcmp(problem, 'katsura10/')
    [f, numOfVars] = sys_katsura10();
elseif strcmp(problem, 'cyclic7/')
    [f, numOfVars] = sys_cyclic7();
elseif strcmp(problem, 'cyclic8/')
    [f, numOfVars] = sys_cyclic8();
elseif strcmp(problem, 'cyclic10/')
    [f, numOfVars] = sys_cyclic10();
elseif strcmp(problem, 'cyclic11/')
    [f, numOfVars] = sys_cyclic11();
elseif strcmp(problem, 'alea6/')
    [f, numOfVars] = sys_alea6();
elseif strcmp(problem, 'eco12/')
    [f, numOfVars] = sys_eco12();
elseif strcmp(problem, 'd1/')
    [f, numOfVars] = sys_d1();
elseif strcmp(problem, 'game6two/')
    [f, numOfVars] = sys_game6two();
elseif strcmp(problem, 'game7two/')
    [f, numOfVars] = sys_game7two();
elseif strcmp(problem, 'pole28sys/')
    [f, numOfVars] = sys_pole28sys();
end

% -- initialize one row of a system matrix for phcpack --
sys_mat = strings(numOfVars+1, 1);

for i = 1:numOfVars
    str_vars = 'x';
    cat_str_vars = strcat(str_vars, num2str(i));
    X(i) = str2sym(cat_str_vars);
end

for p = 1:numOfVars
%for p = 1:1
    % -- 1) find the variables and coefficients of the system --
    [coefficient, variable] = coeffs(f(p), X);
    
    % -- 2) loop over all terms --
    for i = 1:size(variable, 2)
        % -- 2-1) initialize all elements in a row to 0 --
        for k = 1:numOfVars+1
            sys_mat(k,1) = '0';
        end
        
        % -- 2-2) store the coefficient first --
        sys_mat(1,1) = string(coefficient(i));
        
        % -- 2-3) split the unknowns by multiplier * -- 
        unknowns = strsplit(string(variable(i)), '*');
        
        % -- 2-4) loop over all unknowns --
        for j = 1:size(unknowns, 2)
            % -- check if it is a scalar --
            
            % -- check the power of each unknown --
            if contains(unknowns(j), '^')
                degree = extractAfter(unknowns(j), '^');
                x_idx = extractBetween(unknowns(j), 'x', '^');
                sys_mat(str2double(x_idx)+1, 1) = degree;
            else
                x_idx = extractAfter(unknowns(j), 'x');
                sys_mat(str2double(x_idx)+1, 1) = '1';
            end
        end
        
        % -- 2-5) write the result to the file --
        for k = 1:numOfVars+1
            fprintf(outputFileWr_phcpack_sys, sys_mat(k,1));
            fprintf(outputFileWr_phcpack_sys, " ");
        end
        fprintf(outputFileWr_phcpack_sys, ';\n');
    end
    
    % -- 3) write the sperator to the file --
    for k = 1:numOfVars+1
        fprintf(outputFileWr_phcpack_sys, '0');
        fprintf(outputFileWr_phcpack_sys, " ");
    end
    fprintf(outputFileWr_phcpack_sys, ';\n');
    
end


% -- close the files --
fclose(outputFileWr_phcpack_sys);
