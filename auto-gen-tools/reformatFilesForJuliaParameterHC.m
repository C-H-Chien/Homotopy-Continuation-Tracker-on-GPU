%sympref('FloatingPointOutput',false);

fileFolder = '/home/chchien/BrownU/research/22-issac/';
category = 'computer-vision-problems/';
problem = '5pt_rel_pose_w_depth_recon/';
subfolder = 'phc-julia/';
inputSolName = 'julia-start-sols-raw';
%inputSolName = 'start-sols_julia.txt';
inputStartParams = 'start_params.txt';
inputTargetParams = 'target_params.txt';
outputSol_phcjulia = 'phcjulia-start-sols.txt';
outputSol_start_params = 'phcjulia-start-params.txt';
outputSol_target_params = 'phcjulia-target-params.txt';

numOfVars = 16;

fullInputFileName_sols = fullfile(fileFolder, category, problem, inputSolName);
fullOutputFileName_phcjulia = fullfile(fileFolder, category, problem, subfolder, outputSol_phcjulia);

fullInputFileName_start_params = fullfile(fileFolder, category, problem, inputStartParams);
fullOutputFileName_start_params = fullfile(fileFolder, category, problem, subfolder, outputSol_start_params);

fullInputFileName_target_params = fullfile(fileFolder, category, problem, inputTargetParams);
fullOutputFileName_target_params = fullfile(fileFolder, category, problem, subfolder, outputSol_target_params);

solFileRd_sols = fopen(fullInputFileName_sols, 'r');
solFileWr_sols = fopen(fullOutputFileName_phcjulia, 'w');

solFileRd_start_params = fopen(fullInputFileName_start_params, 'r');
solFileWr_start_params = fopen(fullOutputFileName_start_params, 'w');

solFileRd_target_params = fopen(fullInputFileName_target_params, 'r');
solFileWr_target_params = fopen(fullOutputFileName_target_params, 'w');

result = [];
ldata_sols = textscan(solFileRd_sols, '%s', 'Delimiter', {'im', ','});
line_sols = ldata_sols{1,1};
ldata_start_params = textscan(solFileRd_start_params, '%s');
line_start_params = ldata_start_params{1,1};
ldata_target_params = textscan(solFileRd_target_params, '%s');
line_target_params = ldata_target_params{1,1};

fprintf(solFileWr_sols, 's1=[');
indx = 1;
solution_idx = 2;
for i = 1:size(line_sols, 1)

    % -- reformate start solutions --
    cvt_str = string(line_sols{i,1});
    clean_str = strrep(cvt_str, '[', '');
    
    wr_str = strcat(clean_str, 'im');    
    fprintf(solFileWr_sols, wr_str);
    
    if indx == numOfVars
        %wr_str = strcat("],", {' '}, "[");
        wr_str = strcat("];\ns", num2str(solution_idx), "=[");
        fprintf(solFileWr_sols, wr_str);
        indx = 0;
        solution_idx = solution_idx + 1;
    else
        wr_str = strcat(",", {' '});
        fprintf(solFileWr_sols, wr_str);
    end
    indx = indx + 1;   
    
    if mod(i,1000) == 0
        fprintf('. ');
    end
end
%fprintf(solFileWr_sols, ']\n');
fprintf('\n');

fprintf(solFileWr_start_params, '[');
fprintf(solFileWr_target_params, '[');
for i = 1:2:size(line_start_params, 1)
    % -- reformate start params --
    str_real_part = string(line_start_params{i,1});
    str_imag_part = string(line_start_params{i+1,1});
    wr_str = strcat(str_real_part, "+", str_imag_part, "im, ");
    fprintf(solFileWr_start_params, wr_str);
    
    % -- reformate target params --
    str_real_part = string(line_target_params{i,1});
    str_imag_part = string(line_target_params{i+1,1});
    wr_str = strcat(str_real_part, "+", str_imag_part, "im, ");
    fprintf(solFileWr_target_params, wr_str);
end
fprintf(solFileWr_start_params, ']');
fprintf(solFileWr_target_params, ']');

fclose(solFileWr_sols);
fclose(solFileWr_start_params);
fclose(solFileWr_target_params);