% -- Automatic constant index matrix generator --
% -- for Homotopy Continuation on GPU --
% -- Dates: Dec. 9th, 2021 --
% -- Author: Chiang-Heng Chien --
clear;

% -- output write file directory and file name --
fileFolder = '/home/chchien/BrownU/research/22-issac/';
%category = 'computer-vision-problems/';
category = 'benchmark-problems/';
problem = 'alea6-extend/';
wrFolder = 'index-matrices/';
wr2file = 1;

if strcmp(problem, 'katsura6/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_katsura6();
elseif strcmp(problem, 'katsura7/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_katsura7();
elseif strcmp(problem, 'katsura8/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_katsura8();
elseif strcmp(problem, 'katsura9/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_katsura9();
elseif strcmp(problem, 'katsura10/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_katsura10();
elseif strcmp(problem, 'katsura11/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_katsura11();
elseif strcmp(problem, 'katsura12/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_katsura12();
elseif strcmp(problem, 'katsura13/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_katsura13();
elseif strcmp(problem, 'katsura14/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_katsura14();
elseif strcmp(problem, 'katsura15/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_katsura15();
elseif strcmp(problem, 'katsura20/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_katsura20();
elseif strcmp(problem, 'katsura21/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_katsura21();
elseif strcmp(problem, 'd1/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_d1();
elseif strcmp(problem, 'alea6/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_alea6();
elseif strcmp(problem, 'alea6-extend/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_alea6_extend();
elseif strcmp(problem, 'cyclic7/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_cyclic7();
elseif strcmp(problem, 'cyclic8/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_cyclic8();
elseif strcmp(problem, 'cyclic9/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_cyclic9();
elseif strcmp(problem, 'eco12/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_eco12();
elseif strcmp(problem, 'game6two/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_game6two();
elseif strcmp(problem, 'game7two/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_game7two();
elseif strcmp(problem, 'pole28sys/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_pole28sys();
elseif strcmp(problem, '4vTrg/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_4vTrg();
elseif strcmp(problem, '3vTrg/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_3vTrg();
elseif strcmp(problem, '3vTrg_relax/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_3vTrg_relax();
elseif strcmp(problem, '5pt_rel_pose_w_depth_recon/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_5pt_rel_pose_w_depth_recon();
elseif strcmp(problem, 'PnP_wo_principal_point/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_PnP_wo_principal_point();
elseif strcmp(problem, 'optimalPnP_w_quaternion/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_optimalPnP_w_quaternion();
elseif strcmp(problem, '3pt_rel_pose_w_homo_constraint/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_3pt_rel_pose_w_homo_constraint();
elseif strcmp(problem, 'r6p/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_r6p();
elseif strcmp(problem, 'refractive_p5p/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_refractive_p5p();
elseif strcmp(problem, 'refractive_p6p/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_refractive_p6p();
elseif strcmp(problem, '3view_unknownf_pHC/')
    [numOfVars, numOfCoeff, X, C, D, J, Homotopy] = problems_3view_unknownf_pHC();
end

% -- seperate J into Hx and Ht --
partial_H_partial_x = J(1:numOfVars, 1:numOfVars);
partial_H_partial_t = J(1:numOfVars, numOfVars+1);

% -- strings of Hx and Ht --
org_Hx = strings(numOfVars, numOfVars);
org_Ht = strings(numOfVars, 1);
org_H = strings(numOfVars, 1);

% -- organize the system of Hx, Ht, and H --
fprintf("Organizing Hx, Ht, and H ");
for i = 1:size(partial_H_partial_x, 1)
    % -- Hx --
    for j = 1:size(partial_H_partial_x, 2)
        org_Hx(i,j) = organize_sys(partial_H_partial_x(i,j), 'Hx');
    end
    
    % -- Ht --
    org_Ht(i,1) = organize_sys(partial_H_partial_t(i,1), 'Ht');
    
    % -- H --
    org_H(i,1) = organize_sys(Homotopy(i,1), 'Hx');
    fprintf(". ");
end
fprintf("\n");

% -- construct parts in each term of Hx --
fprintf("Constructing parts in each term of Hx and Ht ");
[Hx, Ht] = split_parts_in_each_term(org_Hx, org_Ht, '1');
fprintf("\n");

% -- find the maximal terms and maximal number of parts --
[Hx_maximal_terms, Hx_maximal_parts, Ht_maximal_terms, Ht_maximal_parts] = find_maximal_terms_and_parts(Hx, Ht);

% -- generate the constant matrices for Hx --
fprintf("Generating constant matrices ");
[const_matrix_Hx, const_matrix_Ht] = generate_constant_matrices_improve(numOfVars, numOfCoeff, Hx, Ht, Hx_maximal_terms, Hx_maximal_parts, Ht_maximal_terms, Ht_maximal_parts);

% -- write the constant matrices to txt files --
if wr2file == 1
    constant_mat_write2file_improve(fileFolder, category, problem, wrFolder, const_matrix_Hx, const_matrix_Ht);
end

fprintf("It's finished!\n\n");
fprintf("The maximal terms and parts are:\n");
fprintf("1) Hx.maximal_terms = %d\n", Hx_maximal_terms);
fprintf("2) Hx.maximal_parts = %d\n", Hx_maximal_parts);
fprintf("3) Ht.maximal_terms = %d\n", Ht_maximal_terms);
fprintf("4) Ht.maximal_parts = %d\n", Ht_maximal_parts);
