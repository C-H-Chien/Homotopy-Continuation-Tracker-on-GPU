%% -- reformate the evaluation HxHt used in magmahc: 1st optimization --
fileFolder = '/home/chchien/BrownU/research/22-issac/';
%category = 'computer-vision-problems/';
category = 'benchmark-problems/';
problem = 'katsura21/';
%problem = '3view_unknownf_pHC/coefficient-HxHtH/';

% -- problem parameters --
numOfVar = 22;
numOfCoef = 3;

inputEvalFileName = 'M2-HxH-only-G-raw';
outputEvalFileName = 'output_HxH_for_cpu.txt';

fullInputFileName = fullfile(fileFolder, category, problem, inputEvalFileName);
fullOutputFileName = fullfile(fileFolder, category, problem, outputEvalFileName);

evalFileRd = fopen(fullInputFileName, 'r');
evalFileWr = fopen(fullOutputFileName, 'w');

ldata = textscan(evalFileRd, '%s', 'Delimiter', ';', 'CollectOutput', true);
line = string(ldata{1});

trackStr = '';
startCoefStr = '';
targetCoefStr = '';
xStr = '';
G_Str = '';
% -- strings of r_track registers --
for i = 1:numOfVar
    catStr_track = strcat({'s_track['}, string(i-1), {']'});
    if i > 1
        trackStr = [trackStr catStr_track];
    else
        trackStr = catStr_track;
    end
end
trackStr = [trackStr 't'];

% -- strings of start coefficients registers --
for i = 1:numOfCoef
    catStr_startCoef = strcat({'s_startCoefs['}, string(i-1), {']'});
    if i > 1
        startCoefStr = [startCoefStr catStr_startCoef];
    else
        startCoefStr = catStr_startCoef;
    end
end

% -- strings of target coefficients registers --
for i = 1:numOfCoef
    catStr_targetCoef = strcat({'s_targetCoefs['}, string(i-1), {']'});
    if i > 1
        targetCoefStr = [targetCoefStr catStr_targetCoef];
    else
        targetCoefStr = catStr_targetCoef;
    end
end

% -- strings of X --
for i = 1:(numOfVar+1+numOfCoef+numOfCoef)
    catStr_x = strcat({'X'}, string(i-1));
    if i > 1
        xStr = [xStr catStr_x];
    else
        xStr = catStr_x;
    end
end

neg_sign = 'C1';
neg_lhs = '';
neg_rhs = '';
detect_neg = 0;
valid_Gs = [];
array_Gs = [];
G_lhs_indx = 0;
for i = 1:size(line, 1)
    splitLine = split(line(i));
    if size(splitLine, 1) == 6
        substituteStr_4 = splitLine(4);
        substituteStr_6 = splitLine(6);

%         % -- a negative sign --
%         if strcmp(splitLine(4), neg_sign)
%             neg_lhs = splitLine(2);
%             neg_rhs = splitLine(6);
%             detect_neg = 1;
%             continue;
%         end

        % -- first write the left hand side of the equation --
        fprintf(evalFileWr, 'magmaFloatComplex ');
        %catStr = strcat('G[', string(G_lhs_indx) , {'] '}, splitLine(3), {' '});
        catStr = strcat(splitLine(2), {' '}, splitLine(3), {' '});
        fprintf(evalFileWr, catStr);

        cmpStr_4 = splitLine(4);
        cmpStr_5 = splitLine(5);
        cmpStr_6 = splitLine(6);
        
        % -- extract the negative signed term --
        if detect_neg
            if strcmp(splitLine(4), neg_lhs)
                cmpStr_4 = splitLine(6);
                cmpStr_6 = neg_rhs;
            elseif strcmp(splitLine(6), neg_lhs)
                cmpStr_4 = splitLine(4);
                cmpStr_6 = neg_rhs;
            end

            cmpStr_5 = '-';
            detect_neg = 0;
        end
        substituteStr_4 = cmpStr_4;
        substituteStr_6 = cmpStr_6;

        % -- substitue coefficients and solution variables --
        for j = 1:(numOfVar+1+numOfCoef+numOfCoef)
            if j <= numOfVar+1
                if strcmp(cmpStr_6, xStr(j))
                    substituteStr_6 = trackStr(j);
                end
                if strcmp(cmpStr_4, xStr(j))
                    substituteStr_4 = trackStr(j);
                end
            elseif j > numOfVar+1 && j <= (numOfVar + 1 + numOfCoef)
                p = j - (numOfVar+1);
                if strcmp(cmpStr_6, xStr(j))
                    substituteStr_6 = startCoefStr(p);
                end
                if strcmp(cmpStr_4, xStr(j))
                    substituteStr_4 = startCoefStr(p);
                end
            elseif j > (numOfVar + 1 + numOfCoef)
                q = j - (numOfVar + 1 + numOfCoef);
                if strcmp(cmpStr_6, xStr(j))
                    substituteStr_6 = targetCoefStr(q);
                end
                if strcmp(cmpStr_4, xStr(j))
                    substituteStr_4 = targetCoefStr(q);
                end   
            end
        end
        catStr_rhs = strcat(substituteStr_4, {' '}, cmpStr_5, {' '}, substituteStr_6, {';\n'});
        fprintf(evalFileWr, catStr_rhs);
        
        catStr = strcat('G[', string(G_lhs_indx) , {']'});
        valid_Gs = [valid_Gs splitLine(2)];
        array_Gs = [array_Gs catStr];
        G_lhs_indx = G_lhs_indx + 1;
        
    elseif size(splitLine, 1) == 8
        catStr_rhs = strcat(splitLine(4), splitLine(5), splitLine(6), splitLine(7), splitLine(8), {';\n'});
        fprintf(evalFileWr, catStr_rhs);
        fprintf("G%d\n", i-1);
    else
        strShow = sprintf('Term %d has more than two arithmetic computations.\n', i-1);
        disp(strShow);
    end
end

fclose(evalFileWr);

% % -- open the previous optimized file --
% disp('Optimization Step 2 ...');
% %inputEvalFileName = 'output_gpu_evals_HxHt.txt';
% inputEvalFileName = 'output_gpu_evals_HxH.txt';
% 
% %outputEvalFileName = 'output_gpu_evals_HxHt_v2.txt';
% outputEvalFileName = 'output_gpu_evals_HxH_v2.txt';
% 
% fullInputFileName = fullfile(fileFolder, problem, inputEvalFileName);
% fullOutputFileName = fullfile(fileFolder, problem, outputEvalFileName);
% 
% evalFileRd_v2 = fopen(fullInputFileName, 'r');
% evalFileWr_v2 = fopen(fullOutputFileName, 'w');
% 
% ldata = textscan(evalFileRd_v2, '%s', 'Delimiter', ';', 'CollectOutput', true);
% line = string(ldata{1});
% 
% for i = 1:size(line, 1)
%     splitLine = split(line(i));
%     substituteStr_4 = splitLine(4);
%     substituteStr_6 = splitLine(6);
% 
%     % -- first write the lhs of the equation --
%     catStr_lhs = strcat(splitLine(2), {' '}, splitLine(3), {' '});
%     fprintf(evalFileWr_v2, catStr_lhs);
% 
%     % -- find the matched G's variable and substitue it by its aligned array representation --
%     for a = 1:size(valid_Gs, 2)
%         if strcmp(substituteStr_4, valid_Gs(a))
%             substituteStr_4 = array_Gs(a);
%         elseif strcmp(substituteStr_6, valid_Gs(a))
%             substituteStr_6 = array_Gs(a);
%         end
%     end
% 
%     catStr_rhs = strcat(substituteStr_4, {' '}, splitLine(5), {' '}, substituteStr_6, {';\n'});
%     fprintf(evalFileWr_v2, catStr_rhs);
% end
% 
% fclose(evalFileWr_v2);

% -- reformate the outputs --
% disp('Reformating the outputs ...');
% inputEvalFileName = 'M2-HxHt-Ccode-raw-output-y';
% %inputEvalFileName = 'M2-HxH-Ccode-raw-output-y';
% outputEvalFileName = 'output_gpu_evals_HxHt_y.txt';
% 
% fullInputFileName = fullfile(fileFolder, problem, inputEvalFileName);
% fullOutputFileName = fullfile(fileFolder, problem, outputEvalFileName);
% 
% evalFileRd = fopen(fullInputFileName, 'r');
% evalFileWr_output = fopen(fullOutputFileName, 'w');
% 
% ldata = textscan(evalFileRd, '%s', 'Delimiter', ';', 'CollectOutput', true);
% line = string(ldata{1});
% 
% outA_Str = '';
% startCoefStr = '';
% targetCoefStr = '';
% xStr = '';
% % -- strings of r_track registers --
% for i = 1:numOfVar
%     catStr_track = strcat({'r_cgesvA['}, string(i-1), {']'});
%     if i > 1
%         outA_Str = [outA_Str catStr_track];
%     else
%         outA_Str = catStr_track;
%     end
% end
% 
% fprintf(evalFileWr_output, 'switch(tx) {\n');
% for k = 0:numOfVar-1
%     case_catStr = strcat('case', {' '}, string(k), ':\n');
%     fprintf(evalFileWr_output, case_catStr);
%     write_A_indx = 1;
%     for i = k : numOfVar : numOfVar*numOfVar-1
%         splitLine = split(line(i+1));
%         substituteStr_3 = splitLine(3);
%         
%         % -- find the matched G's variable and substitue it by its aligned array representation --
%         for a = 1:size(valid_Gs, 2)
%             if strcmp(substituteStr_3, valid_Gs(a))
%                 substituteStr_3 = array_Gs(a);
%                 break;
%             end
%         end
%         
%         catStr = strcat(outA_Str(write_A_indx), {' '}, splitLine(2), {' '}, substituteStr_3, {';\n'});
%         fprintf(evalFileWr_output, catStr);
%         write_A_indx = write_A_indx + 1;
%     end
%     fprintf(evalFileWr_output, 'break;\n');
%     fprintf(evalFileWr_output, '\n');
% end
% 
% write_b_indx = 0;
% for i = numOfVar*numOfVar+1:size(line, 1)
%     splitLine = split(line(i));
%     substituteStr_3 = splitLine(3);
%         
%     % -- find the matched G's variable and substitue it by its aligned array representation --
%     for a = 1:size(valid_Gs, 2)
%         if strcmp(substituteStr_3, valid_Gs(a))
%             substituteStr_3 = array_Gs(a);
%             break;
%         end
%     end
%     
%     catStr = strcat({'r_cgesvB '}, splitLine(2), {' '}, substituteStr_3, {';\n'});
%     fprintf(evalFileWr_output, catStr);
%     write_b_indx = write_b_indx + 1;
% end
% 
% fclose(evalFileWr_output);

%% -- stored --
% for i = 1:size(line, 1)
%     decide_cont = 0;
%     splitLine = split(line(i));
%     %fprintf(evalFileWr, 'magmaFloatComplex ');
%     substituteStr_2 = splitLine(2);
%     for g = 1:numOfGvar
%         if strcmp(splitLine(2), G_Str(g))
%             if g == 1
%                 substituteStr_2 = G_Str_shmem(g);
%             else
%                 substituteStr_2 = G_Str_shmem(g-1);
%             end
%         end
%     end
%     %catStr = strcat(splitLine(2), {' '}, splitLine(3), {' '});
%     catStr = strcat(substituteStr_2, {' '}, splitLine(3), {' '});
%     fprintf(evalFileWr, catStr);
%     
%     if size(splitLine, 1) == 6
%         substituteStr_4 = splitLine(4);
%         substituteStr_6 = splitLine(6);
%         for g = 1:numOfGvar
%             if strcmp(splitLine(4), G_Str(g))
%                 if g == 1
%                     substituteStr_4 = G_Str_shmem(g);
%                 else
%                     substituteStr_4 = G_Str_shmem(g-1);
%                 end
%                decide_cont = decide_cont + 1;
%             elseif (strcmp(splitLine(6), G_Str(g)))
%                 if g == 1
%                     substituteStr_6 = G_Str_shmem(g);
%                 else
%                     substituteStr_6 = G_Str_shmem(g-1);
%                 end
%                decide_cont = decide_cont + 1;
%             end
%         end
%         
%         if decide_cont < 2
%             for j = 1:(numOfVar+1+numOfCoef+numOfCoef)
%                 if j <= numOfVar+1
%                     if strcmp(splitLine(6), xStr(j))
%                         substituteStr_6 = trackStr(j);
%                     end
%                     if strcmp(splitLine(4), xStr(j))
%                         substituteStr_4 = trackStr(j);
%                     end
%                 elseif j > numOfVar+1 && j <= (numOfVar + 1 + numOfCoef)
%                     p = j - (numOfVar+1);
%                     if strcmp(splitLine(6), xStr(j))
%                         substituteStr_6 = startCoefStr(p);
%                     end
%                     if strcmp(splitLine(4), xStr(j))
%                         substituteStr_4 = startCoefStr(p);
%                     end
%                 elseif j > (numOfVar + 1 + numOfCoef)
%                     q = j - (numOfVar + 1 + numOfCoef);
%                     if strcmp(splitLine(6), xStr(j))
%                         substituteStr_6 = targetCoefStr(q);
%                     end
%                     if strcmp(splitLine(4), xStr(j))
%                         substituteStr_4 = targetCoefStr(q);
%                     end
%                 end
% 
%             end
%         end
%         catStr_rhs = strcat(substituteStr_4, splitLine(5), substituteStr_6, {';\n'});
%         fprintf(evalFileWr, catStr_rhs);
%         
%     elseif size(splitLine, 1) == 8
%         catStr_rhs = strcat(splitLine(4), splitLine(5), splitLine(6), splitLine(7), splitLine(8), {';\n'});
%         fprintf(evalFileWr, catStr_rhs);
%         fprintf("G%d\n", i-1);
%     else
%         strShow = sprintf('Term %d has more than two arithmetic computations.\n', i-1);
%         disp(strShow);
%     end
% end
