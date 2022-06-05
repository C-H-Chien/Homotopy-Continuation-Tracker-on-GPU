function constant_mat_write2file(fileFolder, problem, category, numOfVars, Ht_numOfX, const_matrix_Hx_scalar, const_matrix_Hx_x, const_matrix_Hx_Y, ...
                                 const_matrix_Ht_scalar, const_matrix_Ht_x, const_matrix_Ht_Y)
    
    % -- maximal parts of X in Hx --
    Hx_X_maximal_parts = size(const_matrix_Hx_x, 3);
                             
    % -- define file names --
    outputFileName_Hx_scalar = 'Hx_scalar.txt';
    outputFileName_Hx_X = 'Hx_X.txt';
    fullOutputFileName_Hx_X = fullfile(fileFolder, problem, category, outputFileName_Hx_X);
    outputFileWr_Hx_X = fopen(fullOutputFileName_Hx_X, 'w');
    outputFileName_Hx_Y = 'Hx_Y.txt';
    outputFileName_Ht_scalar = 'Ht_scalar.txt';
    outputFileName_Ht_X = 'Ht_X.txt';
    outputFileName_Ht_Y = 'Ht_Y.txt';

    fullOutputFileName_Hx_scalar = fullfile(fileFolder, problem, category, outputFileName_Hx_scalar);
    fullOutputFileName_Hx_Y = fullfile(fileFolder, problem, category, outputFileName_Hx_Y);
    fullOutputFileName_Ht_scalar = fullfile(fileFolder, problem, category, outputFileName_Ht_scalar);
    fullOutputFileName_Ht_X = fullfile(fileFolder, problem, category, outputFileName_Ht_X);
    fullOutputFileName_Ht_Y = fullfile(fileFolder, problem, category, outputFileName_Ht_Y);

    outputFileWr_Hx_scalar = fopen(fullOutputFileName_Hx_scalar, 'w');    
    outputFileWr_Hx_Y = fopen(fullOutputFileName_Hx_Y, 'w');
    outputFileWr_Ht_scalar = fopen(fullOutputFileName_Ht_scalar, 'w');
    outputFileWr_Ht_X = fopen(fullOutputFileName_Ht_X, 'w');
    outputFileWr_Ht_Y = fopen(fullOutputFileName_Ht_Y, 'w');
    
    % -- write the constant matrices to the output files --
    for xf = 1:Hx_X_maximal_parts
        for i = 1:size(const_matrix_Hx_scalar, 2)
            for j = 1:size(const_matrix_Hx_scalar, 1)
                if xf == 1
                    % -- Hx scalar matrix --
                    wr_str = strcat('s[', num2str((i-1)*numOfVars*numOfVars + j-1), '] = ', num2str(const_matrix_Hx_scalar(j,i)), ';');
                    fprintf(outputFileWr_Hx_scalar, wr_str);
                    fprintf(outputFileWr_Hx_scalar, "\n");
                    
                    % -- Hx Y --
                    wr_str = strcat('Y[', num2str((i-1)*numOfVars*numOfVars + j-1), '] = ', num2str(const_matrix_Hx_Y(j,i)), ';');
                    fprintf(outputFileWr_Hx_Y, wr_str);
                    fprintf(outputFileWr_Hx_Y, "\n");

                    if mod(j, numOfVars) == 0
                        fprintf(outputFileWr_Hx_scalar, "\n");
                        fprintf(outputFileWr_Hx_Y, "\n");
                    end
                end
                
                % -- Hx X --
                wr_str = strcat('X[', num2str((i-1)*numOfVars*numOfVars + j-1 + size(const_matrix_Hx_scalar, 2)*numOfVars*numOfVars*(xf-1)), '] = ', num2str(const_matrix_Hx_x(j,i,xf)), ';');
                fprintf(outputFileWr_Hx_X, wr_str);
                fprintf(outputFileWr_Hx_X, "\n");
            end
        end
        fprintf(outputFileWr_Hx_X, "\n");
    end

    for numOfX = 1:Ht_numOfX
        for i = 1:size(const_matrix_Ht_scalar, 2)
            for j = 1:size(const_matrix_Ht_scalar, 1)
                if numOfX == 1
                    % -- Ht scalar matrix --
                    wr_str = strcat('s[', num2str((i-1)*numOfVars + j-1), '] = ', num2str(const_matrix_Ht_scalar(j,i)), ';');
                    fprintf(outputFileWr_Ht_scalar, wr_str);
                    fprintf(outputFileWr_Ht_scalar, "\n");

                    % -- Ht Y --
                    wr_str = strcat('Y[', num2str((i-1)*numOfVars + j-1), '] = ', num2str(const_matrix_Ht_Y(j,i)), ';');
                    fprintf(outputFileWr_Ht_Y, wr_str);
                    fprintf(outputFileWr_Ht_Y, "\n");

%                     % -- H Y --
%                     wr_str = strcat('Y[', num2str((i-1)*numOfVars + j-1), '] = ', num2str(const_matrix_H_Y(j,i)), ';');
%                     fprintf(outputFileWr_H_Y, wr_str);
%                     fprintf(outputFileWr_H_Y, "\n");
                end

                % -- Ht X collection --
                wr_str = strcat('Xt[', num2str((i-1)*numOfVars + j-1 + (numOfX-1)*(size(const_matrix_Ht_scalar, 2)*numOfVars)), '] = ', num2str(const_matrix_Ht_x(j,i,numOfX)), ';');
                fprintf(outputFileWr_Ht_X, wr_str);
                fprintf(outputFileWr_Ht_X, "\n");

%                 % -- H X collection --
%                 wr_str = strcat('H_X[', num2str((i-1)*numOfVars + j-1 + (numOfX-1)*(size(const_matrix_Ht_scalar, 2)*numOfVars)), '] = ', num2str(const_matrix_H_x(j,i,numOfX)), ';');
%                 fprintf(outputFileWr_H_X, wr_str);
%                 fprintf(outputFileWr_H_X, "\n");
            end

            if numOfX == 1
                fprintf(outputFileWr_Ht_scalar, "\n");
                fprintf(outputFileWr_Ht_Y, "\n");
                %fprintf(outputFileWr_H_Y, "\n");
            end
            fprintf(outputFileWr_Ht_X, "\n");
            %fprintf(outputFileWr_H_X, "\n");
        end
        fprintf(outputFileWr_Ht_X, "\n");
        %fprintf(outputFileWr_H_X, "\n");
    end
    
    fclose(outputFileWr_Hx_scalar);
    fclose(outputFileWr_Hx_X);
    fclose(outputFileWr_Hx_Y);
    fclose(outputFileWr_Ht_scalar);
    fclose(outputFileWr_Ht_X);
    fclose(outputFileWr_Ht_Y);
%     fclose(outputFileWr_H_X);
%     fclose(outputFileWr_H_Y);
    
end