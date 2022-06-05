function constant_mat_write2file_improve(fileFolder, category, problem, wrFolder, const_matrix_Hx, const_matrix_Ht)
                             
    % -- define file names --
    outputFileName_Hx = 'Hx.txt';
    outputFileName_Ht = 'Ht.txt';

    fullOutputFileName_Hx = fullfile(fileFolder, category, problem, wrFolder, outputFileName_Hx);
    fullOutputFileName_Ht = fullfile(fileFolder, category, problem, wrFolder, outputFileName_Ht);

    outputFileWr_Hx = fopen(fullOutputFileName_Hx, 'w');
    outputFileWr_Ht = fopen(fullOutputFileName_Ht, 'w');
    
    for i = 1:size(const_matrix_Hx, 1)
        for j = 1:size(const_matrix_Hx, 2)
            fprintf(outputFileWr_Hx, num2str(const_matrix_Hx(i,j)));
            fprintf(outputFileWr_Hx, '\t');
        end
        
        for j = 1:size(const_matrix_Ht, 2)
            fprintf(outputFileWr_Ht, num2str(const_matrix_Ht(i,j)));
            fprintf(outputFileWr_Ht, '\t');
        end
        
        fprintf(outputFileWr_Hx, '\n');
        fprintf(outputFileWr_Ht, '\n');
    end
    
    fclose(outputFileWr_Hx);
    fclose(outputFileWr_Ht);
end