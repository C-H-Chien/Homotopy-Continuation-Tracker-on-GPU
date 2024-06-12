function create_GPUHC_kernel_CU_code(problemName, max_order_of_func_t, ...
                                  max_num_of_terms_per_Hx_poly, max_num_of_parts_per_Hx_term, ...
                                  max_num_of_terms_per_poly, max_num_of_parts_per_term, ...
                                  num_of_vars, num_of_coeffs_from_params, outputFileWr_kernel_code)


    fid = fopen('GPUHC_kernel_template.txt', 'rt') ;
    kernel_template = fread(fid);
    fclose(fid);

    kernel_template = char(kernel_template.');
    
    %> Replace problem-specific information
    kernel_template = strrep(kernel_template, "<PROBLEM_NAME>",                 problemName);
    kernel_template = strrep(kernel_template, "<NUM_OF_VARS>",                  string(num_of_vars));
    kernel_template = strrep(kernel_template, "<NUM_OF_COEFFS_FROM_PARAMS>",    string(num_of_coeffs_from_params));
    kernel_template = strrep(kernel_template, "<NUM_OF_DHDX_MAX_TERMS>",        string(max_num_of_terms_per_Hx_poly));
    kernel_template = strrep(kernel_template, "<NUM_OF_DHDX_MAX_PARTS>",        string(max_num_of_parts_per_Hx_term));
    kernel_template = strrep(kernel_template, "<NUM_OF_DHDT_MAX_TERMS>",        string(max_num_of_terms_per_poly));
    kernel_template = strrep(kernel_template, "<NUM_OF_DHDT_MAX_PARTS>",        string(max_num_of_parts_per_term));
    kernel_template = strrep(kernel_template, "<MAX_ORDER_OF_T>",               string(max_order_of_func_t));

    %> Write to the file
    fwrite(outputFileWr_kernel_code, kernel_template);
end