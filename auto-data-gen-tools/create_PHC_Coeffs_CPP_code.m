function max_order_of_func_t = create_PHC_Coeffs_CPP_code(problemName, rep_coeffs_stack, numOfParams, outputFileWr_PHC_script)
    
    fprintf(outputFileWr_PHC_script, "#ifndef P2C_");
    fprintf(outputFileWr_PHC_script, strcat(problemName, "_H\n"));
    fprintf(outputFileWr_PHC_script, "#define P2C_");
    fprintf(outputFileWr_PHC_script, strcat(problemName, "_H\n\n"));
    
    fprintf(outputFileWr_PHC_script, "#include <stdio.h>\n");
    fprintf(outputFileWr_PHC_script, "#include <stdlib.h>\n");
    fprintf(outputFileWr_PHC_script, "#include <cstdio>\n");

    fprintf(outputFileWr_PHC_script, '#include "magma_v2.h"\n');
    fprintf(outputFileWr_PHC_script, '#include "magma_lapack.h"\n');
    fprintf(outputFileWr_PHC_script, '#include "magma_internal.h"\n');
    fprintf(outputFileWr_PHC_script, '#undef max\n');
    fprintf(outputFileWr_PHC_script, '#undef min\n\n');
    
    fprintf(outputFileWr_PHC_script, "void p2c-");
    fprintf(outputFileWr_PHC_script, strcat(problemName, "("));
    fprintf(outputFileWr_PHC_script, "magmaFloatComplex *h_targetParams, magmaFloatComplex *h_startParams, magmaFloatComplex *h_phc_coeffs_Hx, magmaFloatComplex *h_phc_coeffs_Ht)\n");
    fprintf(outputFileWr_PHC_script, "{\n");

    max_order_of_func_t = createUnivPoly4Coeffs(rep_coeffs_stack, numOfParams, outputFileWr_PHC_script);

    fprintf(outputFileWr_PHC_script, "}\n");
    fprintf(outputFileWr_PHC_script, "#endif\n");
end