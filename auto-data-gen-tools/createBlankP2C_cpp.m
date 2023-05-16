function createBlankP2C_cpp(problemName, outputFileWr_P2C_script)
    
    fprintf(outputFileWr_P2C_script, "#ifndef ");
    fprintf(outputFileWr_P2C_script, strcat(problemName, "_h\n"));
    fprintf(outputFileWr_P2C_script, "#define ");
    fprintf(outputFileWr_P2C_script, strcat(problemName, "_h\n\n"));
    
    fprintf(outputFileWr_P2C_script, "#include <stdio.h>\n");
    fprintf(outputFileWr_P2C_script, "#include <stdlib.h>\n");
    fprintf(outputFileWr_P2C_script, "#include <cstdio>\n");
    fprintf(outputFileWr_P2C_script, "#include <iostream>\n");
    fprintf(outputFileWr_P2C_script, "#include <iomanip>\n");
    fprintf(outputFileWr_P2C_script, "#include <cstring>\n");
    fprintf(outputFileWr_P2C_script, "#include <chrono>\n\n");
    
    fprintf(outputFileWr_P2C_script, '#include "flops.h"\n');
    fprintf(outputFileWr_P2C_script, '#include "magma_v2.h"\n');
    fprintf(outputFileWr_P2C_script, '#include "magma_lapack.h"\n');
    fprintf(outputFileWr_P2C_script, '#include "magma_internal.h"\n');
    fprintf(outputFileWr_P2C_script, '#undef max\n');
    fprintf(outputFileWr_P2C_script, '#undef min\n\n');
    
    fprintf(outputFileWr_P2C_script, "namespace magmaHCWrapper{\n\n");
    fprintf(outputFileWr_P2C_script, "void ");
    fprintf(outputFileWr_P2C_script, strcat(problemName, "("));
    fprintf(outputFileWr_P2C_script, "magmaFloatComplex *h_targetParams, magmaFloatComplex *h_startParams, magmaFloatComplex *h_phc_coeffs_Hx, magmaFloatComplex *h_phc_coeffs_Ht)\n");
    fprintf(outputFileWr_P2C_script, "{\n");

end