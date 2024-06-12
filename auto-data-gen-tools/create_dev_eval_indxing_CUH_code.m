function create_dev_eval_indxing_CUH_code(problemName, max_order_of_func_t, max_num_of_dHdX_parts, outputFileWr_dev_eval_indxing_script)
    
    fprintf(outputFileWr_dev_eval_indxing_script, "#ifndef DEV_EVAL_INDXING_");
    fprintf(outputFileWr_dev_eval_indxing_script, strcat(problemName, "_CUH\n"));
    fprintf(outputFileWr_dev_eval_indxing_script, "#define DEV_EVAL_INDXING_");
    fprintf(outputFileWr_dev_eval_indxing_script, strcat(problemName, "_CUH\n\n"));
    
    fprintf(outputFileWr_dev_eval_indxing_script, "#include <stdio.h>\n");
    fprintf(outputFileWr_dev_eval_indxing_script, "#include <stdlib.h>\n");
    fprintf(outputFileWr_dev_eval_indxing_script, "#include <cstdio>\n");

    fprintf(outputFileWr_dev_eval_indxing_script, '#include "magma_v2.h"\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '#include "magma_lapack.h"\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '#include "magma_internal.h"\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '#undef max\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '#undef min\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '#include "magma_templates.h"\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '#include "sync.cuh"\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '#undef max\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '#undef min\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '#include "shuffle.cuh"\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '#undef max\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '#undef min\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '#include "batched_kernel_param.h"\n\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '#include "../definitions.hpp"\n\n');

    %> Device functrion for evaluating the linear interpolations of
    %  parameters of PHC
    fprintf(outputFileWr_dev_eval_indxing_script, 'template< int Num_Of_Vars, int Max_Order_of_t, unsigned Full_Parallel_Offset, \n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t\t unsigned Partial_Parallel_Thread_Offset, unsigned Partial_Parallel_Index_Offset, \n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t\t unsigned Max_Order_of_t_Plus_One, unsigned Partial_Parallel_Index_Offset_Hx, unsigned Partial_Parallel_Index_Offset_Ht >\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '__device__ __inline__ void\n');
    fprintf(outputFileWr_dev_eval_indxing_script, 'eval_parameter_homotopy(\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t const int tx, float t, magmaFloatComplex *s_phc_coeffs_Hx, magmaFloatComplex *s_phc_coeffs_Ht,\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t const magmaFloatComplex __restrict__ *d_phc_coeffs_Hx, const magmaFloatComplex __restrict__ *d_phc_coeffs_Ht )\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '{\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t #pragma unroll 2\n');

    fprintf(outputFileWr_dev_eval_indxing_script, '\t for (int i = 0; i < Full_Parallel_Offset; i++) {\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t   s_phc_coeffs_Hx[ tx + i*Num_Of_Vars ] = d_phc_coeffs_Hx[ tx*Max_Order_of_t_Plus_One + i*Num_Of_Vars*Max_Order_of_t_Plus_One ]\n');
    for i = 1:max_order_of_func_t
        
        fprintf(outputFileWr_dev_eval_indxing_script, strcat('\t\t + d_phc_coeffs_Hx[ tx*Max_Order_of_t_Plus_One + ', string(i),' + i*Num_Of_Vars*Max_Order_of_t_Plus_One ]'));
        for j = 1:i
            fprintf(outputFileWr_dev_eval_indxing_script, ' * t');
        end
        if i < max_order_of_func_t
            fprintf(outputFileWr_dev_eval_indxing_script, '\n');
        else
            fprintf(outputFileWr_dev_eval_indxing_script, ';\n');
        end
    end

    fprintf(outputFileWr_dev_eval_indxing_script, '\t   s_phc_coeffs_Ht[ tx + i*Num_Of_Vars ] = d_phc_coeffs_Ht[ tx*Max_Order_of_t + i*Num_Of_Vars*Max_Order_of_t ]\n');
    for i = 1:max_order_of_func_t-1
        
        fprintf(outputFileWr_dev_eval_indxing_script, strcat('\t\t + d_phc_coeffs_Ht[ tx*Max_Order_of_t + ', string(i),' + i*Num_Of_Vars*Max_Order_of_t ]'));
        for j = 1:i
            fprintf(outputFileWr_dev_eval_indxing_script, ' * t');
        end
        if i < max_order_of_func_t-1
            fprintf(outputFileWr_dev_eval_indxing_script, '\n');
        else
            fprintf(outputFileWr_dev_eval_indxing_script, ';\n');
        end
    end
    fprintf(outputFileWr_dev_eval_indxing_script, '\t}\n\n\t if (tx < Partial_Parallel_Thread_Offset) {\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t   s_phc_coeffs_Hx[ tx + Partial_Parallel_Index_Offset ] = d_phc_coeffs_Hx[ tx*Max_Order_of_t_Plus_One + Partial_Parallel_Index_Offset_Hx ]\n');
    for i = 1:max_order_of_func_t
        fprintf(outputFileWr_dev_eval_indxing_script, strcat('\t\t + d_phc_coeffs_Hx[ tx*Max_Order_of_t_Plus_One + Partial_Parallel_Index_Offset_Hx + ', string(i), ' ] '));
        for j = 1:i
            fprintf(outputFileWr_dev_eval_indxing_script, ' * t');
        end
        if i < max_order_of_func_t
            fprintf(outputFileWr_dev_eval_indxing_script, '\n');
        else
            fprintf(outputFileWr_dev_eval_indxing_script, ';\n');
        end
    end
    fprintf(outputFileWr_dev_eval_indxing_script, '\t   s_phc_coeffs_Ht[ tx + Partial_Parallel_Index_Offset ] = d_phc_coeffs_Ht[ tx*Max_Order_of_t + Partial_Parallel_Index_Offset_Ht ]\n');
    for i = 1:max_order_of_func_t-1
        fprintf(outputFileWr_dev_eval_indxing_script, strcat('\t\t + d_phc_coeffs_Ht[ tx*Max_Order_of_t + Partial_Parallel_Index_Offset_Ht + ', string(i), ' ] '));
        for j = 1:i
            fprintf(outputFileWr_dev_eval_indxing_script, ' * t');
        end
        if i < max_order_of_func_t-1
            fprintf(outputFileWr_dev_eval_indxing_script, '\n');
        else
            fprintf(outputFileWr_dev_eval_indxing_script, ';\n');
        end
    end
    fprintf(outputFileWr_dev_eval_indxing_script, '\t}\n}\n\n');

    %> Device function for evaluating the Jacobian matrix dH/dx
    fprintf(outputFileWr_dev_eval_indxing_script, 'template< int Num_Of_Vars, int dHdx_Max_Terms, int dHdx_Max_Parts, int dHdx_Entry_Offset, int dHdx_Row_Offset >\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '__device__ __inline__ void\n');
    fprintf(outputFileWr_dev_eval_indxing_script, 'eval_Jacobian_Hx(\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t const int tx, magmaFloatComplex *s_track, magmaFloatComplex r_cgesvA[Num_Of_Vars],\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t const int* __restrict__ d_Hx_idx, magmaFloatComplex *s_phc_coeffs )\n{\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t for(int i = 0; i < Num_Of_Vars; i++) {\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t\t r_cgesvA[i] = MAGMA_C_ZERO;\n\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t\t #pragma unroll 2\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t\t for(int j = 0; j < dHdx_Max_Terms; j++) {\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t\t\t r_cgesvA[i] += d_Hx_idx[j*dHdx_Max_Parts + i*dHdx_Entry_Offset + tx*dHdx_Row_Offset] \n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t\t\t              * s_phc_coeffs[ d_Hx_idx[j*dHdx_Max_Parts + 1 + i*dHdx_Entry_Offset + tx*dHdx_Row_Offset] ]\n');
    for i = 2:max_num_of_dHdX_parts-1
        fprintf(outputFileWr_dev_eval_indxing_script, strcat('\t\t\t\t * s_track[      d_Hx_idx[j*dHdx_Max_Parts + ', string(i),' + i*dHdx_Entry_Offset + tx*dHdx_Row_Offset] ]'));
        if i < max_num_of_dHdX_parts-1
            fprintf(outputFileWr_dev_eval_indxing_script, '\n');
        else
            fprintf(outputFileWr_dev_eval_indxing_script, ';\n');
        end
    end
    fprintf(outputFileWr_dev_eval_indxing_script, '\t\t }\n\t }\n}\n\n');

    %> Device function for evaluating the Jacobian matrix dH/dt
    fprintf(outputFileWr_dev_eval_indxing_script, 'template< int dHdt_Max_Terms, int dHdt_Max_Parts, int dHdt_Row_Offset >\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '__device__ __inline__ void\n');
    fprintf(outputFileWr_dev_eval_indxing_script, 'eval_Jacobian_Ht(\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t const int tx, magmaFloatComplex *s_track, magmaFloatComplex &r_cgesvB,\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t const int* __restrict__ d_Ht_idx, magmaFloatComplex *s_phc_coeffs )\n{\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t\t r_cgesvB = MAGMA_C_ZERO;\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t\t #pragma unroll 2\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t\t for (int i = 0; i < dHdt_Max_Terms; i++) {\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t\t\t r_cgesvB -= d_Ht_idx[i*dHdt_Max_Parts + tx*dHdt_Row_Offset]\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t\t\t              * s_phc_coeffs[ d_Ht_idx[i*dHdt_Max_Parts + 1 + tx*dHdt_Row_Offset] ]\n');
    for i = 2:max_num_of_dHdX_parts
        fprintf(outputFileWr_dev_eval_indxing_script, strcat('\t\t\t\t * s_track[      d_Ht_idx[i*dHdt_Max_Parts + ', string(i), ' + tx*dHdt_Row_Offset] ]'));
        if i < max_num_of_dHdX_parts
            fprintf(outputFileWr_dev_eval_indxing_script, '\n');
        else
            fprintf(outputFileWr_dev_eval_indxing_script, ';\n\t\t }\n}\n\n');
        end
    end
    
    %> Device function for evaluating the Homotopy H
    fprintf(outputFileWr_dev_eval_indxing_script, 'template< int dHdt_Max_Terms, int dHdt_Max_Parts, int dHdt_Row_Offset >\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '__device__ __inline__ void\neval_Homotopy(\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t const int tx, magmaFloatComplex *s_track, magmaFloatComplex &r_cgesvB,\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t const int* __restrict__ d_Ht_idx, magmaFloatComplex *s_phc_coeffs)\n{\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t\t r_cgesvB = MAGMA_C_ZERO;\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t\t #pragma unroll 2\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t\t for (int i = 0; i < dHdt_Max_Terms; i++) {\n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t\t\t r_cgesvB += d_Ht_idx[i*dHdt_Max_Parts + tx*dHdt_Row_Offset] \n');
    fprintf(outputFileWr_dev_eval_indxing_script, '\t\t\t * s_phc_coeffs[ d_Ht_idx[i*dHdt_Max_Parts + 1 + tx*dHdt_Row_Offset] ] \n');
    for i = 2:max_num_of_dHdX_parts
        fprintf(outputFileWr_dev_eval_indxing_script, strcat('\t\t\t * s_track[      d_Ht_idx[i*dHdt_Max_Parts + ', string(i), ' + tx*dHdt_Row_Offset] ]'));
        if i < max_num_of_dHdX_parts
            fprintf(outputFileWr_dev_eval_indxing_script, '\n');
        else
            fprintf(outputFileWr_dev_eval_indxing_script, ';\n\t\t }\n}\n\n');
        end
    end

    fprintf(outputFileWr_dev_eval_indxing_script, '#endif');
end