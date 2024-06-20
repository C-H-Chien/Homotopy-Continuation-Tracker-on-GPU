//> ======================================================
//> Macros definitions
//> ======================================================

#define WRITE_FILES_FOLDER                      std::string("Output_Write_Files/")

//> RANSAC Settings
#define NUM_OF_RANSAC_ITERATIONS                (1)
#define IMAG_PART_TOL                           (1e-1)  //(1e-5)
#define ROT_RESIDUAL_TOL                        (1e-1)
#define TRANSL_RESIDUAL_TOL                     (1e-1)
#define TEST_RANSAC_TIMES                       (100)
#define REPROJ_ERROR_INLIER_THRESH              (2) //> in pixels
#define FEED_RANDOM_SEED                        (true)     //> if false, each run is different

//> Settings for GPU-HC Kernel
#define TRUNCATE_HC_PATH_THRESH                 (0.95)

//> Evaluation macros
#define DUPLICATE_SOL_DIFF_TOL                  (1e-4)
#define ZERO_IMAG_PART_TOL_FOR_SP               (1e-1)  //(1e-4)
#define DEBUG_EVALUATOR                         (false)
#define IS_SO3_DET_R_TOL                        (1e-5)

//> Settings for Debugging
#define SHOW_PROBLEM_SETTINGS                   (true)
#define SHOW_EVAL_INDX_DATA_SIZE                (false)
#define GPU_DEBUG                               (true)
#define DATA_READER_DEBUG                       (false)
#define RANSAC_DEBUG                            (false)

//> [DO NOT CHANGE] The following macros are constant. They are used for shuffle operation in a warp level.
#define FULL_MASK                               (0xffffffff)
#define WARP_SIZE                               (32)
#define HALF_WARP_SIZE                          (16)

//> [DO NOT CHANGE] Constant values
#define ONE_OVER_SIX                            (0.166666666667)

//> CUDA error check
#define cudacheck( a )  do { \
                            cudaError_t e = a; \
                            if(e != cudaSuccess) { \
                                printf("\033[1;31m"); \
                                printf("Error in %s:%d %s\n", __func__, __LINE__, cudaGetErrorString(e)); \
                                printf("\033[0m"); \
                            }\
                        } while(0)

#define LOG_INFO_MESG(info_msg)        printf("\033[1;32m[INFO] %s\033[0m\n", std::string(info_msg).c_str() );
#define LOG_FILE_ERROR(err_msg)         printf("\033[1;31m[ERROR] File %s not found!\033[0m\n", std::string(err_msg).c_str() );
#define LOG_ERROR(err_msg)              printf("\033[1;31m[ERROR] %s\033[0m\n", std::string(err_msg).c_str() );
#define LOG_DATA_LOAD_ERROR(err_msg)    printf("\033[1;31m[DATA LOAD ERROR] %s not loaded successfully!\033[0m\n", std::string(err_msg).c_str() );
#define LOG_PRINT_HELP_MESSAGE          printf("Usage: ./magmaHC-main [options] [path]\n\n" \
                                               "options:\n" \
                                               "  -h, --help        show this help message and exit\n" \
                                               "  -d, --directory   repository directory, e.g. /home/chchien/Homotopy-Continuation-Tracker-on-GPU/GPU-HC/\n");


