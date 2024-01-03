//> ======================================================
//> Macros definitions
//> ======================================================

//> Repository directory
#define REPO_PATH                               std::string("/oscar/data/bkimia/cchien3/Homotopy-Continuation-Tracker-on-GPU/GPU-HC/")
#define WRITE_FILES_FOLDER                      std::string("Output_Write_Files/")
#define WRITE_FILES_PATH                        REPO_PATH + WRITE_FILES_FOLDER

//> A list of minimal problems (only one of them is true)
#define TRIFOCAL_2OP1P_30X30                    (false)
#define REL_POSE_5PT_ALG_FORM_QUAT              (false)
#define REL_POSE_5PT_GEO_FORM_QUAT              (true)

//> Homotopy Continuation Hyper-Parameters
#define HC_MAX_STEPS                            (42)
#define HC_MAX_CORRECTION_STEPS                 (4)
#define HC_NUM_OF_STEPS_TO_INCREASE_DELTA_T     (5)
#define APPLY_GAMMA_TRICK                       (false)
#define USE_DOUBLE_PRECISION                    (false)
#define GPU_DEBUG                               (false)

//> Define a random complex numbner gamma used in the gamma-trick
#if USE_DOUBLE_PRECISION
    #define GAMMA                               MAGMA_Z_MAKE(1 / std::sqrt(2), 1 / std::sqrt(2));
#else
    #define GAMMA                               MAGMA_C_MAKE(1 / std::sqrt(2), 1 / std::sqrt(2));
#endif

//> Problem specifications 
#if TRIFOCAL_2OP1P_30X30
    #define HC_PROBLEM                          std::string("trifocal_2op1p_30x30")
    #define NUM_OF_PARAMS                       (33)
    #define NUM_OF_TRACKS                       (312)
    #define NUM_OF_VARS                         (30)
    #define NUM_OF_COEFFS_FROM_PARAMS           (92)          //> don't care for now
    #define HX_MAXIMAL_TERMS                    (8)
    #define HX_MAXIMAL_PARTS                    (5)
    #define HT_MAXIMAL_TERMS                    (16)
    #define HT_MAXIMAL_PARTS                    (6)
    #define MAX_ORDER_OF_T                      (2)
    #define UNDEFINE_HC_PROBLEM                 (false)
#elif REL_POSE_5PT_ALG_FORM_QUAT
    #define HC_PROBLEM                          std::string("5pt_rel_pos_alg_form_quat")
    #define NUM_OF_PARAMS                       (20)
    #define NUM_OF_TRACKS                       (40)
    #define NUM_OF_VARS                         (6)
    #define NUM_OF_COEFFS_FROM_PARAMS           (75)
    #define HX_MAXIMAL_TERMS                    (12)
    #define HX_MAXIMAL_PARTS                    (4)
    #define HT_MAXIMAL_TERMS                    (30)
    #define HT_MAXIMAL_PARTS                    (5)
    #define MAX_ORDER_OF_T                      (2)
    #define UNDEFINE_HC_PROBLEM                 (false)
#elif REL_POSE_5PT_GEO_FORM_QUAT
    #define HC_PROBLEM                          std::string("5pt_rel_pos_geo_form_quat")
    #define NUM_OF_PARAMS                       (20)
    #define NUM_OF_TRACKS                       (40)
    #define NUM_OF_VARS                         (6)
    #define NUM_OF_COEFFS_FROM_PARAMS           (75)
    #define HX_MAXIMAL_TERMS                    (12)
    #define HX_MAXIMAL_PARTS                    (4)
    #define HT_MAXIMAL_TERMS                    (30)
    #define HT_MAXIMAL_PARTS                    (5)
    #define MAX_ORDER_OF_T                      (2)
    #define UNDEFINE_HC_PROBLEM                 (false)
#else
    #define UNDEFINE_HC_PROBLEM                 (true)
#endif



//> [DO NOT CHANGE] The following macros are constant. They are used for shuffle operation in a warp level. 
#define FULL_MASK                               (0xffffffff)
#define WARP_SIZE                               (32)

//> CUDA error check
#define cudacheck( a )  do { \
                            cudaError_t e = a; \
                            if(e != cudaSuccess) { \
                                printf("\033[1;31m"); \
                                printf("Error in %s:%d %s\n", __func__, __LINE__, cudaGetErrorString(e)); \
                                printf("\033[0m"); \
                            }\
                        } while(0)


