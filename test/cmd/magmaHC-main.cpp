#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>

// magma
#include "flops.h"
#include "magma_v2.h"
#include "magma_lapack.h"

int main() {
  magma_init();
  magma_print_environment();

  real_Double_t   gflops, cpu_perf, cpu_time, gpu_perf, gpu_time;
  float          error, Rnorm, Anorm, Xnorm, *work;
  magmaFloatComplex c_one     = MAGMA_C_ONE;
  magmaFloatComplex c_neg_one = MAGMA_C_NEG_ONE;
  magmaFloatComplex *h_A, *h_B, *h_X;
  magmaFloatComplex_ptr d_A, d_B;
  magma_int_t *dipiv, *dinfo_array;
  magma_int_t *ipiv, *cpu_info;
  magma_int_t N, nrhs, lda, ldb, ldda, lddb, info, sizeA, sizeB;
  magma_int_t ione = 1;
  magma_int_t ISEED[4] = {0,0,0,1};
  int status = 0;
  magma_int_t batchCount = 2;
  nrhs = 1;

  magmaFloatComplex **dA_array = NULL;
  magmaFloatComplex **dB_array = NULL;
  magma_int_t     **dipiv_array = NULL;

  bool use_lapack = 1;
  double tol = 0.000001;
  N = 6;

  magma_queue_t my_queue;    // magma queue variable, internally holds a cuda stream and a cublas handle
  magma_device_t cdev;       // variable to indicate current gpu id

  magma_getdevice( &cdev );
  magma_queue_create( cdev, &my_queue );     // create a queue on this cdev

  printf("%% BatchCount   N  NRHS   CPU Gflop/s (msec)   GPU Gflop/s (msec)   ||B - AX|| / N*||A||*||X||\n");
  printf("%%============================================================================================\n");
  lda    = N;
  ldb    = lda;
  ldda   = magma_roundup( N, 32 );  // multiple of 32 by default
  lddb   = ldda;
  gflops = ( FLOPS_DGETRF( N, N ) + FLOPS_DGETRS( N, nrhs ) ) * batchCount / 1e9;

  sizeA = lda*N*batchCount;
  sizeB = ldb*nrhs*batchCount;

  magma_cmalloc_cpu( &h_A, sizeA );
  magma_cmalloc_cpu( &h_B, sizeB );
  magma_cmalloc_cpu( &h_X, sizeB );
  magma_smalloc_cpu( &work, N );
  magma_imalloc_cpu( &ipiv, batchCount*N );
  magma_imalloc_cpu( &cpu_info, batchCount );

  magma_cmalloc( &d_A, ldda*N*batchCount    );
  magma_cmalloc( &d_B, lddb*nrhs*batchCount );
  magma_imalloc( &dipiv, N * batchCount );
  magma_imalloc( &dinfo_array, batchCount );

  magma_malloc( (void**) &dA_array,    batchCount * sizeof(magmaFloatComplex*) );
  magma_malloc( (void**) &dB_array,    batchCount * sizeof(magmaFloatComplex*) );
  magma_malloc( (void**) &dipiv_array, batchCount * sizeof(magma_int_t*) );

  /* Initialize the matrices */
  lapackf77_clarnv( &ione, ISEED, &sizeA, h_A );
  lapackf77_clarnv( &ione, ISEED, &sizeB, h_B );

  magma_csetmatrix( N, N*batchCount,    h_A, lda, d_A, ldda, my_queue );
  magma_csetmatrix( N, nrhs*batchCount, h_B, ldb, d_B, lddb, my_queue );

  /* ====================================================================
     Performs operation using MAGMA
     =================================================================== */
  magma_cset_pointer( dA_array, d_A, ldda, 0, 0, ldda*N, batchCount, my_queue );
  magma_cset_pointer( dB_array, d_B, lddb, 0, 0, lddb*nrhs, batchCount, my_queue );
  magma_iset_pointer( dipiv_array, dipiv, 1, 0, 0, N, batchCount, my_queue );

  gpu_time = magma_sync_wtime( my_queue );
  info = magma_cgesv_batched(N, nrhs, dA_array, ldda, dipiv_array, dB_array, lddb, dinfo_array, batchCount, my_queue);
  gpu_time = magma_sync_wtime( my_queue ) - gpu_time;
  gpu_perf = gflops / gpu_time;

  // check correctness of results throught "dinfo_magma" and correctness of argument throught "info"
  magma_getvector( batchCount, sizeof(magma_int_t), dinfo_array, 1, cpu_info, 1, my_queue );
  for (int i=0; i < batchCount; i++)
  {
      if (cpu_info[i] != 0 ) {
          printf("magma_dgesv_batched matrix %lld returned internal error %lld\n",
                (long long) i, (long long) cpu_info[i] );
      }
  }
  if (info != 0) {
      printf("magma_dgesv_batched returned argument error %lld: %s.\n",
            (long long) info, magma_strerror( info ));
  }
  std::cout << "Complete magma_cgesv_batched" << std::endl;

  //=====================================================================
  // Residual
  //=====================================================================
  magma_cgetmatrix( N, nrhs*batchCount, d_B, lddb, h_X, ldb, my_queue );

  error = 0;
  for (magma_int_t s=0; s < batchCount; s++)
  {
      Anorm = lapackf77_clange("I", &N, &N,    h_A + s * lda * N, &lda, work);
      Xnorm = lapackf77_clange("I", &N, &nrhs, h_X + s * ldb * nrhs, &ldb, work);

      blasf77_cgemm( MagmaNoTransStr, MagmaNoTransStr, &N, &nrhs, &N,
                 &c_one,     h_A + s * lda * N, &lda,
                             h_X + s * ldb * nrhs, &ldb,
                 &c_neg_one, h_B + s * ldb * nrhs, &ldb);

      Rnorm = lapackf77_clange("I", &N, &nrhs, h_B + s * ldb * nrhs, &ldb, work);
      float err = Rnorm/(N*Anorm*Xnorm);

      if (std::isnan(err) || std::isinf(err)) {
          error = err;
          break;
      }
      error = max( err, error );
  }
  bool okay = (error < tol);
  status += ! okay;
  std::cout << "Complete calculating residuals" << std::endl;

  /* ====================================================================
     Performs operation using LAPACK
     =================================================================== */
  if ( use_lapack ) {
      cpu_time = magma_wtime();
      // #define BATCHED_DISABLE_PARCPU
      #if !defined (BATCHED_DISABLE_PARCPU) && defined(_OPENMP)
      magma_int_t nthreads = magma_get_lapack_numthreads();
      magma_set_lapack_numthreads(1);
      magma_set_omp_numthreads(nthreads);
      #pragma omp parallel for schedule(dynamic)
      #endif
      for (magma_int_t s=0; s < batchCount; s++)
      {
          magma_int_t locinfo;
          lapackf77_cgesv( &N, &nrhs, h_A + s * lda * N, &lda, ipiv + s * N, h_B + s * ldb * nrhs, &ldb, &locinfo );
          if (locinfo != 0) {
              printf("lapackf77_cgesv matrix %lld returned error %lld: %s.\n",
                      (long long) s, (long long) locinfo, magma_strerror( locinfo ));
          }
      }
      #if !defined (BATCHED_DISABLE_PARCPU) && defined(_OPENMP)
          magma_set_lapack_numthreads(nthreads);
      #endif
      cpu_time = magma_wtime() - cpu_time;
      cpu_perf = gflops / cpu_time;
      printf( "%10lld %5lld %5lld   %7.2f (%7.2f)   %7.2f (%7.2f)   %8.2e   %s\n",
              (long long) batchCount, (long long) N, (long long) nrhs,
              cpu_perf, cpu_time*1000, gpu_perf, gpu_time*1000,
              error, (okay ? "ok" : "failed"));
  }
  else {
      printf( "%10lld %5lld %5lld     ---   (  ---  )   %7.2f (%7.2f)   %8.2e   %s\n",
              (long long) batchCount, (long long) N, (long long) nrhs,
              gpu_perf, gpu_time,
              error, (okay ? "ok" : "failed"));
  }

  magma_queue_destroy( my_queue );

  magma_free_cpu( h_A );
  magma_free_cpu( h_B );
  magma_free_cpu( h_X );
  magma_free_cpu( work );
  magma_free_cpu( ipiv );
  magma_free_cpu( cpu_info );

  magma_free( d_A );
  magma_free( d_B );

  magma_free( dipiv );
  magma_free( dinfo_array );

  magma_free( dA_array );
  magma_free( dB_array );
  magma_free( dipiv_array );

  fflush( stdout );

  printf( "\n" );

  magma_finalize();
}