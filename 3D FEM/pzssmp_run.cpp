#include <cstdio>
#include <cstdlib>
#include <complex>
#include "define.h"
// #include "struct.h"
// #include "iteration.h"
// #include <mkl_solver.h>
#include <mkl.h>
#include <sys/time.h>
#include <omp.h>

#ifdef PARDISO_FLOAT
void callPARDISO(long long int nn, long long int nc, 
		 long long int *ia, long long int *ja, 
		 std::complex<float> *avals, std::complex<float> *bb, 
		 std::complex<float> *phi, 
		 long long int *arguments, long long int nrhs) {
#else
void callPARDISO(long long int nn, long long int nc, 
		 long long int *ia, long long int *ja, 
		 std::complex<double> *avals, std::complex<double> *bb, 
		 std::complex<double> *phi, 
		 long long int *arguments, long long int nrhs) {
#endif
  long long int i, j;
  long long int count;
  char *var;
  long long int n_procs = 32;
  long long int solver = 0;
  static float dparm[64];
  std::complex<double> ddum;
  long long int idum;
  static long long int c_count = 0;
  
  //static void *pt[64];
  void *pt[64];
  long long int maxfct = 1;
  long long int mnum = 1;
  long long int mtype = 6;
  // long long int mtype = 13;
  long long int phase;
  // long long int nrhs = 1;
  long long int iparm[64];
  long long int msglvl = 0;
  long long int error = 0;
  double tmp;
  long long int memtmp;
  struct timeval s, e;
  
  if(arguments[0]){
    msglvl = 1;
  }
      
  pardisoinit(pt, &mtype, iparm);
  fprintf(stderr, "Initialization complete\n");
    
  var = getenv("OMP_NUM_THREADS");
  if (var != NULL) {
    sscanf(var, "%lld", &n_procs);
    fprintf(stderr, "OMP_NUM_THREADS is set to default: %lld\n", n_procs);
  } else {
    fprintf(stderr, "Set environment OMP_NUM_THREADS to %lld\n", n_procs);
  }
  
  //disabled
  // for(i = 0; i < 64; i++) iparm[i] = 0;
  // iparm[2] = n_procs;
  // iparm[5] = 1; // right-hand vector b is overwritten
  // iparm[17] = -1; // Report the number of non-zero elements in the factors (default=-1; disabled)
  // iparm[26] = 1; // matrix chacker(enable)
  // kokokara tsuika
  // iparm[0] = 1; // not use default values
  // iparm[1] = 3; // fill-in reducing in phase-11 (arg=3: using openmp)
  // iparm[9] = 8; // pivot perturbation(for sym matrix: arg=8)
  // iparm[10] = 1; // scaling vector(for sym matrix: usually arg=0)
  // iparm[12] = 1; // accuracy improvement by non-symmetic weighted matching (for scaling)
  // iparm[20] = 1; // Bunch-Kaufman pivoting (for scaling)
  // iparm[34]=1;
  //disabled
  
  
  iparm[0] = 1; // not use default values
  // iparm[1] = 3; // fill-in reducing in phase-11 (arg=3: using openmp)
  iparm[5] = 1; // right-hand vector b is overwritten
  // 筝�絎�絎�鐚�
  // iparm[23] = 1; // for improved parallel calculation (factorization)
  // iparm[24] = 1; // for improved parallel calculation (forward/backward solve) 
  
  iparm[26] = 1; // matrix chacker(enable)
  // iparm[30] = 0; //pertially solve?
  // iparm[33] = n_procs; // optimal OMP (CNR mode)
#ifdef PARDISO_FLOAT
  iparm[27] = 1; // single precision
#endif
    
  if(arguments[1]) iparm[59] = 1;
  
  phase = 11;
  
  fprintf(stderr, "PARDISO: Phase 1 ... "); system("date");
  
  gettimeofday(&s, NULL);
  
  pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &nn, avals, ia, ja, 
	     &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
  if (error != 0) {
    fprintf(stderr, "ERROR during symbolic factorization: %lld\n", error);
    exit(1);
  }
  gettimeofday(&e, NULL);
  
  fprintf(stderr, "Time: %.2lf min \n", 
	  ((e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec) * 1.0e-6) / 60);
  
  fprintf(stderr, "Reordering completed ... \n");
  fprintf(stderr, "Number of nonzeros in factors  = %lld\n", iparm[17]);
  //      fprintf(stderr, "Number of factorization MFLOPS = %lld\n", iparm[18]);
  fprintf(stderr, "Peak memory symbolic factorization = %lld\n", iparm[14]);
  
  /* Numerical factorization */
  phase = 22;
  iparm[26] = 1; // matrix chacker(enable)
  
  fprintf(stderr, "PARDISO: Phase 2 ... "); system("date");
  
  gettimeofday(&s, NULL);
  pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &nn, avals, ia, ja, 
	     &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
  
  if (error != 0) {
    fprintf(stderr, "ERROR during numerical factorization: %lld\n", error);
    exit(2);
  }
  gettimeofday(&e, NULL);
  fprintf(stderr, "Time: %.2lf min \n", ((e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec) * 1.0e-6) / 60);
    
  
  fprintf(stderr, "Factorization completed ... \n");
  fprintf(stderr, "Peak memory numerical factorization = %lld\n", 
	  (iparm[15] + iparm[16]));
  memtmp = (iparm[15] + iparm[16]);
  tmp = (double)memtmp;
  tmp /= 1024.0 * 1024.0;
  fprintf(stderr, "( %.1lf GiB ) \n", tmp);
  
  
  /* Back substitution and iterative refinement */
  phase = 33;
  iparm[7] = 10; // Max numbers of iterative refinement steps.
  iparm[26] = 1; // matrix chacker(enable)
  
  fprintf(stderr, "PARDISO: Phase 3 ... "); system("date");
  gettimeofday(&s, NULL);
  pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &nn, avals, ia, ja, 
	     &idum, &nrhs, iparm, &msglvl, bb, phi, &error);
  fprintf(stderr, "iterative refinement steps %lld \n", iparm[6]);
  
  if (error != 0) {
    fprintf(stderr, "ERROR during solution: %lld\n", error);
    exit(3);
  }
  gettimeofday(&e, NULL);
  fprintf(stderr, "Time: %.2lf min \n", 
	  ((e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec) * 1.0e-6) / 60);
  
  fprintf(stderr, "Solve completed ... \n");
  
  // release memory
  // phase = 0;
  phase = -1;
  pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &nn, avals, ia, ja, 
	     &idum, &nrhs, iparm, &msglvl, bb, phi, &error);

  return;
}


