#include "optfem.h"
#include <cstdio>
#include <cstdlib>
#include <complex>

#include <mkl_solvers_ee.h>
#include <mkl_pardiso.h>

extern int callPARDISO(DataTable *data, std::complex<double> **A, 
		       std::complex<double> *bb, int nn, 
		       std::complex<double> *phi, CGtable *cg_table, 
		       int c_count) {
    int i, j;
    int count;
    char *var;
    int n_procs = 4;
    int solver = 0;
    static float dparm[64];
    std::complex<double> ddum;
    int idum;
    // static int c_count = 0;

    static void *pt[64];
    int maxfct = 1;
    int mnum = 1;
    int mtype;
    
    int phase;
    static std::complex<double> *avals;
    static int *ia;
    static int *ja;
    int nrhs = 1;
    static int iparm[64];
    int msglvl = 0;
    int error = 0;
    
    
    if (c_count < 0) {
      fprintf(stderr, "PARDISO free ... \n");
      
      // release memory
      // phase = 0;
      phase = -1;
      pardiso(pt, &maxfct, &mnum, &mtype, &phase, &nn, avals, ia, ja, 
	      &idum, &nrhs, iparm, &msglvl, bb, phi, &error);

      free(ia);
      free(ja);
      free(avals);
      return 0;
    }
    
    //    if (data->flag.anisotropy == 0) {
    if (0) {

        //  2: Real and symmetric positive definite
        // -2: Real and symmetric indefinite
        //  3: Complex and structurally symmetric
        //  4: Complex and Hermitian positive definite
        // -4: Complex and Hermitian indefinite
        //  6: Complex and symmetric matrix
        // 11: Real and unsymmetric matrix
        // 13: Complex and unsymmetric matrix
        if (c_count == 0) fprintf(stderr, "(Complex and Symmetric matrix)\n");
        mtype = 6;
    }else {
        if (c_count == 0) fprintf(stderr, "(Complex and Unsymmetric matrix)\n");
        mtype = 13;
    }
    // iparm[0] = 0; // iparm(2) - iparm(64) are filled with default values.

    // iparm[0] = 1; // You must supply all values in components iparm(2) - iparm(64).
    // iparm[1] = 2; // Fill-in reducing ordering for the input matrix.
    // iparm[2] = 0; // Reserved. Set to zero.
    // iparm[3] = 61; // Preconditioned CGS/CG.

 
    if (c_count == 0) {
    
        /* Memory allocation */
        for (i = 0, count = 0; i < nn; i++) {
            for (j = 0; j < cg_table[i].n; j++) {
                count++;
            }
        }

		ALLOCATION(ia, int, nn+1);
		ALLOCATION(ja, int, count);
		ALLOCATION(avals, std::complex<double>, count);
    
        /* Setup Pardiso control parameters and initialize the solvers */
        /* internal adress pointers */
    
        pardisoinit(pt, &mtype, iparm); // iparm[64]をデフォルト値にset
        fprintf(stderr, "Initialization complete\n");
    
                
        // Numbers of processors, value of OMP_NUM_THREADS */
        var = getenv("OMP_NUM_THREADS");
        if (var != NULL) {
        sscanf(var, "%d", &n_procs);
        } else {
        fprintf(stderr, "Set environment OMP_NUM_THREADS to 1\n");
        }
        
        iparm[2] = n_procs;
    
        /* Fill-in reduction reordering */
        //      iparm[1] = 0;
    
        /* Direct-iterative preconditioning */
        //      iparm[3] = 61;

        /* Argument substitution */
        count = 0;
        for (i = 0; i < nn; i++) {
            ia[i] = count + 1;
            for (j = 0; j < cg_table[i].n; j++) {
	      // if (data->flag.anisotropy == 0) {
                if (0) {
                    if (cg_table[i].col[j] >= i) {
                        ja[count] = cg_table[i].col[j] + 1;
                        avals[count] = A[i][j];
                        count++;
                    }
                }else {
                    ja[count] = cg_table[i].col[j] + 1;
                    avals[count] = A[i][j];
                    count++;
                }
            }
        }
        ia[nn] = count + 1;

    
        /* Reordering and Symbolic Factorization. This step also allocates */
        /* all memory that is necessary for the factorization. */
        phase = 11;
    
        pardiso(pt, &maxfct, &mnum, &mtype, &phase, &nn, avals, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
        if (error != 0) {
            fprintf(stderr, "ERROR during symbolic factorization: %d\n", error);
            exit(1);
        }
        fprintf(stderr, "Reordering completed ... \n");
        fprintf(stderr, "Number of nonzeros in factors  = %d\n", iparm[17]);
        // fprintf(stderr, "Number of factorization MFLOPS = %d\n", iparm[18]);
        fprintf(stderr, "Peak memory symbolic factorization = %d\n", iparm[14]);
    
        /* Numerical factorization */
        phase = 22;
    
        pardiso(pt, &maxfct, &mnum, &mtype, &phase, &nn, avals, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    
        if (error != 0) {
            fprintf(stderr, "ERROR during numerical factorization: %d\n", error);
            exit(2);
        }
        fprintf(stderr, "Factorization completed ... \n");
        fprintf(stderr, "Peak memory numerical factorization = %d\n", (iparm[15] + iparm[16]));
    
    }

    /* Back substitution and iterative refinement */
    phase = 33;
    iparm[7] = 16; // Max numbers of iterative refinement steps.
    iparm[19] = 1; // examination

    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &nn, avals, ia, ja, &idum, &nrhs, iparm, &msglvl, bb, phi, &error);

    if (error != 0) {
        fprintf(stderr, "ERROR during solution: %d\n", error);
        exit(3);
    }

    fprintf(stderr, "Solve completed ... ");
    fprintf(stderr, "(refinement %d -> %d, exam(%d)\n", iparm[7], iparm[6], iparm[19]);

    return 0;
}
