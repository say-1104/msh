#include <cstdio>
#include <cstdlib>
#include <complex>
#include "fem3d.h"

void pwsmp_pre(std::complex<double> **A, std::complex<double> *bb, 
	       long long int nf, DataTable *data, long long int ns, 
	       long long int ne, long long int nrhs) {
  long long int nn = ne - ns;
  long long int i, j;
  long long int count;
    
  CGtable *cg_table = data->fem.cg_table;
    
  count = 0;
  for (i = ns; i < ne; i++) {
    for (j = 0; j < cg_table[i].n; j++) {
      if (cg_table[i].col[j] >= i)
	count++;
    }
  }
    
  data->matrix.nn = ne - ns;
  data->matrix.nc = count;
    
  ALLOCATION(data->matrix.ia, long long int, nn+1);
  ALLOCATION(data->matrix.ja, long long int, count);
  ALLOCATION(data->matrix.avals, std::complex<double>, count);
  ALLOCATION(data->matrix.b, std::complex<double>, nn*nrhs);

  count = 0;
  for (i = ns; i < ne; i++) {
    data->matrix.ia[i - ns] = count + 1;
    for (j = 0; j < cg_table[i].n; j++) {
      if (cg_table[i].col[j] >= i) {
	data->matrix.ja[count] = cg_table[i].col[j] + 1;
	data->matrix.avals[count] = A[i][j];
	count++;
      }
    }
    for (j = 0; j < nrhs; j++) {
      data->matrix.b[j*nf + i - ns] = bb[j*nf + i];
    }
  }
  data->matrix.ia[nn] = count + 1;
  
  return;
}

long long int wsmp_pre(std::complex<double> **A, std::complex<double> *bb, 
		       long long int nf, DataTable *data, 
		       long long int nrhs) {
  long long int ne = nf;
  long long int ns = 0;
  long long int nn = ne - ns;
  long long int i, j;
  long long int count;
    
  CGtable *cg_table = data->fem.cg_table;
  
  count = 0;
  // wsmp
  /*
    for (i = ns; i < ne; i++) {
    for (j = 0; j < cg_table[i].n; j++) {
    if (cg_table[i].col[j] >= i)
    count++;
    }
    }
  */
  // pardiso
  
  for (i = 0, count = 0; i < nn; i++) {
    for (j = 0; j < cg_table[i].n; j++) {
      count++;
    }
  }
  
  data->matrix.nn = ne - ns;
  data->matrix.nc = count;
    
  ALLOCATION(data->matrix.ia, long long int, nn+1);
  ALLOCATION(data->matrix.ja, long long int, count);
  ALLOCATION(data->matrix.avals, std::complex<double>, count);
  ALLOCATION(data->matrix.b, std::complex<double>, nn*nrhs);

  count = 0;
  for (i = ns; i < ne; i++) {
    data->matrix.ia[i - ns] = count + 1;
    for (j = 0; j < cg_table[i].n; j++) {
      if (cg_table[i].col[j] >= i) {
	data->matrix.ja[count] = cg_table[i].col[j] + 1;
	data->matrix.avals[count] = A[i][j];
	count++;
      }
    }
    for (j = 0; j < nrhs; j++) {
      data->matrix.b[j*nf + i - ns] = bb[j*nf + i];
    }
  }
  data->matrix.ia[nn] = count + 1;
    
  return 0;
}
