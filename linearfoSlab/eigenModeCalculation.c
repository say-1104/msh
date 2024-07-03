#include "fem.h"

// function for calculating eigenvalue and eigenvector by using eigv4p

#ifndef CG
extern int eigenModeCalculation(DataTable *data, int ii)
{
  int i, j, k;
  int n, mu, ne, ne0, ialg, n1;
  FEMparam *fem = &(data->fem);
  double *eig, *eigv, *rn, *wk, eval;
  double xc_temp, xc_max;
  int *iwk, ier = 0;

  n = data->fem.np;
  mu = data->fem.nbw;
  ne = data->input.number;
  ne0 = ne+10;
  ialg = 1;
  ier = 0;
  n1 = n+1;
  eval = 1.0;

  ALLOCATION(eig, double, n)
  ALLOCATION(eigv, double, n*ne0)
  ALLOCATION(rn, double, n)
  ALLOCATION(wk, double, (3*mu+50)*n1)
  ALLOCATION(iwk, int, 3*n1)

  fprintf(stderr, "using eigv4p_ No.%d\n", ne);
  eigv4p_(fem->aa_half, fem->bb_half, &n, &mu, &ne0, &ne, &ialg, &eval,
	  eig, eigv, rn, wk, iwk, &ier);
  fprintf(stderr, "finish eigv4p_ No.%d\n", ne);

  for (i = 0; i < data->input.number; i++){
    data->par.beta[ii][i] = sqrt(1.0/eig[i]);
    data->par.n0 = data->par.beta[ii][i]/data->par.k0;
    data->par.neff[ii][i] = data->par.n0;
  }

  for (i = 0; i < data->input.number; i++){
    for (j = 0; j < data->fem.np; j++) data->fem.vv[j] = eigv[j+n*(i)];
    for (j = data->fem.nr; j < data->fem.np; j++) data->fem.vv[j] = 0.0;

    // if the max is minus, multiply -1
    xc_max = 0.0;
    for (j = 0; j < data->fem.np; j++){
      if(fabs(data->fem.vv[j]) > xc_max || j == 0){
	xc_max = data->fem.vv[j];
      }
    }

    if(xc_max < 0.0){
      for (j = 0; j < data->fem.np; j++){
	data->fem.vv[j] *= -1.0;
      }
    }

    for(j = 0; j < data->fem.np; j++){
      data->fem.field[j] = data->fem.vv[j];
      data->par.field[ii][i][j] = data->fem.vv[j];
    }

  }

  free(eig); free(eigv);
  free(rn); free(wk); 
  free(iwk); 

}
#endif
