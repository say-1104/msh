#include "fem.h"

// function for calculating effective index

#ifdef CG
extern int calcNeff(DataTable *data)
{
  int		i, j, is, ie, it;
  double_complex c1, bunbo, bunsi;
  int		nbw, nbw2;
  double  beta, neff;
  CGtable *cgt = data->fem.cg_table;
	
  for(i = 0; i < data->fem.np; i++){
    data->fem.amat[i] = data->fem.bmat[i] = 0.0;
    for (j = 0; j < cgt[i].n; j++) {
      c1 = data->fem.aa[i][j]*data->fem.field[cgt[i].col[j]];
      data->fem.amat[i] += c1;
      c1 = data->fem.Pminus[i][j]*data->fem.field[cgt[i].col[j]];
      data->fem.bmat[i] += c1;
    }
  }
  bunbo = bunsi = 0.0;
  for(i = 0; i < data->fem.np; i++){
    c1 = data->fem.amat[i]*conj(data->fem.field[i]);
    bunbo += c1;
    c1 = data->fem.bmat[i]*conj(data->fem.field[i]);
    bunsi += c1;
  }
  c1 = bunsi/bunbo;
  beta = sqrt(data->input.c_beta*data->input.c_beta-real(c1));
  neff = beta/data->par.k0;
  data->input.p_neff = data->input.c_neff;
  data->input.c_neff = data->par.n0 = neff;
  data->input.c_beta = beta;
  fprintf(stderr, "neff = %15.10lf\n", data->par.n0);
}

#else
extern int calcNeff(DataTable *data)
{
  int		i, j, is, ie, it;
  double_complex c1, bunbo, bunsi;
  int		nbw, nbw2;
  double  beta, neff;
  
  nbw = data->fem.nbw;
  nbw2 = 2*(data->fem.nbw)+1;
  for(i = 0; i < data->fem.np; i++){
    data->fem.amat[i] = data->fem.bmat[i] = 0.0;
    is = ((i-nbw) > 0)?(i-nbw):0;
    ie = ((i+nbw) < (data->fem.np-1))?(i+nbw):(data->fem.np-1);
    for(j = is; j <= ie; j++){
      it = j-i+nbw+nbw2*i;
      c1 = data->fem.aa[it]*data->fem.field[j];
      data->fem.amat[i] += c1;
      c1 = data->fem.Pminus[it]*data->fem.field[j];
      data->fem.bmat[i] += c1;
    }
  }
  bunbo = bunsi = 0.0;
  for(i = 0; i < data->fem.np; i++){
    c1 = data->fem.amat[i]*conj(data->fem.field[i]);
    bunbo += c1;
    c1 = data->fem.bmat[i]*conj(data->fem.field[i]);
    bunsi += c1;
  }
  c1 = bunsi/bunbo;
  beta = sqrt(data->input.c_beta*data->input.c_beta-real(c1));
  neff = beta/data->par.k0;
  data->input.p_neff = data->input.c_neff;
  data->input.c_neff = data->par.n0 = neff;
  data->input.c_beta = beta;
  fprintf(stderr, "neff = %15.10lf\n", data->par.n0);

}
#endif
