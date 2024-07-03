#include "fem.hpp"

extern int propagation(DataTable *data)
{
  int		i, j, is, ie, it, ier;
  int		nbw, nbw2;
  double_complex	c1, c2;
  double beta2 = data->input.c_beta*data->input.c_beta, eps;   

#ifdef CG
  CGtable *cgt = data->fem.cg_table;
  
  for (i = 0; i < data->fem.np; i++) {
    for (j = 0; j < cgt[i].n; j++) {
      data->fem.Pplus[i][j] = beta2*data->fem.aa[i][j]
	-data->fem.bb[i][j];
      data->fem.Pminus[i][j] = beta2*data->fem.aa[i][j]
	-data->fem.bb[i][j];
    }
  }
#else
  ie = (data->fem.nbw*2+1)*data->fem.np;
  
  for(i = 0; i < ie; i++){
    data->fem.Pplus[i] = beta2*data->fem.aa[i]-data->fem.bb[i];
    data->fem.Pminus[i] = beta2*data->fem.aa[i]-data->fem.bb[i];        
  }

#endif


#ifdef CG
  /* LU decomposition */
  multifrontalF(data->fem.Pplus, data->fem.fieldn, data->fem.np,
					data->fem.field, cgt);

#else

  /*	µÍ¤áÄ¾¤·	*/
  reArrangment(data);

  eps = 1.0e-16;

  /* LU decomposition */
  bhlu1c_(data->fem.Pplus, data->fem.cac,
	  &(data->fem.np), &(data->fem.nbw), &(data->fem.nbw), &eps,
	  data->fem.cwk, data->fem.iwp, &ier);
#endif
  for(i = 0; i < data->fem.np; i++)
    data->fem.fieldn[i] = 0.0;

#ifdef CG
  for(i = 0; i < data->fem.np; i++){
    for(j = 0; j < cgt[i].n; j++){
      data->fem.fieldn[i]
	+= data->fem.aa[i][j]*data->fem.field[cgt[i].col[j]];
    }
  }
#else
  nbw = data->fem.nbw;
  nbw2 = 2*(data->fem.nbw)+1;
  for(i = 0; i < data->fem.np; i++){
    is = ((i-nbw) > 0)?(i-nbw):0;
    ie = ((i+nbw) < (data->fem.np-1))?(i+nbw):(data->fem.np-1);
    for(j = is; j <= ie; j++){
      it = (j-i+nbw)+nbw2*i;
      data->fem.fieldn[i] += data->fem.aa[it]*data->fem.field[j];
    }
  }
  
#endif

#ifdef CG

  /* forward and backward substitution */
  multifrontalS(data->fem.Pplus, data->fem.fieldn, data->fem.np,
		data->fem.field, cgt);
  
  fprintf(stderr, "free -------------\n");
  multifrontalFree();

  for(i = 0; i < data->fem.np; i++)
    data->fem.field[i] = abs(data->fem.fieldn[i]);
  
#else
  bhslv1c_(data->fem.Pplus, data->fem.cac,
	   &(data->fem.np), &(data->fem.nbw), &(data->fem.nbw),
	   data->fem.fieldn, data->fem.iwp);

  for(i = 0; i < data->fem.np; i++)
    data->fem.field[i] = abs(data->fem.fieldn[i]);

#endif
}
