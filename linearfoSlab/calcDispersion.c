#include "fem.h"

// function for calculating chromatic dispersion
// only active if inverse iterative method is used

extern int calcDispersion(DataTable *data)
{
  int i, j, k, tmp;
  double *D1, *D2;
  double delta, temp, micro = 1.0e-6;
  FILE *fp;

  fp = fopen("Dispersion", "w");
  
  tmp = (int)((data->par.lamda2-data->par.lamda1)/data->par.step);
  
  ALLOCATION(D1, double, tmp-1)
  ALLOCATION(D2, double, tmp-2)

  for(i = 0; i < tmp; i++){
    D1[i] = (data->par.neff[i+1]-data->par.neff[i])/(data->par.step);
  }

  for(i = 0; i < tmp-1; i++){
    D2[i] = (D1[i+1]-D1[i])/(data->par.step);
  }

  for(i = 0; i < tmp-1; i++){
    D2[i] = -1.0*1.0e6*D2[i]*data->par.wl[i+1]/(300.0);
  }

  for(i = 0; i < tmp-1; i++){
    fprintf(fp, "%lf %15.10lf\n", data->par.wl[i+1], D2[i]);
  }
  fclose(fp);
}
