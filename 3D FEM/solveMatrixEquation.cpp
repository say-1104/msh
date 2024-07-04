#include "fem3d.h"

void readField(DataTable *data) {
  long long int i, j;
  FILE *fp;
  double re, im;
  std::complex<double> cj = std::complex<double>(0, 1);
  std::complex<double> *cwkG = NULL;
    
  FEM *fem = &(data->fem);
    
    
  fprintf(stderr, "nr: %lld \n",fem->nr_edge+fem->nr);
  if ((fp = fopen("FIELD", "r")) != NULL) {
    ALLOCATION(cwkG, std::complex<double>, fem->nr_edge+fem->nr);
    for (i = 0; i < fem->nr_edge + fem->nr; i++) {
      fscanf(fp, "%lf %lf", &re, &im);
      cwkG[i] = re + cj * im;
    }
  }
  fclose(fp);
    
  for (i = 0; i < fem->nr_edge; i++) {
    fem->phi[i] = cwkG[i];
  }
    
  for (i = 0; i < fem->nr_edge; i++) {
    j = fem->Gtable[0][i];
    if (j != -1)
      fem->phi[i] += fem->G[0][i] * cwkG[j + fem->nr_edge];
    j = fem->Gtable[1][i];
    if (j != -1)
      fem->phi[i] += fem->G[1][i] * cwkG[j + fem->nr_edge];
  }
  
  return;
}

void solveMatrixEquation(DataTable *data) {
  long long int i, j;
  // long long int nrhs = data->port[0].modenum;
  std::complex<double> *B = NULL;
  
  FEM *fem = &(data->fem);
  
  // B作らなくてもいい気がするが？
  ALLOCATION(B, std::complex<double>, fem->nr_edge * data->nrhs);
    
  for (j = 0; j < data->nrhs; j++) {
    for (i = 0; i < fem->nr_edge; i++) {
      B[j*fem->nr_edge +i] = fem->phi[j*fem->nr_edge + i];
    }
  }
  fprintf(stderr, "solveMatrixEquation() : wsmp_pre()\n"); 
  // [A]行列を2次元から1次元に，cgからia,jaに
  wsmp_pre(fem->ok, B, fem->nr_edge, data, data->nrhs);

  return;
}
