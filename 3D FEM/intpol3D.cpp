#include "fem3d.h"
#define DIV_X 20
#define DIV_Y 40
#define DIV_Z 200
#define NEARLY_ZERO 1.0e-6


CXYZdouble intpol3D(DataTable *data, long long int n_element, 
		    long long int n_point) {
  CXYZdouble phi;
  long long int i, k, l, m, n;
  long long int itmp;
  double a[4], b[4], c[4], d[4], alpha[4];
  double U[24], V[24], W[24];
  double ss[6], aa[6], bb[6], cc[6], ll[6];
  static long long int edge[6][2] = { { 0, 1 }, { 0, 2 }, { 0, 3 }, 
				      { 1, 2 }, { 3, 1 }, { 2, 3 } };
  double Ve;
  double L1, L2, L3, L4;
  double L1x, L2x, L3x, L4x;
  double L1y, L2y, L3y, L4y;
  double L1z, L2z, L3z, L4z;
  double x[10], y[10], z[10];
  std::complex<double> vv[24];
  double x0, y0, z0;
  
  double area[4];
  double co[24];
  
  FEM *fem = &(data->fem);
  
  phi.x = phi.y = phi.z = 0.0;
  
  alpha[0] = alpha[2] = 1.0;
  alpha[1] = alpha[3] = -1.0;
  
  x0 = fem->node[n_point].x;
  y0 = fem->node[n_point].y;
  z0 = fem->node[n_point].z;
  
  for (i = 0; i < 4; i++) {
    itmp = fem->element[n_element].kn[i];
    x[i] = fem->node[itmp].x;
    y[i] = fem->node[itmp].y;
    z[i] = fem->node[itmp].z;
  }
  
  for (i = 0; i < 4; i++) {
    k = i;
    l = (i + 1) % 4;
    m = (i + 2) % 4;
    n = (i + 3) % 4;
    a[k] = alpha[k] * (x[l] * (y[m] * z[n] - y[n] * z[m]) 
		       + x[m] * (y[n] * z[l] - y[l] * z[n]) 
		       + x[n] * (y[l] * z[m] - y[m] * z[l]));
    b[k] = alpha[k] * (y[l] * (z[n] - z[m]) 
		       + y[m] * (z[l] - z[n]) + y[n] * (z[m] - z[l]));
    c[k] = alpha[k] * (z[l] * (x[n] - x[m]) 
		       + z[m] * (x[l] - x[n]) + z[n] * (x[m] - x[l]));
    d[k] = alpha[k] * (x[l] * (y[n] - y[m]) 
		       + x[m] * (y[l] - y[n]) + x[n] * (y[m] - y[l]));
  }
  Ve = (a[0] + a[1] + a[2] + a[3]) / 6.0;
  
  L1 = (a[0] + b[0] * x0 + c[0] * y0 + d[0] * z0) / (6.0 * Ve);
  L2 = (a[1] + b[1] * x0 + c[1] * y0 + d[1] * z0) / (6.0 * Ve);
  L3 = (a[2] + b[2] * x0 + c[2] * y0 + d[2] * z0) / (6.0 * Ve);
  L4 = (a[3] + b[3] * x0 + c[3] * y0 + d[3] * z0) / (6.0 * Ve);

  if ((L1 >= -NEARLY_ZERO) && (L2 >= -NEARLY_ZERO) 
      && (L3 >= -NEARLY_ZERO)
      && (L4 >= -NEARLY_ZERO)) {

    for (i = 0; i < 6; i++) {
      aa[i] = x[edge[i][1]] - x[edge[i][0]];
      bb[i] = y[edge[i][1]] - y[edge[i][0]];
      cc[i] = z[edge[i][1]] - z[edge[i][0]];
      ll[i] = sqrt(aa[i] * aa[i] + bb[i] * bb[i] + cc[i] * cc[i]);
    }
    
    area[0] = sqrt(pow(bb[3] * cc[4] - bb[4] * cc[3], 2.0) 
		   + pow(cc[3] * aa[4] - cc[4] * aa[3], 2.0) 
		   + pow(aa[3] * bb[4] - aa[4] * bb[3], 2.0));
    area[1] = sqrt(pow(bb[1] * cc[2] - bb[2] * cc[1], 2.0) 
		   + pow(cc[1] * aa[2] - cc[2] * aa[1], 2.0) 
		   + pow(aa[1] * bb[2] - aa[2] * bb[1], 2.0));
    area[2] = sqrt(pow(bb[0] * cc[2] - bb[2] * cc[0], 2.0) 
		   + pow(cc[0] * aa[2] - cc[2] * aa[0], 2.0) 
		   + pow(aa[0] * bb[2] - aa[2] * bb[0], 2.0));
    area[3] = sqrt(pow(bb[0] * cc[1] - bb[1] * cc[0], 2.0) 
		   + pow(cc[0] * aa[1] - cc[1] * aa[0], 2.0) 
		   + pow(aa[0] * bb[1] - aa[1] * bb[0], 2.0));
    
    for (i = 0; i < 6; i++) {
      if (cc[i] > 0.0 || (cc[i] == 0.0 && bb[i] > 0.0) 
	  || (cc[i] == 0.0 && bb[i] == 0.0 && aa[i] > 0.0)) {
	ss[i] = 1.0;
      } else {
	ss[i] = -1.0;
      }
    }

    L1x = b[0] / (6.0 * Ve);
    L1y = c[0] / (6.0 * Ve);
    L1z = d[0] / (6.0 * Ve);
    L2x = b[1] / (6.0 * Ve);
    L2y = c[1] / (6.0 * Ve);
    L2z = d[1] / (6.0 * Ve);
    L3x = b[2] / (6.0 * Ve);
    L3y = c[2] / (6.0 * Ve);
    L3z = d[2] / (6.0 * Ve);
    L4x = b[3] / (6.0 * Ve);
    L4y = c[3] / (6.0 * Ve);
    L4z = d[3] / (6.0 * Ve);
    if (data->fem.order == 1) {
      
      U[0] = L1 * L2x - L2 * L1x;
      V[0] = L1 * L2y - L2 * L1y;
      W[0] = L1 * L2z - L2 * L1z;
      U[1] = L1 * L3x - L3 * L1x;
      V[1] = L1 * L3y - L3 * L1y;
      W[1] = L1 * L3z - L3 * L1z;
      U[2] = L1 * L4x - L4 * L1x;
      V[2] = L1 * L4y - L4 * L1y;
      W[2] = L1 * L4z - L4 * L1z;
      U[3] = L2 * L3x - L3 * L2x;
      V[3] = L2 * L3y - L3 * L2y;
      W[3] = L2 * L3z - L3 * L2z;
      U[4] = L4 * L2x - L2 * L4x;
      V[4] = L4 * L2y - L2 * L4y;
      W[4] = L4 * L2z - L2 * L4z;
      U[5] = L3 * L4x - L4 * L3x;
      V[5] = L3 * L4y - L4 * L3y;
      W[5] = L3 * L4z - L4 * L3z;
      
      for (i = 0; i < 6; i++) {
	U[i] *= ss[i] * ll[i];
	V[i] *= ss[i] * ll[i];
	W[i] *= ss[i] * ll[i];
      }
    } else {
      co[12] = 4.0 * area[0] / ll[5];
      co[13] = 4.0 * area[0] / ll[4];
      co[14] = 4.0 * area[0] / ll[3];
      co[15] = 4.0 * area[1] / ll[5];
      co[16] = 4.0 * area[1] / ll[1];
      co[17] = 4.0 * area[1] / ll[2];
      co[18] = 4.0 * area[2] / ll[0];
      co[19] = 4.0 * area[2] / ll[4];
      co[20] = 4.0 * area[2] / ll[2];
      co[21] = 4.0 * area[3] / ll[0];
      co[22] = 4.0 * area[3] / ll[1];
      co[23] = 4.0 * area[3] / ll[3];
      
      U[0] = L1 * L2x;
      V[0] = L1 * L2y;
      W[0] = L1 * L2z;
      U[1] = L1 * L3x;
      V[1] = L1 * L3y;
      W[1] = L1 * L3z;
      U[2] = L1 * L4x;
      V[2] = L1 * L4y;
      W[2] = L1 * L4z;
      U[3] = L2 * L3x;
      V[3] = L2 * L3y;
      W[3] = L2 * L3z;
      U[4] = L4 * L2x;
      V[4] = L4 * L2y;
      W[4] = L4 * L2z;
      U[5] = L3 * L4x;
      V[5] = L3 * L4y;
      W[5] = L3 * L4z;
      
      U[6] = L1x * L2;
      V[6] = L1y * L2;
      W[6] = L1z * L2;
      U[7] = L1x * L3;
      V[7] = L1y * L3;
      W[7] = L1z * L3;
      U[8] = L1x * L4;
      V[8] = L1y * L4;
      W[8] = L1z * L4;
      U[9] = L2x * L3;
      V[9] = L2y * L3;
      W[9] = L2z * L3;
      U[10] = L4x * L2;
      V[10] = L4y * L2;
      W[10] = L4z * L2;
      U[11] = L3x * L4;
      V[11] = L3y * L4;
      W[11] = L3z * L4;
      
      U[12] = L3 * L4 * L2x;
      V[12] = L3 * L4 * L2y;
      W[12] = L3 * L4 * L2z;
      U[13] = L4 * L2 * L3x;
      V[13] = L4 * L2 * L3y;
      W[13] = L4 * L2 * L3z;
      U[14] = L2 * L3 * L4x;
      V[14] = L2 * L3 * L4y;
      W[14] = L2 * L3 * L4z;
      U[15] = L4 * L3 * L1x;
      V[15] = L4 * L3 * L1y;
      W[15] = L4 * L3 * L1z;
      U[16] = L3 * L1 * L4x;
      V[16] = L3 * L1 * L4y;
      W[16] = L3 * L1 * L4z;
      U[17] = L1 * L4 * L3x;
      V[17] = L1 * L4 * L3y;
      W[17] = L1 * L4 * L3z;
      U[18] = L1 * L2 * L4x;
      V[18] = L1 * L2 * L4y;
      W[18] = L1 * L2 * L4z;
      U[19] = L2 * L4 * L1x;
      V[19] = L2 * L4 * L1y;
      W[19] = L2 * L4 * L1z;
      U[20] = L4 * L1 * L2x;
      V[20] = L4 * L1 * L2y;
      W[20] = L4 * L1 * L2z;
      U[21] = L2 * L1 * L3x;
      V[21] = L2 * L1 * L3y;
      W[21] = L2 * L1 * L3z;
      U[22] = L1 * L3 * L2x;
      V[22] = L1 * L3 * L2y;
      W[22] = L1 * L3 * L2z;
      U[23] = L3 * L2 * L1x;
      V[23] = L3 * L2 * L1y;
      W[23] = L3 * L2 * L1z;
      
      for (i = 0; i < 12; i++) {
	U[i] *= ll[i % 6];
	V[i] *= ll[i % 6];
	W[i] *= ll[i % 6];
      }
      for (i = 12; i < fem->n_en; i++) {
	U[i] *= co[i];
	V[i] *= co[i];
	W[i] *= co[i];
      }
    }
    
    for (i = 0; i < data->fem.n_en; i++) {
      itmp = fem->element[n_element].kk[i];
      if (itmp >= 0) {
	vv[i] = fem->phi[itmp];
      } else {
	vv[i] = 0.0;
      }
    }
    for (i = 0; i < data->fem.n_en; i++) {
      phi.x += U[i] * vv[i];
      phi.y += V[i] * vv[i];
      phi.z += W[i] * vv[i];
    }
    // fprintf(stderr, "%10.6lf %10.6lf %10.6e\n", x0, z0, std::abs(phix));
  }
  
  return phi;
}


void outputField_for_Gid(DataTable *data) {
  long long int i, j;
  long long int num;
  long long int *flag = NULL;
  CXYZdouble *phi = NULL;
  FILE *fp;
  char string[256];
  
  FEM *fem = &(data->fem);
  
  ALLOCATION(flag, long long int, fem->np);
  ALLOCATION(phi, CXYZdouble, fem->np);

  for (i = 0; i < fem->np; i++) {
    flag[i] = 0;
    phi[i].x = phi[i].y = phi[i].z = 0.0;
  }
  
  for (i = 0; i < fem->ne; i++) {
    for (j = 0; j < fem->n_node; j++) {
      num = fem->element[i].kn[j];
      if (flag[num] == 0) {
	phi[num] = intpol3D(data, i, num);
	flag[num] = 1;
      }
    }
  }
  
  // sprintf(string, "%3.3lf.flavia.res", par->parameter);
  sprintf(string, "%s.flavia.res", data->name);
  fprintf(stderr, "creating file : %s.flavia.res\n", data->name);
  
  if ((fp = fopen(string, "w")) == NULL) {
    fprintf(stderr, "can't open file (%s).\n", string);
    exit(EXIT_FAILURE);
  }
  
  fprintf(fp, "Electric-Field  3 0.00E+00 2 1 0\n");
  
  fprintf(stderr, "num = %lld\n", fem->np);
  for (i = 0; i < fem->np; i++) {
    fprintf(fp, "%lld %10.6e %10.6e %10.6e\n", i + 1, std::abs(phi[i].x),
	    std::abs(phi[i].y), std::abs(phi[i].z));
  }
  
  if (fclose(fp) == EOF) {
    fprintf(stderr, "can't close file ( %s )\n", string);
    exit(EXIT_FAILURE);
  }
  
  // sprintf(string, "%3.3lf.flavia.msh", par->parameter);
  sprintf(string, "%s.flavia.msh", data->name);
  fprintf(stderr, "creating file : %s.flavia.msh\n", data->name);
  
  if ((fp = fopen(string, "w")) == NULL) {
    fprintf(stderr, "can't open file (%s).\n", string);
    exit(EXIT_FAILURE);
  }
  
  fprintf(fp, "MESH Dimension 3 ElemType Tetrahedra Nnode 10\n");
  fprintf(fp, "Coordinates\n");
  
  for (i = 0; i < fem->np; i++) {
    fprintf(fp, "  %lld %10.6E %10.6E %10.6E\n", i + 1, fem->node[i].x,
	    fem->node[i].y, fem->node[i].z);
  }
  fprintf(fp, "End Coordinates \n");
  fprintf(fp, "Elements \n");
  
  fprintf(stderr, "num = %lld\n", fem->ne);
  for (i = 0; i < fem->ne; i++) {
    fprintf(fp, "%lld ", i + 1);//%5d
    for (j = 0; j < fem->n_node; j++) {
      fprintf(fp, "%lld ", fem->element[i].kn[j] + 1);//%5d
    }
    fprintf(fp, "%lld", fem->element[i].matID + 1);//%5d
    fprintf(fp, "\n");
  }
  
  fprintf(fp, " End Elements\n");
  
  if (fclose(fp) == EOF) {
    fprintf(stderr, "can't close file ( %s )\n", string);
    exit(EXIT_FAILURE);
  }
  
  return;
}
