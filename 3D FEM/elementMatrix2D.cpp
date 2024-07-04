#include "fem3d.h"

#define LEE

#define AG1 0.05971587
#define BG1 0.47014206
#define AG2 0.79742669
#define BG2 0.10128651

long long int checkLongEdge(double *x1, double *y1, double al[]);

void setAreaCoordinate2d(long long int i, double L[], double *W);

void setAreaCoordinateDiff2d(double jj[][2], double Lx[], double Ly[]);

void setMaterialInfo2d(DataTable *data, long long int matID, 
		       double *pp, double *qq);

void setNodalShapeFunc2d(DataTable *data, double L[], double N[], double nx[],
			 double ny[]);

void setNodalShapeFuncDiff2d(DataTable *data, double jj[][2], double nx[],
			     double ny[], double Nx[], double Ny[]);

void setJacobianMatrix2d(DataTable *data, double *x1, double *y1, double nx[],
			 double ny[], double jj[][2]);

void setEdgeShapeFunc2d(DataTable *data, double L[], double Lx[], double Ly[],
			double jj[][2], double al[], long long int kaki, 
			double ss[], double U[],
			double V[]);

void setEdgeShapeFuncDiff2d(DataTable *data, double L[], double Lx[],
			    double Ly[], double jj[][2], double al[], 
			    long long int kaki, double ss[],
			    double Uy[], double Vx[]);

void elementMatrix2D(DataTable *data, double *x1, double *y1, 
		     long long int matID,
		     std::complex<double> BB[N_EN][N_EN],
		     std::complex<double> BBeh[N_EN][N_EN]) {
  long long int i, j, l, j1, j2;
  double nx[N_NODE], ny[N_NODE];
  double N[N_NODE], Nx[N_NODE], Ny[N_NODE];
  double U[N_EDGE], V[N_EDGE], Uy[N_EDGE], Vx[N_EDGE];
  double UU[N_EDGE][N_EDGE], VV[N_EDGE][N_EDGE];
  // double UyUy[N_EDGE][N_EDGE], VxVx[N_EDGE][N_EDGE];
  // double VxUy[N_EDGE][N_EDGE], UyVx[N_EDGE][N_EDGE];
  double UV[N_EDGE][N_EDGE], VU[N_EDGE][N_EDGE];
  double UNx[N_EDGE][N_NODE], VNy[N_EDGE][N_NODE];
  // double NN[N_NODE][N_NODE], NxNx[N_NODE][N_NODE], NyNy[N_NODE][N_NODE];
  
  double b[3], c[3], ll[3], ss[3];
  double al[3];
  double W;
  double L[3];
  double Lx[3], Ly[3];
  // double k0;
  // double px, py, pz, qx, qy, qz;
  // double s12, s23, s31;
  double pp, qq;
  
  // std::complex<double> cj = std::complex<double>(0, 1);
  
  long long int l_edge = 8, l_node = 6;
  
  double jj[2][2];
  double dj;
  long long int kaki;
  
  FEM *fem = &(data->fem);
  
  switch (fem->order) {
  case 1:
    l_edge = 3;
    l_node = 3;
    break;
  case 2:
    l_edge = 6;
    l_node = 6;
    break;
  case 3:
    l_edge = 8;
    l_node = 6;
    break;
  }
  
  for (i = 0; i < N_NODE; i++)
    N[i] = Nx[i] = Ny[i] = nx[i] = ny[i] = 0.0;
  for (i = 0; i < N_EDGE; i++)
    U[i] = V[i] = Uy[i] = Vx[i] = 0.0;
  
  for (i = 0; i < N_EN; i++) {
    for (j = 0; j < N_EN; j++) {
      BB[i][j] = 0.0;
      BBeh[i][j] = 0.0;
    }
  }
  
  // k0 = 2.0 * M_PI / par->wavelength;
  kaki = checkLongEdge(x1, y1, al);
  
  create2dElementParam(x1, y1, b, c, ll, ss);
  /*
    if (fem->order == 1) {
    for (i = 0; i < 3; i++) {
    if (b[i] < -1.0e-4 || (fabs(b[i]) < 1.0e-4 && c[i] > 0.0)) {
    ss[i] = 1.0;
    } else {
    ss[i] = -1.0;
    }
    }
    }
    s23 = ss[0]; s31 = ss[1]; s12 = ss[2];
  */
  
  for (l = 0; l < 7; l++) {
    setMaterialInfo2d(data, matID, &pp, &qq);
    
    setAreaCoordinate2d(l, L, &W);
    setNodalShapeFunc2d(data, L, N, nx, ny);
    
    setJacobianMatrix2d(data, x1, y1, nx, ny, jj);
    dj = jj[0][0] * jj[1][1] - jj[0][1] * jj[1][0];
    
    setAreaCoordinateDiff2d(jj, Lx, Ly);
    setNodalShapeFuncDiff2d(data, jj, nx, ny, Nx, Ny);
    
    setEdgeShapeFunc2d(data, L, Lx, Ly, jj, al, kaki, ss, U, V);
    setEdgeShapeFuncDiff2d(data, L, Lx, Ly, jj, al, kaki, ss, Uy, Vx);
    
    for (j1 = 0; j1 < l_edge; ++j1) {
      for (j2 = 0; j2 < l_edge; ++j2) {
	UU[j1][j2] = dj * (0.5 * W) * U[j1] * U[j2];
	VV[j1][j2] = dj * (0.5 * W) * V[j1] * V[j2];
	UV[j1][j2] = dj * (0.5 * W) * U[j1] * V[j2];
	VU[j1][j2] = dj * (0.5 * W) * V[j1] * U[j2];
	/*
	  UyUy[j1][j2] = dj*(0.5*W)*Uy[j1]*Uy[j2];
	  VxVx[j1][j2] = dj*(0.5*W)*Vx[j1]*Vx[j2];
	  UyVx[j1][j2] = dj*(0.5*W)*Uy[j1]*Vx[j2];
	  VxUy[j1][j2] = dj*(0.5*W)*Vx[j1]*Uy[j2];
	*/
      }
    }
    
    for (j1 = 0; j1 < l_edge; ++j1) {
      for (j2 = 0; j2 < l_node; ++j2) {
	VNy[j1][j2] = dj * (0.5 * W) * V[j1] * Ny[j2];
	UNx[j1][j2] = dj * (0.5 * W) * U[j1] * Nx[j2];
      }
    }
    
    /*
      for (j1 = 0; j1 < l_node; ++j1) {
      for (j2 = 0; j2 < l_node; ++j2) {
      NN[j1][j2] = dj*(0.5*W)*N[j1]*N[j2];
      NxNx[j1][j2] = dj*(0.5*W)*Nx[j1]*Nx[j2];
      NyNy[j1][j2] = dj*(0.5*W)*Ny[j1]*Ny[j2];
      }
      }
    */
    
    for (j1 = 0; j1 < l_edge; ++j1) {
      for (j2 = 0; j2 < l_edge; ++j2) {
	BB[j1][j2] += pp * (UU[j1][j2] + VV[j1][j2]);
      }
    }
    for (j1 = 0; j1 < l_edge; ++j1) {
      for (j2 = 0; j2 < l_node; ++j2) {
	if (fem->order == 1) {
	  BB[j1][j2 + l_edge] += pp * (UNx[j1][j2] + VNy[j1][j2]);
	} else {
	  BB[j1][j2 + l_edge] += pp * (UNx[j1][j2] + VNy[j1][j2]);
	}
      }
    }
    
    for (j1 = 0; j1 < l_edge; ++j1) {
      for (j2 = 0; j2 < l_edge; ++j2) {
	BBeh[j1][j2] += (UV[j1][j2] - VU[j1][j2]);
      }
    }
  }
  
  return;
}

void elementMatrix2DP(DataTable *data, double *x1, double *y1, 
		      long long int matID,
		      std::complex<double> BB[N_EN][N_EN],
		      std::complex<double> BBeh[N_EN][N_EN]) {
  long long int i, j, l, j1, j2;
  double nx[N_NODE], ny[N_NODE];
  double N[N_NODE], Nx[N_NODE], Ny[N_NODE];
  double U[N_EDGE], V[N_EDGE], Uy[N_EDGE], Vx[N_EDGE];
  double UU[N_EDGE][N_EDGE], VV[N_EDGE][N_EDGE];
  // double UyUy[N_EDGE][N_EDGE], VxVx[N_EDGE][N_EDGE];
  // double VxUy[N_EDGE][N_EDGE], UyVx[N_EDGE][N_EDGE];
  double UV[N_EDGE][N_EDGE], VU[N_EDGE][N_EDGE];
  double UNx[N_EDGE][N_NODE], VNy[N_EDGE][N_NODE];
  // double NN[N_NODE][N_NODE], NxNx[N_NODE][N_NODE], NyNy[N_NODE][N_NODE];
  
  double b[3], c[3], ll[3], ss[3];
  double al[3];
  double W;
  double L[3];
  double Lx[3], Ly[3];
  // double k0;
  // double px, py, pz, qx, qy, qz;
  double pp, qq;
  
  // std::complex<double> cj = std::complex<double>(0, 1);
  
  long long int l_edge = 8, l_node = 6;
  
  double jj[2][2];
  double dj;
  long long int kaki;
  
  FEM *fem = &(data->fem);
  Param *par = &(data->par);
  
  switch (fem->order) {
  case 1:
    l_edge = 3;
    l_node = 3;
    break;
  case 2:
    l_edge = 6;
    l_node = 6;
    break;
  case 3:
    l_edge = 8;
    l_node = 6;
    break;
  }
  
  for (i = 0; i < N_NODE; i++)
    N[i] = Nx[i] = Ny[i] = nx[i] = ny[i] = 0.0;
  for (i = 0; i < N_EDGE; i++)
    U[i] = V[i] = Uy[i] = Vx[i] = 0.0;
  
  for (i = 0; i < N_EN; i++) {
    for (j = 0; j < N_EN; j++) {
      BB[i][j] = 0.0;
      BBeh[i][j] = 0.0;
    }
  }
  
  // k0 = 2.0 * M_PI / par->wavelength;
  kaki = checkLongEdge(x1, y1, al);
  
  create2dElementParam(x1, y1, b, c, ll, ss);
  
  
  for (l = 0; l < 7; l++) {
    setMaterialInfo2d(data, matID, &pp, &qq);
    
    setAreaCoordinate2d(l, L, &W);
    setNodalShapeFunc2d(data, L, N, nx, ny);
    
    setJacobianMatrix2d(data, x1, y1, nx, ny, jj);
    dj = jj[0][0] * jj[1][1] - jj[0][1] * jj[1][0];
    
    setAreaCoordinateDiff2d(jj, Lx, Ly);
    setNodalShapeFuncDiff2d(data, jj, nx, ny, Nx, Ny);
    
    setEdgeShapeFunc2d(data, L, Lx, Ly, jj, al, kaki, ss, U, V);
    setEdgeShapeFuncDiff2d(data, L, Lx, Ly, jj, al, kaki, ss, Uy, Vx);
    
    for (j1 = 0; j1 < l_edge; ++j1) {
      for (j2 = 0; j2 < l_edge; ++j2) {
	UU[j1][j2] = dj * (0.5 * W) * U[j1] * U[j2];
	VV[j1][j2] = dj * (0.5 * W) * V[j1] * V[j2];
	UV[j1][j2] = dj * (0.5 * W) * U[j1] * V[j2];
	VU[j1][j2] = dj * (0.5 * W) * V[j1] * U[j2];
	
      }
    }
    
    for (j1 = 0; j1 < l_edge; ++j1) {
      for (j2 = 0; j2 < l_node; ++j2) {
	VNy[j1][j2] = dj * (0.5 * W) * V[j1] * Ny[j2];
	UNx[j1][j2] = dj * (0.5 * W) * U[j1] * Nx[j2];
      }
    }
    
    
    for (j1 = 0; j1 < l_edge; ++j1) {
      for (j2 = 0; j2 < l_edge; ++j2) {
	BB[j1][j2] += pp * (UU[j1][j2] + VV[j1][j2]);
      }
    }
    for (j1 = 0; j1 < l_edge; ++j1) {
      for (j2 = 0; j2 < l_node; ++j2) {
	if (fem->order == 1) {
	  BB[j1][j2 + l_edge] += pp * (UNx[j1][j2] + VNy[j1][j2]);
	} else {
	  BB[j1][j2 + l_edge] += pp * (UNx[j1][j2] + VNy[j1][j2]);
	}
      }
    }
    
    for (j1 = 0; j1 < l_edge; ++j1) {
      for (j2 = 0; j2 < l_edge; ++j2) {
	BBeh[j1][j2] += (UV[j1][j2] - VU[j1][j2]);
      }
    }
  }
  
  return;
}

long long int checkLongEdge(double *x1, double *y1, double al[]) {
  long long int kaki;
  
  al[0] = sqrt((x1[1] - x1[0]) * (x1[1] - x1[0]) 
	       + (y1[1] - y1[0]) * (y1[1] - y1[0]));
  al[1] = sqrt((x1[2] - x1[1]) * (x1[2] - x1[1]) 
	       + (y1[2] - y1[1]) * (y1[2] - y1[1]));
  al[2] = sqrt((x1[0] - x1[2]) * (x1[0] - x1[2]) 
	       + (y1[0] - y1[2]) * (y1[0] - y1[2]));
  
  if (al[0] > al[1] && al[0] > al[2]) {
    kaki = 0;
  } else if (al[1] > al[0] && al[1] > al[2]) {
    kaki = 1;
  } else {
    kaki = 2;
  }
  
  return kaki;
}

void create2dElementParam(double *x1, double *y1, double *b, double *c,
			  double *ll, double *ss) {
  long long int i;
  
  b[0] = y1[1] - y1[2];
  b[1] = y1[2] - y1[0];
  b[2] = y1[0] - y1[1];
  c[0] = x1[2] - x1[1];
  c[1] = x1[0] - x1[2];
  c[2] = x1[1] - x1[0];
  
  ll[0] = sqrt((x1[2] - x1[1]) * (x1[2] - x1[1]) 
	       + (y1[2] - y1[1]) * (y1[2] - y1[1]));
  ll[1] = sqrt((x1[0] - x1[2]) * (x1[0] - x1[2]) 
	       + (y1[0] - y1[2]) * (y1[0] - y1[2]));
  ll[2] = sqrt((x1[1] - x1[0]) * (x1[1] - x1[0]) 
	       + (y1[1] - y1[0]) * (y1[1] - y1[0]));
  
  for (i = 0; i < 3; i++) {
    if (b[i] > ZERO || (std::abs(b[i]) < ZERO && c[i] < ZERO)) {
      ss[i] = 1.0;
    } else {
      ss[i] = -1.0;
    }
  }
  
  return;
}

void setNodalShapeFunc2d(DataTable *data, double L[], double N[], double nx[],
			 double ny[]) {
  FEM *fem = &(data->fem);
  
  if (fem->order == 1) {
    N[0] = L[0];
    nx[0] = 1.0;
    ny[0] = 0.0;
    N[1] = L[1];
    nx[1] = 0.0;
    ny[1] = 1.0;
    N[2] = L[2];
    nx[2] = -1.0;
    ny[2] = -1.0;
  } else {
    N[0] = L[0] * (2.0 * L[0] - 1.0);
    nx[0] = 4.0 * L[0] - 1.0;
    ny[0] = 0.0;
    N[1] = L[1] * (2.0 * L[1] - 1.0);
    nx[1] = 0.0;
    ny[1] = 4.0 * L[1] - 1.0;
    N[2] = L[2] * (2.0 * L[2] - 1.0);
    nx[2] = 1.0 - 4.0 * L[2];
    ny[2] = 1.0 - 4.0 * L[2];
    N[3] = 4.0 * L[0] * L[1];
    nx[3] = 4.0 * L[1];
    ny[3] = 4.0 * L[0];
    N[4] = 4.0 * L[1] * L[2];
    nx[4] = -4.0 * L[1];
    ny[4] = 4.0 * (L[2] - L[1]);
    N[5] = 4.0 * L[2] * L[0];
    nx[5] = 4.0 * (L[2] - L[0]);
    ny[5] = -4.0 * L[0];
  }
  
  return;
}

void setNodalShapeFuncDiff2d(DataTable *data, double jj[][2], double nx[],
			     double ny[], double Nx[], double Ny[]) {
  long long int i;
  long long int l_node = 6;
  double dj;
  
  FEM *fem = &(data->fem);
  
  switch (fem->order) {
  case 1:
    l_node = 3;
    break;
  case 2:
    l_node = 6;
    break;
  case 3:
    l_node = 6;
    break;
  }
  
  dj = jj[0][0] * jj[1][1] - jj[0][1] * jj[1][0];
  
  for (i = 0; i < l_node; i++) {
    Nx[i] = (jj[1][1] * nx[i] - jj[0][1] * ny[i]) / dj;
    Ny[i] = (-jj[1][0] * nx[i] + jj[0][0] * ny[i]) / dj;
  }
  
  return;
}

void setJacobianMatrix2d(DataTable *data, double *x1, double *y1, double nx[],
			 double ny[], double jj[][2]) {
  long long int i;
  long long int l_node = 6;
  
  FEM *fem = &(data->fem);
  
  switch (fem->order) {
  case 1:
    l_node = 3;
    break;
  case 2:
    l_node = 6;
    break;
  case 3:
    l_node = 6;
    break;
  }
  
  for (jj[0][0] = 0.0, i = 0; i < l_node; i++)
    jj[0][0] += nx[i] * x1[i];
  for (jj[0][1] = 0.0, i = 0; i < l_node; i++)
    jj[0][1] += nx[i] * y1[i];
  for (jj[1][0] = 0.0, i = 0; i < l_node; i++)
    jj[1][0] += ny[i] * x1[i];
  for (jj[1][1] = 0.0, i = 0; i < l_node; i++)
    jj[1][1] += ny[i] * y1[i];
  
  return;
}

void setAreaCoordinate2d(long long int i, double L[], double *W) {
  switch (i) {
  case 0:
    L[0] = 1.0 / 3.0;
    L[1] = 1.0 / 3.0;
    L[2] = 1.0 / 3.0;
    *W = 0.225;
    break;
  case 1:
    L[0] = AG1;
    L[1] = BG1;
    L[2] = BG1;
    *W = 0.13239415;
    break;
  case 2:
    L[0] = BG1;
    L[1] = AG1;
    L[2] = BG1;
    *W = 0.13239415;
    break;
  case 3:
    L[0] = BG1;
    L[1] = BG1;
    L[2] = AG1;
    *W = 0.13239415;
    break;
  case 4:
    L[0] = AG2;
    L[1] = BG2;
    L[2] = BG2;
    *W = 0.12593918;
    break;
  case 5:
    L[0] = BG2;
    L[1] = AG2;
    L[2] = BG2;
    *W = 0.12593918;
    break;
  case 6:
    L[0] = BG2;
    L[1] = BG2;
    L[2] = AG2;
    *W = 0.12593918;
    break;
  }
  
  return;
}

void setMaterialInfo2d(DataTable *data, long long int matID, double *pp, double *qq) {
  FEM *fem = &(data->fem);
  Param *par = &(data->par);
  
  if (fem->unknown == 1) {
    *pp = 1.0 / par->mr[matID];
    *qq = par->er[matID];
  } else {
    *pp = 1.0 / par->er[matID];
    *qq = par->mr[matID];
  }
  
  return;
}

void setAreaCoordinateDiff2d(double jj[][2], double Lx[], double Ly[]) {
  double dj;

  dj = jj[0][0] * jj[1][1] - jj[0][1] * jj[1][0];
  
  Lx[0] = jj[1][1] / dj;
  Lx[1] = -jj[0][1] / dj;
  Lx[2] = (-jj[1][1] + jj[0][1]) / dj;
  Ly[0] = -jj[1][0] / dj;
  Ly[1] = jj[0][0] / dj;
  Ly[2] = (jj[1][0] - jj[0][0]) / dj;
  
  return;
}

void setEdgeShapeFunc2d(DataTable *data, double L[], double Lx[], double Ly[],
			double jj[][2], double al[], long long int kaki, 
			double ss[], double U[], double V[]) {
  long long int k;
  double s12, s23, s31;
  double dj;
  FEM *fem = &(data->fem);
  
  k = kaki;
  s23 = ss[0];
  s31 = ss[1];
  s12 = ss[2];
  dj = jj[0][0] * jj[1][1] - jj[0][1] * jj[1][0];
  
  if (fem->order == 1) {
    U[0] = s12 * al[0] * (L[0] * Lx[1] - L[1] * Lx[0]);
    U[1] = s23 * al[1] * (L[1] * Lx[2] - L[2] * Lx[1]);
    U[2] = s31 * al[2] * (L[2] * Lx[0] - L[0] * Lx[2]);
    
    V[0] = s12 * al[0] * (L[0] * Ly[1] - L[1] * Ly[0]);
    V[1] = s23 * al[1] * (L[1] * Ly[2] - L[2] * Ly[1]);
    V[2] = s31 * al[2] * (L[2] * Ly[0] - L[0] * Ly[2]);
  } else {
    U[0] = al[0] * L[0] * Lx[1];
    U[1] = al[1] * L[1] * Lx[2];
    U[2] = al[2] * L[2] * Lx[0];
    U[3] = al[0] * L[1] * Lx[0];
    U[4] = al[1] * L[2] * Lx[1];
    U[5] = al[2] * L[0] * Lx[2];
#ifdef LEE
    U[6] = 4.0 * dj * (L[(1 + k) % 3] * L[(2 + k) % 3] * Lx[(0 + k) % 3])
      / al[(1 + k) % 3];
    U[7] = 4.0 * dj * (L[(2 + k) % 3] * L[(0 + k) % 3] * Lx[(1 + k) % 3])
      / al[(2 + k) % 3];
    U[8] = 4.0 * dj * (L[(0 + k) % 3] * L[(1 + k) % 3] * Lx[(2 + k) % 3])
      / al[(0 + k) % 3];
#else
    U[6] = 4.0*al[0]*(L[1]*L[2]*Lx[0] - L[2]*L[0]*Lx[1]);
    U[7] = 4.0*al[1]*(L[2]*L[0]*Lx[1] - L[0]*L[1]*Lx[2]);
    U[8] = 4.0*al[2]*(L[0]*L[1]*Lx[2] - L[1]*L[2]*Lx[0]);
#endif
    V[0] = al[0] * L[0] * Ly[1];
    V[1] = al[1] * L[1] * Ly[2];
    V[2] = al[2] * L[2] * Ly[0];
    V[3] = al[0] * L[1] * Ly[0];
    V[4] = al[1] * L[2] * Ly[1];
    V[5] = al[2] * L[0] * Ly[2];
#ifdef LEE
    V[6] = 4.0 * dj * (L[(1 + k) % 3] * L[(2 + k) % 3] * Ly[(0 + k) % 3])
      / al[(1 + k) % 3];
    V[7] = 4.0 * dj * (L[(2 + k) % 3] * L[(0 + k) % 3] * Ly[(1 + k) % 3])
      / al[(2 + k) % 3];
    V[8] = 4.0 * dj * (L[(0 + k) % 3] * L[(1 + k) % 3] * Ly[(2 + k) % 3])
      / al[(0 + k) % 3];
#else
    V[6] = 4.0*al[0]*(L[1]*L[2]*Ly[0] - L[2]*L[0]*Ly[1]);
    V[7] = 4.0*al[1]*(L[2]*L[0]*Ly[1] - L[0]*L[1]*Ly[2]);
    V[8] = 4.0*al[2]*(L[0]*L[1]*Ly[2] - L[1]*L[2]*Ly[0]);
#endif
  }
  
  return;
}

void setEdgeShapeFuncDiff2d(DataTable *data, double L[], double Lx[],
			    double Ly[], double jj[][2], double al[], 
			    long long int kaki, double ss[],
			    double Uy[], double Vx[]) {
  long long int k;
  double s12, s23, s31;
  double dj;
  FEM *fem = &(data->fem);
  
  k = kaki;
  s23 = ss[0];
  s31 = ss[1];
  s12 = ss[2];
  dj = jj[0][0] * jj[1][1] - jj[0][1] * jj[1][0];

  if (fem->order == 1) {
    Uy[0] = s12 * al[0] * (Ly[0] * Lx[1] - Ly[1] * Lx[0]);
    Uy[1] = s23 * al[1] * (Ly[1] * Lx[2] - Ly[2] * Lx[1]);
    Uy[2] = s31 * al[2] * (Ly[2] * Lx[0] - Ly[0] * Lx[2]);
    
    Vx[0] = s12 * al[0] * (Lx[0] * Ly[1] - Lx[1] * Ly[0]);
    Vx[1] = s23 * al[1] * (Lx[1] * Ly[2] - Lx[2] * Ly[1]);
    Vx[2] = s31 * al[2] * (Lx[2] * Ly[0] - Lx[0] * Ly[2]);
  } else {
    Uy[0] = al[0] * Ly[0] * Lx[1];
    Uy[1] = al[1] * Ly[1] * Lx[2];
    Uy[2] = al[2] * Ly[2] * Lx[0];
    Uy[3] = al[0] * Ly[1] * Lx[0];
    Uy[4] = al[1] * Ly[2] * Lx[1];
    Uy[5] = al[2] * Ly[0] * Lx[2];
#ifdef LEE
    Uy[6] = 4.0 * dj * ((Ly[(1 + k) % 3] * L[(2 + k) % 3] + L[(1 + k) % 3]
			 * Ly[(2 + k) % 3]) * Lx[(0 + k) % 3]) / al[(1 + k) % 3];
    Uy[7] = 4.0 * dj * ((Ly[(2 + k) % 3] * L[(0 + k) % 3] + L[(2 + k) % 3]
			 * Ly[(0 + k) % 3]) * Lx[(1 + k) % 3]) / al[(2 + k) % 3];
    Uy[8] = 4.0 * dj * ((Ly[(0 + k) % 3] * L[(1 + k) % 3] + L[(0 + k) % 3]
			 * Ly[(1 + k) % 3]) * Lx[(2 + k) % 3]) / al[(0 + k) % 3];
#else
    Uy[6] = 4.0*al[0]*((Ly[1]*L[2] + L[1]*Ly[2])*Lx[0]
		       - (Ly[2]*L[0] + L[2]*Ly[0])*Lx[1]);
    Uy[7] = 4.0*al[1]*((Ly[2]*L[0] + L[2]*Ly[0])*Lx[1]
		       - (Ly[0]*L[1] + L[0]*Ly[1])*Lx[2]);
    Uy[8] = 4.0*al[2]*((Ly[0]*L[1] + L[0]*Ly[1])*Lx[2]
		       - (Ly[1]*L[2] + L[1]*Ly[2])*Lx[0]);
#endif
    Vx[0] = al[0] * Lx[0] * Ly[1];
    Vx[1] = al[1] * Lx[1] * Ly[2];
    Vx[2] = al[2] * Lx[2] * Ly[0];
    Vx[3] = al[0] * Lx[1] * Ly[0];
    Vx[4] = al[1] * Lx[2] * Ly[1];
    Vx[5] = al[2] * Lx[0] * Ly[2];
#ifdef LEE
    Vx[6] = 4.0 * dj * ((Lx[(1 + k) % 3] * L[(2 + k) % 3] + L[(1 + k) % 3]
			 * Lx[(2 + k) % 3]) * Ly[(0 + k) % 3]) / al[(1 + k) % 3];
    Vx[7] = 4.0 * dj * ((Lx[(2 + k) % 3] * L[(0 + k) % 3] + L[(2 + k) % 3]
			 * Lx[(0 + k) % 3]) * Ly[(1 + k) % 3]) / al[(2 + k) % 3];
    Vx[8] = 4.0 * dj * ((Lx[(0 + k) % 3] * L[(1 + k) % 3] + L[(0 + k) % 3]
			 * Lx[(1 + k) % 3]) * Ly[(2 + k) % 3]) / al[(0 + k) % 3];
#else
    Vx[6] = 4.0*al[0]*((Lx[1]*L[2] + L[1]*Lx[2])*Ly[0]
		       - (Lx[2]*L[0] + L[2]*Lx[0])*Ly[1]);
    Vx[7] = 4.0*al[1]*((Lx[2]*L[0] + L[2]*Lx[0])*Ly[1]
		       - (Lx[0]*L[1] + L[0]*Lx[1])*Ly[2]);
    Vx[8] = 4.0*al[2]*((Lx[0]*L[1] + L[0]*Lx[1])*Ly[2]
		       - (Lx[1]*L[2] + L[1]*Lx[2])*Ly[0]);
#endif
  }
  
  return;
}
