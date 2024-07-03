#include "fem.hpp"

#define GZAI1  -0.90617985
#define GZAI2  -0.53846931
#define GZAI3  0.0
#define GZAI4  0.53846931
#define GZAI5  0.90617985
 
// function for making finite-element matrix

extern int matrix(DataTable *data, double *xx, double_complex *f,
		  double_complex sk[3][3], double_complex sm[3][3], int matID, 
		  double_complex *Eps, int element)
{
  int i, j, k;
  double xl, yl, zl, w;
  double u, v, r0, r1, r2;
  double tmp1, tmp2;
  double A11, A12, A13;
  double A21, A22, A23;
  double A31, A32, A33;
  double a1, b1, a2, b2, le, gzai, rr, ll;
  double L1, L2;
  double f1a[3], f2a[3], nn[3], nr[3], ny[3];
  double j11, j12, j21, j22, dj, dj_inv;
  double f1[3][3], f2[3][3], f3[3][3], f4[3][3], f5[3][3], f6[3][3];
  double g1[3][3], g2[3][3], g3[3][3], g4[3][3], g5[3][3], g6[3][3];
  double_complex pp, qq;
  double Eps_re, Eps_im, ff2, kerr, sat, Kerr, n_L;
  double_complex  ff, fx, fy, EPS, eps, Ey, Ez;

  /* PML */
  double s_max = 1.0/(2.0*data->par.k0)*log(1.0/data->pml.tanD);
  double_complex s0, sx, sy, sz;
  double_complex Sx, Sy, Sz;
  double dd, rho;
  double mm = data->pml.m;
  double xg, yg, zg, zmax, zmin;
  double_complex cj(0.0, 1.0);

  if( data->input.modeID == 1 ){
    pp = 1.0;
    qq = data->fem.epsilon[matID];
  } else {
    pp = data->fem.epsilonInv[matID];
    qq = 1.0;
  }

  for(i = 0; i < 3; i++){
    for(j = 0; j < 3; j++){
      g1[j][i] = 0.0;
      g2[j][i] = 0.0;
      g3[j][i] = 0.0;
      g4[j][i] = 0.0;
      g5[j][i] = 0.0;
      g6[j][i] = 0.0;
      sk[j][i] = 0.0;
      sm[j][i] = 0.0;
    }
  }
	
  le = xx[1]-xx[0];

  if(le < 0.0){
    printf("stop in matrix\n");
    printf("element length became negative!\n");
    printf("element number = %d, value = %lf\n", element, dj);
    exit(-1);
  }

  for(k = 0; k < 5; k++){
    switch(k){
    case 0:   
      gzai = GZAI1; w = 0.23692689; break;
    case 1:   
      gzai = GZAI2; w = 0.47862867; break;
    case 2:   
      gzai = GZAI3; w = 0.56888889; break;
    case 3:   
      gzai = GZAI4; w = 0.47862867; break;
    case 4:   
      gzai = GZAI5; w = 0.23692689; break;
    }

    L1 = (1-gzai)/2;  L2 = (1+gzai)/2;
    
    nn[0] = L1*(2*L1-1);
    nn[1] = L2*(2*L2-1);
    nn[2] = 4*L1*L2;

    ny[0] = -2.0*L1+0.5;
    ny[1] = 2.0*L2-0.5;
    ny[2] = 2.0*(L1-L2);

    /*
      sx = sy = sz = 1.0;
      for (i = 0; i < data->pml.nx; i++) {
      if ((data->pml.x_data[i].x0-xg)*(data->pml.x_data[i].x1-xg) < 0.0 &&
      (data->pml.x_data[i].y0-yg)*(data->pml.x_data[i].y1-yg) < 0.0) {
      dd = abs(data->pml.x_data[i].dx);
      rho = abs(xg-data->pml.x_data[i].x0);
      s0 = 1.0-cj*s_max/dd*pow(rho/dd, mm)*(1.0+mm);
      sx = s0;
      }
      }
      for (i = 0; i < data->pml.ny; i++) {
      if ((data->pml.y_data[i].x0-xg)*(data->pml.y_data[i].x1-xg) < 0.0 &&
      (data->pml.y_data[i].y0-yg)*(data->pml.y_data[i].y1-yg) < 0.0) {
      dd = abs(data->pml.y_data[i].dy);
      rho = abs(yg-data->pml.y_data[i].y0);
      s0 = 1.0-cj*s_max/dd*pow(rho/dd, mm)*(1.0+mm);
      sy = s0;
      }
      }
      Sx = sy*sz/sx;
      Sy = sz*sx/sy;
      Sz = sx*sy/sz;
    */

    for(i = 0; i < 3; i++){
      for(j = 0; j < 3; j++) f1[i][j] = 0.0;
      for(j = 0; j < 3; j++) f2[i][j] = 0.0;
      for(j = 0; j < 3; j++) f3[i][j] = 0.0;
      for(j = 0; j < 3; j++) f4[i][j] = 0.0;
    }
    
    for(i = 0; i < 3; i++){
      for(j = 0; j < 3; j++) f1[i][j] = 0.0;
      for(j = 0; j < 3; j++) f1[i][j] += nn[j]*nn[i];
    }
    for(i = 0; i < 3; i++){
      for(j = 0; j < 3; j++) f2[i][j] = 0.0;
      for(j = 0; j < 3; j++) f2[i][j] += ny[j]*ny[i];
    }
    
    for(i = 0; i < 3; i++){
      for(j = 0; j < 3; j++) g1[i][j] = 0.5*le*w*f1[i][j];
		}
    for(i = 0; i < 3; i++){
      for(j = 0; j < 3; j++) g2[i][j] = 2.0*w*f2[i][j]/le;
      
    }
    
    for(j = 0; j < 3; j++){
      for(i = 0; i < 3; i++){
	sk[j][i] += (-pp*g2[j][i]
		     +data->par.k02*qq*g1[j][i]);
	
	sm[j][i] += pp*g1[j][i];
      }
    }
  }

}

#ifndef CG

// function for rearranging element of grobal matrix for bhlu1c.c

extern int reArrangment(DataTable *data)
{
  int	i, j;
  int	it1, it2;
  int	nbw, nbw2;
  
  nbw = data->fem.nbw;
  nbw2 = 2*nbw+1;

  for(i = 0; i < nbw; i++){
    for(j = 0; j < nbw+1+i; j++){
      it1 = j+nbw2*i;
      it2 = (j+nbw-i)+nbw2*i;
      data->fem.Pplus[it1] = data->fem.Pplus[it2];
    }
    for(j = nbw+1+i; j < nbw2; j++){
      it1 = j+nbw2*i;
      data->fem.Pplus[it1] = 0.0;
    }
  }

}
#endif

// function for creating grobal finite-element matrix

extern int FEMmatrix(DataTable *data)
{
  int		i, j, k, ie, matID, tmp;
  double	x[3];
  double_complex sk[3][3], sm[3][3];
  double	mdk0n0, hdz, mk02n02, qk02n02, dz;
  double_complex	c1, c2, f[3];
  int		ii, jj, it;
  double	eps;
  int		ier;
  int		count, itmp;
  double_complex cj(0.0, 1.0), Eps[5];
  double optVelocity = 0.299792458;

#ifdef CG
  CGtable	*cg_table;
#endif

#ifdef CG
  cg_table = data->fem.cg_table;
  for (i = 0; i < data->fem.np; i++) {
    for (j = 0; j < cg_table[i].n; j++) {
      data->fem.aa[i][j] = 0.0;
      data->fem.bb[i][j] = 0.0;
      data->fem.Pplus[i][j] = 0.0;
      data->fem.Pminus[i][j] = 0.0;
    }
  }
#else
  ie = (data->fem.nbw*2+1)*data->fem.np;
  
  for(i = 0; i < ie; i++){
    data->fem.aa[i] = 0.0;
    data->fem.bb[i] = 0.0;
  }

  ie = (data->fem.nbw+1)*data->fem.np;

  for(i = 0; i < ie; i++){
    data->fem.aa_half[i] = 0.0;
    data->fem.bb_half[i] = 0.0;
  }
#endif

  for(k = 0; k < data->fem.ne; k++){
    for(i = 0; i < 3; i++){
#ifdef CG
      x[i] = data->fem.xx[data->fem.element[k].kk[i]];
#else
      x[i] = data->fem.xx[data->fem.element[k][i]];
#endif
      f[i] = 0.0;
    }
    
    matID = data->fem.matID[k];
    
    matrix(data, x, f, sk, sm, matID, Eps, k);
    
    for(i = 0; i < 3; i++){
      for(j = 0; j < 3; j++){
#ifdef CG
	ii = data->fem.element[k].kk[i];
	jj = data->fem.element[k].kk[j];
	it = columnInCGmatrix(ii, jj, cg_table, data->fem.np);
	data->fem.aa[ii][it] += real(sm[i][j]);
	data->fem.bb[ii][it] += real(sk[i][j]);
#else
	ii = data->fem.element[k][i];
	jj = data->fem.element[k][j];
	it = (jj-ii+data->fem.nbw)+(data->fem.nbw*2+1)*ii;
	data->fem.aa[it] += real(sm[i][j]);
	data->fem.bb[it] += real(sk[i][j]);
	if (ii <= jj && jj < data->fem.np) {
	  it = (data->fem.nbw+1)*ii+(jj-ii);
	  data->fem.aa_half[it] += real(sm[i][j]);
	  data->fem.bb_half[it] += real(sk[i][j]);
	}
#endif
      }
    }
  }

}
