#include "fem.h"

// function for interpolating field for output
// for TE mode, Ex, Hy, Hz
// for TM mode, Hx, Ey, Ez

extern int intpol(DataTable *data, int ii, int number, double *field)
{

  int i, j, k, ie, matID;
  int np, ne;
  double  x0 = 0.0;
  double xmin, xmax, dx;
  double phi, phi_diff, phi_z, epsilon;
  double *phix, *phiy, *phiz;
  double le, b1, b2, L1, L2, N[3], Ny[3], x[3], f[3];

  FILE *fp;
  char file[256];

  //  fprintf(stderr, "%d\n", data->par.div);

  b1 = -1.0; b2 = 1.0;

  np = data->fem.np;
  ne = data->fem.ne;
  
  ALLOCATION(data->par.x, double, data->par.div)
  ALLOCATION(phix, double, data->par.div)
  ALLOCATION(phiy, double, data->par.div)
  ALLOCATION(phiz, double, data->par.div)

  xmin = data->fem.xx[0];
  xmax = data->fem.xx[np-1];

  dx = (xmax-xmin)/(double)data->par.div;

  for(i = 0; i < data->par.div; i++){
    data->par.x[i] = xmin+0.5*dx+(double)i*dx;
  }

  // interpolate field at point x0
  for(i = 0; i < data->par.div; i++){
    x0 = data->par.x[i];
    
    // search an element including target point x0
    for(k = 0; k < ne; k++){

      for(j = 0; j < 3; j++){

#ifdef CG
	x[j] = data->fem.xx[data->fem.element[k].kk[j]];
	f[j] = field[data->fem.element[k].kk[j]];
#else
	x[j] = data->fem.xx[data->fem.element[k][j]];
	f[j] = field[data->fem.element[k][j]];
#endif

      }

      // if x0 is in element k
      if(x[0] <= x0 && x0 < x[1]){
	le = x[1]-x[0];	
	matID = data->fem.matID[k];
	epsilon = real(data->fem.epsilon[matID]);

	L1 = (x[1]-x0)/le;
	L2 = (-x[0]+x0)/le;

	N[0] = L1*(2.0*L1-1.0);
	N[1] = L2*(2.0*L2-1.0);
	N[2] = 4.0*L1*L2;

	Ny[0] = b1*(4.0*L1-1.0)/le;
	Ny[1] = b2*(4.0*L2-1.0)/le;
	Ny[2] = 4.0*(b1*L2+b2*L1)/le;

	// field at x0
	phi = N[0]*f[0]+N[1]*f[1]+N[2]*f[2];
	// differentiation of the field at x0
	phi_diff = Ny[0]*f[0]+Ny[1]*f[1]+Ny[2]*f[2];

	if(data->input.modeID == 1){
	  // TE mode: Ex, Hy, Hz
	  phix[i] = phi;
	  phiy[i] = (data->par.beta[ii][number]/(data->par.k0*Z0))*phi;
	  phiz[i] = (-1.0/(data->par.k0*Z0))*phi_diff;
	}else{
	  // TM mode: Hx, Ey, Ez
	  phix[i] = phi;
	  phiy[i] = -1.0*((data->par.beta[ii][number]*Z0)
			  /(data->par.k0*epsilon))*phi;
	  phiz[i] = (Z0/(data->par.k0*epsilon))*phi_diff;
	}

	break;

      }

    }
    
  }

  // output to file
  if(data->input.modeID == 1){
    // TE mode: Ex, Hy, Hz
    sprintf(file, "Ex-%1.3f-%d.slv\0", data->par.wl[ii], number);
    fp = fopen(file, "w");

    fprintf(fp, "# %d %15.10lf\n", data->par.div, data->par.beta[ii][number]);
    for(i = 0; i < data->par.div; i++){
      fprintf(fp, "%15.10lf %15.10lf\n", data->par.x[i], phix[i]);
    }

    fclose(fp);

    sprintf(file, "Hy-%1.3f-%d.slv\0", data->par.wl[ii], number);
    fp = fopen(file, "w");

    fprintf(fp, "# %d %15.10lf\n", data->par.div, data->par.beta[ii][number]);
    for(i = 0; i < data->par.div; i++){
      fprintf(fp, "%15.10lf %15.10lf\n", data->par.x[i], phiy[i]);
    }

    fclose(fp);

    sprintf(file, "Hz-%1.3f-%d.slv\0", data->par.wl[ii], number);
    fp = fopen(file, "w");

    fprintf(fp, "# %d %15.10lf\n", data->par.div, data->par.beta[ii][number]);
    for(i = 0; i < data->par.div; i++){
      fprintf(fp, "%15.10lf %15.10lf\n", data->par.x[i], phiz[i]);
    }

    fclose(fp);

  }else{
    // TM mode: Hx, Ey, Ez
    sprintf(file, "Hx-%1.3f-%d.slv\0", data->par.wl[ii], number);
    fp = fopen(file, "w");

    fprintf(fp, "# %d %15.10lf\n", data->par.div, data->par.beta[ii][number]);
    for(i = 0; i < data->par.div; i++){
      fprintf(fp, "%15.10lf %15.10lf\n", data->par.x[i], phix[i]);
    }

    fclose(fp);

    sprintf(file, "Ey-%1.3f-%d.slv\0", data->par.wl[ii], number);
    fp = fopen(file, "w");

    fprintf(fp, "# %d %15.10lf\n", data->par.div, data->par.beta[ii][number]);
    for(i = 0; i < data->par.div; i++){
      fprintf(fp, "%15.10lf %15.10lf\n", data->par.x[i], phiy[i]);
    }

    fclose(fp);

    sprintf(file, "Ez-%1.3f-%d.slv\0", data->par.wl[ii], number);
    fp = fopen(file, "w");

    fprintf(fp, "# %d %15.10lf\n", data->par.div, data->par.beta[ii][number]);
    for(i = 0; i < data->par.div; i++){
      fprintf(fp, "%15.10lf %15.10lf\n", data->par.x[i], phiz[i]);
    }

    fclose(fp);

  }


  free(data->par.x);
  free(phix);
  free(phiy);
  free(phiz);

}
