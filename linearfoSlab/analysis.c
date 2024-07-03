#include "fem.h"

extern int analysis(DataTable *data, int argc, char **argv)
{
  int i, j, k, ii, nf, count = 0, flag, tmp;
  int matID, kk;
  double xc_max, xc_temp, delta, wl;
  double p_total, p_core, total, temp;
  double_complex eps;
  char string[256];
  FILE *fp;

  fp = fopen("NEFF", "w");
  fclose(fp);
  fp = fopen("Beta", "w");
  fclose(fp);
  fp = fopen("Gamma", "w");
  fclose(fp);
  
  // check command line options
  nf = checkArgument(data, argc, argv);

  // input data
  strcpy(data->file, argv[nf]);
  inputDataMsh(data, argv[nf]);
	
  data->par.numoflamda 
    = (int)((data->par.lamda2-data->par.lamda1)/data->par.step);

  if(data->par.numoflamda == 0) data->par.numoflamda = 1;

  // mesh generation
  MeshGeneration(data);
	
  ALLOCATION(data->par.wl, double, data->par.numoflamda);
  ALLOCATION(data->par.neff, double *, data->par.numoflamda);
  ALLOCATION(data->par.beta, double *, data->par.numoflamda);
  ALLOCATION(data->par.field, double **, data->par.numoflamda);
  ALLOCATION(data->par.Gamma, double **, data->par.numoflamda);

  for(i = 0; i < data->par.numoflamda; i++){
    ALLOCATION(data->par.neff[i], double, data->input.number);
    ALLOCATION(data->par.beta[i], double, data->input.number);
    ALLOCATION(data->par.field[i], double *, data->input.number);
    ALLOCATION(data->par.Gamma[i], double *, data->input.number);
  }

  for(i = 0; i < data->par.numoflamda; i++){
    for(j = 0; j < data->input.number; j++){
      ALLOCATION(data->par.field[i][j], double, data->fem.np);
      ALLOCATION(data->par.Gamma[i][j], double, data->par.numOflayer);
    }
  }

  // set initial field guess
  Initialfield(data);

  // wavelength sweep
  for(ii = 0; ii < data->par.numoflamda; ii++){

    wl = data->par.lamda1+(double)ii*data->par.step;
    data->par.wl[ii] = wl;
    data->par.wavelength = wl;
    data->par.k0 = 2.0*PI/data->par.wavelength;
    data->par.k02 = data->par.k0*data->par.k0;	
    fprintf(stderr, "current wavelength = %lf, k0 = %lf\n", wl, data->par.k0); 

    // set refractive index if specified
    //	  material(data);

    fprintf(stderr, "current core index  = %lf\n", 
	    sqrt(real(data->fem.epsilon[0]))); 
    fprintf(stderr, "current cladding index  = %lf\n", 
	    sqrt(real(data->fem.epsilon[1]))); 

    // create finite-element matrix
    FEMmatrix(data);

    data->input.c_beta = data->par.k0*data->input.c_neff;
    fprintf(stderr, "beta = %lf\n", data->input.c_beta); 

    if(data->input.eigv == 0){
      // eigenmode calculation by inverse iterative method
      flag = 0;
      for(j = 0; j < 100; j++){
	
	propagation(data);
	
	for(i = 0; i < data->fem.np; i++){
	  xc_temp = abs(data->fem.field[i]);
	  if(i == 0 || xc_max < xc_temp) xc_max = xc_temp;
	}
	  
	for(i = 0; i < data->fem.np; i++){
	  data->fem.field[i] /= xc_max;
	}
	
	// calculate effective index
	calcNeff(data);
	
	if(flag != 0){
	  delta = fabs((data->input.c_neff-data->input.p_neff)
		       /data->input.c_neff);
	  if (delta < 1.0e-9) break;
	}
	flag = 1;
      }

      data->par.neff[ii][0] = data->par.n0;
      data->par.beta[ii][0] = data->par.n0*data->par.k0;

      for(i = 0; i < data->fem.np; i++){
	data->par.field[ii][0][i] = abs(data->fem.field[i]);
      }

      fp = fopen("NEFF", "a+");
      fprintf(fp, "%lf %15.10lf\n", wl, data->par.n0);
      fclose(fp);
      fp = fopen("Beta", "a+");
      fprintf(fp, "%lf %15.10lf\n", wl, data->par.beta);
      fclose(fp);

    }else{
#ifndef CG
      // eigenmode calculation by eigv4p
      eigenModeCalculation(data, ii);

      for(i = 0; i < data->input.number; i++){      
	fp = fopen("NEFF", "a+");
	fprintf(fp, "%lf %d %15.10lf\n", wl, i, data->par.neff[ii][i]);
	fclose(fp);
	fp = fopen("Beta", "a+");
	fprintf(fp, "%lf %d %15.10lf\n", wl, i, data->par.beta[ii][i]);
	fclose(fp);
      }

#endif
    }

    /*
    fp = fopen("InitialField", "w");
    for(i = 0; i < data->fem.np; i++){
      fprintf(fp, "%lf %lf\n", (data->par.field[0][0][i]), 0.0);
    }
    fprintf(fp, "%lf\n", data->par.beta[0][0]);  
    fclose(fp);
    */

    // calc confinement factor
    fp = fopen("Gamma", "a+");
    fprintf(fp, "Wavelength =  %15.10lf\n", wl);

    for(i = 0; i < data->input.number; i++){ 

      // scaling amplitude according to the input power
      NormalizeField(data, data->par.field[ii][i], data->par.neff[ii][i]);

      total = 0.0;
      p_total = normalize(data, data->par.field[ii][i], -1);

      for(j = 0; j < data->par.numOflayer; j++){
	p_core = normalize(data, data->par.field[ii][i], j);
	data->par.Gamma[ii][i][j] = p_core/p_total;
	total += data->par.Gamma[ii][i][j];
	fprintf(fp, "%d %15.10lf\n", j, data->par.Gamma[ii][i][j]);
      }

      fprintf(fp, "total = %15.10lf\n", total);
      
    }
    
    fclose(fp);

    // field interpolation
    for(i = 0; i < data->input.number; i++){ 
      intpol(data, ii, i, data->par.field[ii][i]);
    }

  }

  // output refractive index distributions

  fp = fopen("Index.data", "w");

  temp = 0.0;
  for(i = 0; i < data->par.numOflayer; i++){
    fprintf(fp, "%15.10lf %15.10lf\n", temp, sqrt(real(data->fem.epsilon[i])));
    temp += data->par.thickness[i];
    fprintf(fp, "%15.10lf %15.10lf\n", temp, sqrt(real(data->fem.epsilon[i])));
  }

  fclose(fp);

  // output InitialFiled for BPM

  /*
  fp = fopen("InitialField", "w");
  */
  /*
  for(i = 0; i < data->fem.np; i++){
    fprintf(fp, "%15.10lf %15.10lf\n", (data->par.field[0][0][i]), 0.0);
  }
  fprintf(fp, "%15.10lf\n", data->par.beta[0][0]);
  */
  /*
  for(i = 0; i < data->fem.np; i++){
    fprintf(fp, "%15.10lf %15.10lf\n", (data->par.field[0][0][i]), 
	    (data->par.field[0][1][i]));
  }
  fprintf(fp, "%15.10lf\n", data->par.beta[0][0]);

  fclose(fp);
  */

  // confinement factor for matID material

  temp = 0.0;
  matID = 6;
  for(i = 0; i < data->par.numOflayer; i++){
    if(data->fem.epsilon[i] == data->fem.epsilon[matID]){
      temp += data->par.Gamma[0][0][i];
    }
  }
  fprintf(stderr, "Confinmentfactor [%d] = %15.10lf\n", matID, temp);

  // calculate chromatic dispersion
  if(data->input.eigv == 0){
    calcDispersion(data);
  }

  fp = fopen("NEFF-plot", "w");
  for(i = 0; i < data->input.number; i++){    
    for(ii = 0; ii < data->par.numoflamda; ii++){

      wl = data->par.lamda1+(double)ii*data->par.step;
      fprintf(fp, "%lf %15.10lf\n", wl, data->par.neff[ii][i]);
    }

    fprintf(fp, "\n");

  }

  fclose(fp);

}
