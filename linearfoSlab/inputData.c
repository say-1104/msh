#include "fem.h"

// function for checking command line options

extern int checkArgument(DataTable *data, int argc, char **argv)
{
  
  int i, nf;
  
  data->problem = OPTICAL;
  data->input.eigv = 0;
  data->input.number = 1;
  data->par.lamda1 = 1.55;
  data->par.lamda2 = 1.55;
  data->par.step = 1.0;
  data->par.x0 = 0.0;
  
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      if (argv[i][1] == 'a') {
	data->input.c_neff = atof(&(argv[i][2]));
      }else if (argv[i][1] == 'e') {
	data->input.eigv = 1;
      }else if (argv[i][1] == 'n') {
	sscanf(&(argv[i][2]), "%d", &data->input.number);
      }else if (argv[i][1] == 'l') {
	sscanf(&(argv[i][2]), "%lf,%lf,%lf", &data->par.step, 
	       &data->par.lamda1, &data->par.lamda2);
      }else if (argv[i][1] == 'o' && argv[i][2] == 'f' && argv[i][3] == 'f') {
	sscanf(&(argv[i][4]), "%lf", &data->par.x0);
      }
    } else {
      nf = i;
    }
  }

  fprintf(stderr, "desired number of modes = %d\n", data->input.number);
  fprintf(stderr, "coordinate offset = %lf\n", data->par.x0);

  return nf;

}

// function for inputting data

extern int inputDataMsh(DataTable *data, char *fileName)
{
  FILE	*fp;
  char	string[256], f_file[256];
  int		i, j, k, matNum, matID;
  double	check, re_tmp, im_tmp, re, im;
  double_complex cj(0.0, 1.0);
  double	xx, yy;
  
  sprintf(string, "%s.msh", fileName);
  if( (fp = fopen(string, "r")) == NULL ){
    fprintf(stderr, "can't open input file ( %s ).\n", string);
    exit(0);
  }
  
  fscanf(fp, "%*s %*s");
  // 1st initial wavelength
  // 2nd eng wavelength
  // 3rd wavelength step
  fscanf(fp, "%lf %lf %lf", &(data->par.lamda1), &(data->par.lamda2), 
	 &data->par.step);

  data->par.wavelength = data->par.lamda1;
  data->par.k0 = 2.0*PI/data->par.wavelength;
  data->par.k02 = data->par.k0*data->par.k0;

  // input power
  fscanf(fp, "%*s = %lf", &data->par.power);

  // modeID, modeNo, and symFlag
  fscanf(fp, "%*s = %d", &(data->input.modeID));
  fscanf(fp, "%*s = %d", &(data->input.modeNo));
  fscanf(fp, "%*s = %d", &(data->par.symFlag));	

  // the number of material
  fscanf(fp, "%*s %*s %*s %*s = %d", &data->par.numOflayer);
  
  fscanf(fp, "%*s %*s %*s %*s");
  // 1st real part of refractive index
  // 2nd imaginary part of refractive index
  // 3rd Thickness of each layer in micro meter
  // 4th the number of division 

  ALLOCATION(data->par.division, int, data->par.numOflayer)
  ALLOCATION(data->par.thickness, double, data->par.numOflayer)

  for(i = 0; i < data->par.numOflayer; i++){
    fscanf(fp, "%lf %lf %lf %d", &re_tmp, &im_tmp, 
	   &(data->par.thickness[i]), &(data->par.division[i]));
    data->fem.epsilon[i] = (re_tmp+cj*im_tmp)*(re_tmp+cj*im_tmp);
    data->fem.epsilonInv[i] = 1.0/data->fem.epsilon[i];
 
    fprintf(stderr, "%lf, %d\n", data->par.thickness[i], 
	    data->par.division[i]);

  }
  
  // the number of division for output field
  fscanf(fp, "%*s %*s %*s");
  fscanf(fp, "%*s = %d", &(data->par.div));

  /*  
  // initial coordinate of x
  fscanf(fp, "%*s %*s");
  fscanf(fp, "%*s = %lf", &(data->par.x0));
  */

  fprintf(stderr, "start wavelength = %lf\n", data->par.lamda1);
  fprintf(stderr, "end wavelength = %lf\n", data->par.lamda2);
  fprintf(stderr, "wavelength step = %lf\n", data->par.step);

  fprintf(stderr, "Input power = %lf [W/m]\n", data->par.power);

  fprintf(stderr, "division for output field = %d\n", data->par.div);

}

// function for generating mesh

extern int MeshGeneration(DataTable *data)
{

  FILE	*fp;
  char	string[256], f_file[256];
  int		i, j, k, matNum, matID;
  double	temp, thick, offset = 0.0;
  double_complex cj(0.0, 1.0);
  double	xx, yy, deltax;
  
  offset = data->par.x0;

  data->fem.np = 1;
  for(i = 0; i < data->par.numOflayer; i++){
    data->fem.np += 2*data->par.division[i];
  }
  for(i = 0; i < data->par.numOflayer; i++){
    data->fem.ne += data->par.division[i];
  }

  fprintf(stderr, "np = %d, ne = %d\n", data->fem.np, data->fem.ne);
  
  ALLOCATION(data->fem.xx, double, data->fem.np)
  ALLOCATION(data->fem.matID, int, data->fem.ne)

#ifdef CG
  data->fem.element = (Element *)malloc(data->fem.ne*sizeof(Element));
#else
  ALLOCATION(data->fem.element, int *, data->fem.ne)
  for(i = 0; i < data->fem.ne; i++){
    ALLOCATION(data->fem.element[i], int, 3)
  }
#endif

  k = 0;
  temp = 0.0;
  for(i = 0; i < data->par.numOflayer; i++){
    deltax = data->par.thickness[i]/(2.0*(double)data->par.division[i]);

    for(j = 0; j < data->par.division[i]*2; j++){
      data->fem.xx[k] = offset+deltax*(double)j+temp;
      k++;
      //      fprintf(stderr, "%d, %d, xx = %lf\n", k-1, j, data->fem.xx[k-1]);
    }
    temp += data->par.thickness[i];
    data->fem.xx[k] = temp+offset;
  }
  //  fprintf(stderr, "%d, %d, xx = %lf\n", k, j, data->fem.xx[k]);

  for(i = 0; i < data->fem.ne; i++){
#ifdef CG
    data->fem.element[i].kk[0] = 2*i;
    data->fem.element[i].kk[1] = 2*i+2;
    data->fem.element[i].kk[2] = 2*i+1;
#else
    data->fem.element[i][0] = 2*i;
    data->fem.element[i][1] = 2*i+2;
    data->fem.element[i][2] = 2*i+1;

#endif
  }
	
  k = 0;
  thick = data->par.thickness[0]+offset;
  for(i = 0; i < data->fem.ne; i++){
#ifdef CG
    if(data->fem.xx[data->fem.element[i].kk[2]] < thick){
#else
    if(data->fem.xx[data->fem.element[i][2]] < thick){
#endif
      data->fem.matID[i] = k;
    }else{
      data->fem.matID[i] = k+1;
      k++;
      thick += data->par.thickness[k];
    }
  }
    
  data->fem.nbw = checkBandWidth(data);
  allocWorkSpace(data);
}

extern int maxValue(int *a, int n)
{
  int		i, imax;
  
  imax = a[0];
  for(i = 1; i < n; i++) if( a[i] > imax ) imax = a[i];
  
  return(imax);
}

extern int minValue(int *a, int n)
{
  int		i, imin;
  
  imin = a[0];
  for(i = 1; i < n; i++) if( a[i] < imin ) imin = a[i];
  
  return(imin);
}

extern int checkBandWidth(DataTable *data)
{
  int		i, imax, itmp;

  imax = 0;
  for(i = 0; i < data->fem.ne; i++){
#ifdef CG
    itmp = maxValue(data->fem.element[i].kk, 3)
      -minValue(data->fem.element[i].kk, 3);
#else
    itmp = maxValue(data->fem.element[i], 3)
      -minValue(data->fem.element[i], 3);
#endif
    if( itmp > imax ) imax = itmp;
  }
  fprintf(stderr, "band width = %d\n", imax);
  
  return(imax);
}

// fucntion for allocating memory

extern int allocWorkSpace(DataTable *data)
{
  int		i, j;
  int		np, ne, nr, nbw, nbw2;
  double	ftmp;
#ifdef CG
  CGtable	*cg_table;
#endif
  
  np = data->fem.np;
  data->fem.nr = nr = np;
  ne = data->fem.ne;
  nbw = data->fem.nbw;
  nbw2 = 2*(data->fem.nbw)+1;
  fprintf(stderr, "nbw = %d\n", data->fem.nbw);
  
#ifdef CG
  data->fem.cg_table = cg_table
    = (CGtable *)malloc(data->fem.np*sizeof(CGtable));
  makeCGmatrixTable(data->fem.element, data->fem.ne, data->fem.np,
		    cg_table, 3);
  ALLOCATION(data->fem.aa, double_complex *, data->fem.np)
  ALLOCATION(data->fem.bb, double_complex *, data->fem.np)
  ALLOCATION(data->fem.Pplus, double_complex *, data->fem.np)
  ALLOCATION(data->fem.Pminus, double_complex *, data->fem.np)
    
  for (i = 0; i < data->fem.np; i++) {
    ALLOCATION(data->fem.aa[i], double_complex, cg_table[i].n)
    ALLOCATION(data->fem.bb[i], double_complex, cg_table[i].n)
    ALLOCATION(data->fem.Pplus[i], double_complex, cg_table[i].n)
    ALLOCATION(data->fem.Pminus[i], double_complex, cg_table[i].n)
  }
#else
  ALLOCATION(data->fem.aa, double, np*nbw2)
  ALLOCATION(data->fem.bb, double, np*nbw2)
  ALLOCATION(data->fem.aa_half, double, np*(nbw+1))
  ALLOCATION(data->fem.bb_half, double, np*(nbw+1))

  ALLOCATION(data->fem.vv, double, np*nbw2)
  ALLOCATION(data->fem.Pplus, double_complex, np*nbw2)
  ALLOCATION(data->fem.Pminus, double_complex, np*nbw2)
  ALLOCATION(data->fem.cac, double_complex, np*nbw)
  ALLOCATION(data->fem.cwk, double_complex, np)
  ALLOCATION(data->fem.iwp, int, np)
#endif

  ALLOCATION(data->fem.field, double_complex, np)
  ALLOCATION(data->fem.fieldn, double_complex, np)

  ALLOCATION(data->fem.amat, double_complex, np)
  ALLOCATION(data->fem.bmat, double_complex, np)

}
