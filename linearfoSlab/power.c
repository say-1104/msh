#include "fem.h"

// function for normalizing amplitude of guided mode
// according to given optical power

extern int NormalizeField(DataTable *data, double *field, double neff)
{
  int i, j, k;
  double scalingfactor;
  double power, scale;

  scalingfactor = sqrt(1.0/normalize(data, field, -1));
  fprintf(stderr, "scalingfactor = %e\n", scalingfactor);
  
  for(i = 0; i < data->fem.np; i++){
    field[i] = scalingfactor*field[i];
  }  

  if(data->par.symFlag == 2 || data->par.symFlag == 3){
    power = data->par.power/2.0;
  }else if(data->par.symFlag == 4){
    power = data->par.power/4.0;
  }

  if(data->input.modeID == 1){
    scale = 1.0e3*sqrt(2.0*(Z0)*data->par.power/neff);
  }else{
    scale = 1.0e3*sqrt(2.0*data->par.power/(neff*(Z0)));
  }
  
  for(i = 0; i < data->fem.np; i++){
    field[i] = field[i]*scale;
  }   

}

// function for calculating scaling factor

extern double normalize(DataTable *data, double *field, int flag)
{
  int i, j, k;
  int matID, ne;
  int kk[3];
  double xx[3], le;
  double_complex FF[3];
  double epsilon;
  double power, elementpower;

  power = 0.0;
  ne = data->fem.ne;

  for(i = 0; i < ne; i++){
    elementpower = 0.0;

    matID = data->fem.matID[i];

    if(matID == flag || flag == -1){

#ifdef CG
	for(j = 0; j < 3; j++){
	  kk[j] = data->fem.element[i].kk[j];
	}
#else
	for(j = 0; j < 3; j++){
	  kk[j] = data->fem.element[i][j];
	}
#endif

      for(j = 0; j < 3; j++){
	xx[j] = data->fem.xx[kk[j]];
	FF[j] = field[kk[j]];
      }
      le = (data->fem.xx[kk[1]]-data->fem.xx[kk[0]]);

      epsilon = real(data->fem.epsilon[matID]);
      
      integr(data, xx, FF, kk, le, epsilon, &elementpower);     
      power += elementpower;

    }
  }

  //  fprintf(stderr, "power = %e\n", power);

  return(power);

}

#define GZAI1  -0.90617985
#define GZAI2  -0.53846931
#define GZAI3  0.0
#define GZAI4  0.53846931
#define GZAI5  0.90617985
 
extern int integr(DataTable *data, double xx[3], double_complex FF[3], 
                  int kk[3], double le, double NN, double *elementpower)
{
  int     i,j,k;
  double  gzai, w, L1, L2, pp;
  double  N[3];
  double  Ex_square;
  double_complex Ex;

  if(data->input.modeID == 1){
    pp = 1.0;
  }else{
    pp = 1.0/(NN*NN);
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
    
    N[0] = L1*(2*L1-1);
    N[1] = L2*(2*L2-1);
    N[2] = 4*L1*L2;

    Ex = FF[0]*N[0]+FF[1]*N[1]+FF[2]*N[2];
    
    Ex_square = abs(Ex)*abs(Ex);
	 
    *elementpower += 0.5*le*w*pp*Ex_square;

  }

}
