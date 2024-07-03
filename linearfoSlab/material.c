#include "fem.h"

#define  PCF

// function for setting wavelength depenedent refractive index

// Delta : relative index difference
// return core index as out[0] and cladding index as out[1]

extern double material(DataTable *data, double Delta, double out[2])
{
  int i;
  double wl;
  double rn, n_core, n_clad;

#ifdef PCF
/*
ピュアシリカ屈折率 出典不明
*/

  static double a[3] = {0.696750, 0.408218, 0.890815};
  static double l[3] = {0.069066, 0.115662, 9.900559};

  wl = data->par.wavelength;
  rn = 1.0;

  /*出典不明のシリカ屈折率の場合*/
  for (i = 0; i < 3; i++) {
    rn += a[i]*wl*wl/(wl*wl-l[i]*l[i]);
  }

  n_clad = rn;
  n_core = n_clad/(1.0-2.0*Delta);
  
#else


  // GeO2 6.3mol% ドープのシリカ屈折率 from 光導波路の基礎 by 岡本勝就

  static double a_Ge[3] = {0.7083952, 0.4203993, 0.8663412};
  static double b_Ge[3] = {7.290464e-3, 1.050294e-2, 97.93428};


  //  ピュアシリカ屈折率 from 光導波路の基礎 by 岡本勝就

  static double a[3] = {0.6965325, 0.4083099, 0.8968766};
  static double b[3] = {4.368309e-3, 1.394999e-2, 97.93399};
  
  wl = data->par.wavelength;
  n_core = n_clad = 1.0;
  
  /*光導波路の基礎のGeO2ドープシリカ屈折率の場合*/

  for (i = 0; i < 3; i++) {
    n_core += a_Ge[i]*wl*wl/(wl*wl-b_Ge[i]);
  }

  data->fem.epsilon[0] = n_core;
  data->fem.epsilonInv[0] = 1.0/(n_core);
	
  /*光導波路の基礎のシリカ屈折率の場合*/

  for (i = 0; i < 3; i++) {
    n_clad += a[i]*wl*wl/(wl*wl-b[i]);
  }

  data->fem.epsilon[1] = n_clad;
  data->fem.epsilonInv[1] = 1.0/(n_clad);

#endif

  out[0] = n_core;
  out[1] = n_clad;

  fprintf(stderr, "core index = %15.10lf\n", n_core);
  fprintf(stderr, "cladding index = %15.10lf\n", n_clad);

}







