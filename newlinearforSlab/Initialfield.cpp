
#include "fem.hpp"

// function for setting initial field guess 
// for inverse iterative method

extern int Initialfield(DataTable *data)
{

  int	  i, j, k, ie, matID;
  double  x0 = 0.0, y0 = 0.5, xx, yy;

  for(i = 0; i < data->fem.np; i++){
    data->fem.field[i] = 1.0;
  }

  /*	
    // Gaussian
    for(i = 0; i < data->fem.np; i++){
    xx = data->fem.xx[i];
    data->fem.field[i] = exp(-(xx*xx)/(x0*x0));
    }
  */	

}
