#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <vector>
#include <algorithm>

#define double_complex std::complex<double>
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

// #define SYMMETRY
// definition below indicates that yz-symmetry would be 0x3 (= y_sym & z_sym)
#define Y_SYMMETRY 0x1
#define Z_SYMMETRY 0x2
#define YeqZ_SYMMETRY 0x4     // symmetric about y=z


#define XY_COORDINATE

#define SINGLE_DIRECTION 1
#define BOTH_DIRECTION   0

typedef enum _SOLVER {
    SOLVER_NONE,
    INPUTEIGEN,
    POINTSOURCE,
    SHEETSOURCE
} SOLVER;

typedef enum _DEVICETYPE {
  DEVICE_NONE,
  MODE_DIVIDER,
  MODE_MUX,
  POWER_SPLITTER,
  WAVEGUIDE_CROSSING,
} DEVICETYPE;

#define OFF 0
#define ON  1

#define ZERO -1.0e-6

#define PI		    3.141592653589793
#define Z0		    (120.0*PI)
#define LENGTH(a,b)	sqrt(pow(a.x-b.x, 2.0)+pow(a.y-b.y, 2.0))

#define MAX_COMPONENT  24

#include "struct.h"
#include "function.h"

#define ALLOCATION(a,c,n)  if ((a = (c *)malloc((n)*sizeof(c))) == NULL) {\
                                fprintf(stderr, "can't allocate memory\n");\
                                exit(0);\
                            }

