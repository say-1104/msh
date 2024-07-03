#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define double_complex std::complex<double>

#include "multifrontal.h"

// #define REAL_VALUE
#define SVEA
// #define CG
#define _NO_ADAPTIVE
// #define ADAPTIVE
#define TIME_DOMAIN

#define Element Element3D

#define PI 		3.1415926535897
#define Z0      (PI*120.0)
#define	True	1
#define	False	0
#define N_PML	100

#define SLAB_MODE    0
#define PERIODIC_FILE 1
#define PLANE_WAVE   2
#define POINT_WAVE   3
#define INPUTFILE    4

#define CW_WAVE   0

#define OPTICAL  0
#define ACOUSTIC 1

#define ALLOCATION(a,c,n)  if ((a = (c *)malloc((n)*sizeof(c))) == NULL) {\
                                fprintf(stderr, "can't allocate memory\n");\
                                exit(0);\
                            }

#include "struct.h"
#include "function.h"
