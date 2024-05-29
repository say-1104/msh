#ifndef MAIN_HPP
#define MAIN_HPP

#include </home/okazaki/Solver/Eigen/Dense>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>

constexpr std::complex<double> cj(0.0, 1.0);
#ifndef PI
#define PI 3.1415926535897
#endif
#define Z0 376.73036

#define ALLOCATION(a,c,n)  if ((a = (c *)malloc((n)*sizeof(c))) == NULL) {\
								fprintf(stderr, "can't allocate memory\n");\
								exit(0);\
							}

#endif