#ifndef MAIN_HPP
#define MAIN_HPP

#include </usr/include/eigen3/Eigen/Dense>
//#include </home/okazaki/Solver/Eigen/Dense>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>
#include "LI.hpp"
#include "function.hpp"
#include "struct.hpp"

constexpr std::complex<double> cj(0.0, 1.0);

#define ALLOCATION(a,c,n) \
		if ((a = (c *)malloc((n)*sizeof(c))) == NULL) {\
				fprintf(stderr, "can't allocate memory\n");\
				exit(EXIT_FAILURE);\
		}

#define SI -1
#define CPCM 0
#define APCM 1



#endif