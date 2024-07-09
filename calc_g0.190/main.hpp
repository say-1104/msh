#ifndef MAIN_HPP
#define MAIN_HPP

#include </usr/include/eigen3/Eigen/Dense>
//#include </home/okazaki/Solver/Eigen/Dense>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <array>
#include <vector>
#include <string>
#include <tuple>

constexpr std::complex<double> cj(0.0, 1.0);

#define ALLOCATION(a,c,n) \
		if ((a = (c *)malloc((n)*sizeof(c))) == NULL) {\
				fprintf(stderr, "can't allocate memory\n");\
				exit(EXIT_FAILURE);\
		}

#define SI -1
#define CPCM 0
#define APCM 1

typedef struct _Param {
    double wl;  //波長
    int N_dset; //データセットの個数
    double dz;  //伝搬方向の分割数
    double Leff;    //cPCMの長さ
    double wst;     //始点のwidth
    int N_taper;    //テーパの接続個数
    std::vector<std::pair<double, double>> taper;   //テーパの構造パラメータ
    std::vector<std::tuple<double, double, int>> ZtoW;
} Param;

typedef struct _Dataset {
    int N;  //widthの個数
    LI<double, double> beta_even;
    LI<double, double> beta_odd;
    LI<double, double> beta_1;
    LI<double, double> beta_2;
} Dataset;

typedef struct _Flag {
    int pcm;    //0なら堆積無し、1なら片側堆積、2なら両側堆積
} Flag;

typedef struct _DataTable {
    Flag flag;
    Param par;
	std::array<Dataset, 4> dset;
} DataTable;



#endif