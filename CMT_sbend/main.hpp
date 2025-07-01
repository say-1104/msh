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
#if __cplusplus < 201703L

#include <sys/types.h>
#include <sys/stat.h>

bool IsFileExist(const std::string& name) {
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0 && S_ISREG(buffer.st_mode));
}

#else // C++17以降の場合

#include <filesystem>

bool IsFileExist(const std::string& name) {
    return std::filesystem::is_regular_file(name);
}

#endif

constexpr std::complex<double> cj(0.0, 1.0);

#define ALLOCATION(a,c,n) \
		if ((a = (c *)malloc((n)*sizeof(c))) == NULL) {\
				fprintf(stderr, "can't allocate memory\n");\
				exit(EXIT_FAILURE);\
		}

#define SI 0
#define CPCM 1
#define APCM 2
#define PI 3.14159265359
typedef struct _Element {
    std::string elem_name;
    double L;
    double gfi;
    int pcm;
} Element;

typedef struct _Param {
    double wl;                                          // 波長
    double dz;                                          // 伝搬方向の分割数
    int N_dset;                                         // データセットの個数
    double gst;                                         // 始端の導波路間隔
    double w;                                           // 導波路幅
    int N_elem;                                         // 素子(dc or sbend)の接続個数
    std::vector<Element> elem;                          // 素子の構造パラメータ
    std::vector<std::tuple<double, double, int, int>> ZtoW;  // zに対する導波路間隔 int1 = 構造(dc = 1, sbend = 2) int2 = pcm
} Param;

typedef struct _Dataset {
    int N;  // widthの個数
    LI<double, double> beta_even;
    LI<double, double> beta_odd;
    LI<double, double> beta_1;
    LI<double, double> beta_2;
} Dataset;

typedef struct _Flag {
    bool zw;    // trueならZtoWファイルを出力
    bool sbend; 
} Flag;

typedef struct _DataTable {
    Flag flag;
    Param par;
	std::array<Dataset, 9> dset;
	std::array<Dataset, 3> dset_sbend;
} DataTable;



#endif