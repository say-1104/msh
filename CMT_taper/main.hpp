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

typedef struct _Param {
    double wl;      // 波長
    double dz;      // 伝搬方向の分割数
    int N_dset;     // データセットの個数
    double wst;     // テーパの先端幅
    int N_taper;    // テーパの接続個数
    int N_pcm;      // PCMの接続個数
    std::vector<std::pair<double, double>> taper;       // テーパの構造パラメータ(z, w)
    std::vector<std::pair<double, int>> pcm;            // PCM層の結晶情報(z, PCM)
    std::vector<std::tuple<double, double, int>> ZtoW;
} Param;

typedef struct _Dataset {
    int N;  // widthの個数
    LI<double, double> beta_even;
    LI<double, double> beta_odd;
    LI<double, double> beta_1;
    LI<double, double> beta_2;
} Dataset;

typedef struct _Flag {
    int pcm;    // 0なら堆積無し、1なら片側堆積、2なら両側堆積
    bool zw;    // trueならZtoWファイルを出力
} Flag;

typedef struct _DataTable {
    Flag flag;
    Param par;
	std::array<Dataset, 5> dset;
} DataTable;



#endif