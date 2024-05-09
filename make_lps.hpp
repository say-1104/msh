#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <omp.h>
#include <vector>
using namespace std;

#define PI		    3.141592653589793

#define ZERO -1.0e-6

// struct
typedef struct _Parameter{
	double W;			//MMI幅
	double L;			//MMI長
	double Wport;		//導波路幅(狭)
	double Wtp;			//導波路幅(広)
	double g;			//導波路間隔
	double Ltp;			//テーパ長
	double hpcm; 		//PCM層の厚さ

	double lm;			//入出力導波路とPML間の距離(x方向)
	double wm;			//MMIとPML間の距離(y方向)
	double pml;			//PML幅
	double hco;			//コアの厚さ
	double ycent;		//y方向の中心
} Parameter;

typedef struct _Mosaic{
	double pix;			//ピクセル幅
	int divz;			//z方向のピクセル数
	int divx;			//x方向のピクセル数
	double x;
	double y;			//モザイク領域でx,yが共に最も小さい座標
} Mosaic;

typedef struct _PixData{
	int z;			//列
	int x;			//行
	vector<int> l;		//構成線番号
	vector<int> p;		//構成点番号
} PixData;

typedef struct _TP{
	double y;
	int num;
} TP;


//function
int Div(double x, double y){
	return((int)(round(x*10000)/round(y*10000)));
} // 小数の商計算

double REM(double x, double y){
	int tmp = (int)(round(x*10000))%(int)(round(y*10000));
	return(tmp/10000.0);
}// 小数の剰余計算

Parameter ParamInit(double W, double L, double Wport, double Wtp, double g, double Ltp, double hpcm, double lm, double wm, double pml, double hco, double ycent){
	Parameter par; 
	par.W = W;
	par.L = L;
	par.Wport = Wport;
	par.Wtp = Wtp;
	par.g = g;
	par.Ltp = Ltp;
	par.hpcm = hpcm;
	par.lm = lm;
	par.wm = wm;
	par.pml = pml;
	par.hco = hco;
	par.ycent = ycent;
	return(par);
}

Mosaic MosInit(double pix, int divx, int divz, double x, double y){
	Mosaic mos;
	mos.pix = pix;
	mos.divz = divz;
	mos.divx = divx;
	mos.x = x;
	mos.y = y;
	return(mos);
}

TP TPInit(double y, int num){
	TP tp;
	tp.y = y;
	tp.num = num;
	return(tp);
}

void printTP(vector<TP> tp){
	for(int i=0; i<tp.size();i++){
		printf("y : %f, ピクセル : %d番目\n", tp[i].y, tp[i].num);
	}
}

void Printvec(vector<int> a){
	for(int i=0; i<a.size(); i++){
		printf("%d ", a[i]);
	}
	printf("\n");
}