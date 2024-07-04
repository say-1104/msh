#include "fem3d.h"

// #define STRAIGHT_EDGE
#define CURVE

static const long long int ELECTRIC_FIELD = 1;

static const double COORD_1 = 1.0 / 4.0;
static const double COORD_2 = (7.0 - sqrt(15.0)) / 34.0;
static const double COORD_3 = (7.0 + sqrt(15.0)) / 34.0;
static const double COORD_4 = (13.0 + 3.0 * sqrt(15.0)) / 34.0;
static const double COORD_5 = (13.0 - 3.0 * sqrt(15.0)) / 34.0;
static const double COORD_6 = (5.0 - sqrt(15.0)) / 20.0;
static const double COORD_7 = (5.0 + sqrt(15.0)) / 20.0;

static const double WEIGHT1 = 8.0 / 405.0;
static const double WEIGHT2 = (2665.0 + 14.0 * sqrt(15.0)) / 226800.0;
static const double WEIGHT3 = (2665.0 - 14.0 * sqrt(15.0)) / 226800.0;
static const double WEIGHT4 = 5.0 / 567.0;

long long int elementMatrix(DataTable *data, double *x, double *y, double *z, 
				long long int matID,
				std::complex<double> sk[24][24], 
				std::complex<double> sm[24][24]) {
	double a[4], b[4], c[4], d[4], alpha[4];
	double ss[6], aa[6], bb[6], cc[6], ll[6];
	double area[4];
	static const long long int edge[6][2] = { { 0, 1 }, { 0, 2 }, { 0, 3 }, 
						{ 1, 2 }, { 3, 1 }, { 2, 3 } };
	double Ve, WW;
	double L1, L2, L3, L4;
	double L1x, L2x, L3x, L4x;
	double L1y, L2y, L3y, L4y;
	double L1z, L2z, L3z, L4z;
	double U[24], V[24], W[24];
	double Vx[24], Wx[24];
	double Wy[24], Uy[24];
	double Uz[24], Vz[24];
	double co[24];
	
	double pp[3][3], qq[3][3];
	std::complex<double> pp_s[3][3], qq_s[3][3];
	std::complex<double> t_sk[15][24][24], t_sm[15][24][24];
	
	/* --------PML--------------------------------------- */
	double tanD = data->pml.tanD;
	double mm = data->pml.m;
	double dx1 = data->pml.dx1, dx2 = data->pml.dx2;
	double dy1 = data->pml.dy1, dy2 = data->pml.dy2;
	double dz1 = data->pml.dz1, dz2 = data->pml.dz2;
	double pml_x1 = data->fem.min.x + dx1;
	double pml_x2 = data->fem.max.x - dx2;
	double pml_y1 = data->fem.min.y + dy1;
	double pml_y2 = data->fem.max.y - dy2;
	double pml_z1 = data->fem.min.z + dz1;
	double pml_z2 = data->fem.max.z - dz2;
	std::complex<double> sx, sy, sz;
	double xg, yg, zg;
	/* --------PML END----------------------------------- */

#ifdef STRAIGHT_EDGE
	/* - mod kakki 直辺要素になるように座標を上書き(3次元) - */
	x[4] = 0.5 * (x[0] + x[1]);
	y[4] = 0.5 * (y[0] + y[1]);
	z[4] = 0.5 * (z[0] + z[1]);
	x[5] = 0.5 * (x[1] + x[2]);
	y[5] = 0.5 * (y[1] + y[2]);
	z[5] = 0.5 * (z[1] + z[2]);
	x[6] = 0.5 * (x[0] + x[2]);
	y[6] = 0.5 * (y[0] + y[2]);
	z[6] = 0.5 * (z[0] + z[2]);
	x[7] = 0.5 * (x[0] + x[3]);
	y[7] = 0.5 * (y[0] + y[3]);
	z[7] = 0.5 * (z[0] + z[3]);
	x[8] = 0.5 * (x[1] + x[3]);
	y[8] = 0.5 * (y[1] + y[3]);
	z[8] = 0.5 * (z[1] + z[3]);
	x[9] = 0.5 * (x[2] + x[3]);
	y[9] = 0.5 * (y[2] + y[3]);
	z[9] = 0.5 * (z[2] + z[3]);
	/* -------------------------------------------------- */
#endif

	for (long long int i = 0; i < 3; ++i) {
		for (long long int j = 0; j < 3; ++j) {
			if (data->fem.unknown == ELECTRIC_FIELD) {
				pp[i][j] = data->par.mrmatrix[matID][i][j];
				qq[i][j] = data->par.ermatrix[matID][i][j];
			} else {
				pp[i][j] = data->par.ermatrix[matID][i][j];
				qq[i][j] = data->par.mrmatrix[matID][i][j];
			}
		}
	}
	InverseMatrix(pp);

#ifdef CURVE
	/* ----- 曲辺要素に対応した計算手法をとる ----------- */
	double J;
	double N[10], nx[10], ny[10], nz[10];
	double j11, j12, j13, j21, j22, j23, j31, j32, j33;
	double j11_i, j12_i, j13_i;
	double j21_i, j22_i, j23_i;
	double j31_i, j32_i, j33_i;
#endif

	alpha[0] = alpha[2] = 1.0;
	alpha[1] = alpha[3] = -1.0;
	
	// 全体座標から体積座標系への変換行列の係数に由来？ (有限要素法の基礎 p.57)
	for (long long int i = 0; i < 4; ++i) {
		long long int k = i, l = (i + 1) % 4, m = (i + 2) % 4, n = (i + 3) % 4;
		a[k] = alpha[k] * (x[l] * (y[m] * z[n] - y[n] * z[m]) 
					+ x[m] * (y[n] * z[l] - y[l] * z[n]) 
					+ x[n] * (y[l] * z[m] - y[m] * z[l]));
		b[k] = alpha[k] * (y[l] * (z[n] - z[m]) 
					+ y[m] * (z[l] - z[n]) 
					+ y[n] * (z[m] - z[l]));
		c[k] = alpha[k] * (z[l] * (x[n] - x[m]) 
					+ z[m] * (x[l] - x[n]) 
					+ z[n] * (x[m] - x[l]));
		d[k] = alpha[k] * (x[l] * (y[n] - y[m]) 
					+ x[m] * (y[l] - y[n]) 
					+ x[n] * (y[m] - y[l]));
	}

	Ve = (a[0] + a[1] + a[2] + a[3]) / 6.0;
	// if (Ve <= 0.0) {
	if (Ve < 1.0e-12) {
		fprintf(stderr, "volume is negative (%e)\n", Ve);
		for (long long int i = 0; i < 4; ++i) {
			fprintf(stderr, "     (%lf, %lf, %lf)\n", x[i], y[i], z[i]);
		}
		exit(EXIT_FAILURE);
	}
	
	for (long long int i = 0; i < 6; ++i) {
		aa[i] = x[edge[i][1]] - x[edge[i][0]];
		bb[i] = y[edge[i][1]] - y[edge[i][0]];
		cc[i] = z[edge[i][1]] - z[edge[i][0]];
		ll[i] = sqrt(aa[i] * aa[i] + bb[i] * bb[i] + cc[i] * cc[i]);
	}

	// 頂点[k]と向かい合う面の面積[k]
	area[0] = sqrt(pow(bb[3] * cc[4] - bb[4] * cc[3], 2.0) 
		 + pow(cc[3] * aa[4] - cc[4] * aa[3], 2.0) 
		 + pow(aa[3] * bb[4] - aa[4] * bb[3], 2.0));
	area[1] = sqrt(pow(bb[1] * cc[2] - bb[2] * cc[1], 2.0) 
		 + pow(cc[1] * aa[2] - cc[2] * aa[1], 2.0) 
		 + pow(aa[1] * bb[2] - aa[2] * bb[1], 2.0));
	area[2] = sqrt(pow(bb[0] * cc[2] - bb[2] * cc[0], 2.0) 
		 + pow(cc[0] * aa[2] - cc[2] * aa[0], 2.0) 
		 + pow(aa[0] * bb[2] - aa[2] * bb[0], 2.0));
	area[3] = sqrt(pow(bb[0] * cc[1] - bb[1] * cc[0], 2.0) 
		 + pow(cc[0] * aa[1] - cc[1] * aa[0], 2.0) 
		 + pow(aa[0] * bb[1] - aa[1] * bb[0], 2.0));

	for (long long int i = 0; i < 6; ++i) {
		if (cc[i] > ZERO 
		|| (std::abs(cc[i]) < ZERO && bb[i] > ZERO)
		|| (std::abs(cc[i]) < ZERO && std::abs(bb[i]) < ZERO && aa[i] > ZERO)) {
			ss[i] = 1.0;
		} else {
			ss[i] = -1.0;
		}
	}

#ifndef CURVE
	// 全体→体積の座標系変換行列の各係数は、x微分(要素B)、y微分(要素C)、z微分(要素D)に一致する
	L1x = b[0] / (6.0 * Ve);
	L1y = c[0] / (6.0 * Ve);
	L1z = d[0] / (6.0 * Ve);
	L2x = b[1] / (6.0 * Ve);
	L2y = c[1] / (6.0 * Ve);
	L2z = d[1] / (6.0 * Ve);
	L3x = b[2] / (6.0 * Ve);
	L3y = c[2] / (6.0 * Ve);
	L3z = d[2] / (6.0 * Ve);
	L4x = b[3] / (6.0 * Ve);
	L4y = c[3] / (6.0 * Ve);
	L4z = d[3] / (6.0 * Ve);

	/* ------------------------------ */
	if (data->fem.order == 1) {
		Vx[0] = L1x * L2y - L2x * L1y;
		Wx[0] = L1x * L2z - L2x * L1z;
		Vx[1] = L1x * L3y - L3x * L1y;
		Wx[1] = L1x * L3z - L3x * L1z;
		Vx[2] = L1x * L4y - L4x * L1y;
		Wx[2] = L1x * L4z - L4x * L1z;
		Vx[3] = L2x * L3y - L3x * L2y;
		Wx[3] = L2x * L3z - L3x * L2z;
		Vx[4] = L4x * L2y - L2x * L4y;
		Wx[4] = L4x * L2z - L2x * L4z;
		Vx[5] = L3x * L4y - L4x * L3y;
		Wx[5] = L3x * L4z - L4x * L3z;
		
		Uy[0] = L1y * L2x - L2y * L1x;
		Wy[0] = L1y * L2z - L2y * L1z;
		Uy[1] = L1y * L3x - L3y * L1x;
		Wy[1] = L1y * L3z - L3y * L1z;
		Uy[2] = L1y * L4x - L4y * L1x;
		Wy[2] = L1y * L4z - L4y * L1z;
		Uy[3] = L2y * L3x - L3y * L2x;
		Wy[3] = L2y * L3z - L3y * L2z;
		Uy[4] = L4y * L2x - L2y * L4x;
		Wy[4] = L4y * L2z - L2y * L4z;
		Uy[5] = L3y * L4x - L4y * L3x;
		Wy[5] = L3y * L4z - L4y * L3z;
		
		Uz[0] = L1z * L2x - L2z * L1x;
		Vz[0] = L1z * L2y - L2z * L1y;
		Uz[1] = L1z * L3x - L3z * L1x;
		Vz[1] = L1z * L3y - L3z * L1y;
		Uz[2] = L1z * L4x - L4z * L1x;
		Vz[2] = L1z * L4y - L4z * L1y;
		Uz[3] = L2z * L3x - L3z * L2x;
		Vz[3] = L2z * L3y - L3z * L2y;
		Uz[4] = L4z * L2x - L2z * L4x;
		Vz[4] = L4z * L2y - L2z * L4y;
		Uz[5] = L3z * L4x - L4z * L3x;
		Vz[5] = L3z * L4y - L4z * L3y;
		
		for (long long int i = 0; i < 6; ++i) {
			Uy[i] *= ss[i] * ll[i];
			Uz[i] *= ss[i] * ll[i];
			Vz[i] *= ss[i] * ll[i];
			Vx[i] *= ss[i] * ll[i];
			Wy[i] *= ss[i] * ll[i];
			Wx[i] *= ss[i] * ll[i];
		}
	} else {
		co[12] = 4.0 * area[0] / ll[5];
		co[13] = 4.0 * area[0] / ll[4];
		co[14] = 4.0 * area[0] / ll[3];
		co[15] = 4.0 * area[1] / ll[5];
		co[16] = 4.0 * area[1] / ll[1];
		co[17] = 4.0 * area[1] / ll[2];
		co[18] = 4.0 * area[2] / ll[0];
		co[19] = 4.0 * area[2] / ll[4];
		co[20] = 4.0 * area[2] / ll[2];
		co[21] = 4.0 * area[3] / ll[0];
		co[22] = 4.0 * area[3] / ll[1];
		co[23] = 4.0 * area[3] / ll[3];
		
		Vx[0] = L1x * L2y;
		Wx[0] = L1x * L2z;
		Vx[1] = L1x * L3y;
		Wx[1] = L1x * L3z;
		Vx[2] = L1x * L4y;
		Wx[2] = L1x * L4z;
		Vx[3] = L2x * L3y;
		Wx[3] = L2x * L3z;
		Vx[4] = L4x * L2y;
		Wx[4] = L4x * L2z;
		Vx[5] = L3x * L4y;
		Wx[5] = L3x * L4z;
		Vx[6] = L1y * L2x;
		Wx[6] = L1z * L2x;
		Vx[7] = L1y * L3x;
		Wx[7] = L1z * L3x;
		Vx[8] = L1y * L4x;
		Wx[8] = L1z * L4x;
		Vx[9] = L2y * L3x;
		Wx[9] = L2z * L3x;
		Vx[10] = L4y * L2x;
		Wx[10] = L4z * L2x;
		Vx[11] = L3y * L4x;
		Wx[11] = L3z * L4x;
		
		Uy[0] = L1y * L2x;
		Wy[0] = L1y * L2z;
		Uy[1] = L1y * L3x;
		Wy[1] = L1y * L3z;
		Uy[2] = L1y * L4x;
		Wy[2] = L1y * L4z;
		Uy[3] = L2y * L3x;
		Wy[3] = L2y * L3z;
		Uy[4] = L4y * L2x;
		Wy[4] = L4y * L2z;
		Uy[5] = L3y * L4x;
		Wy[5] = L3y * L4z;
		Uy[6] = L1x * L2y;
		Wy[6] = L1z * L2y;
		Uy[7] = L1x * L3y;
		Wy[7] = L1z * L3y;
		Uy[8] = L1x * L4y;
		Wy[8] = L1z * L4y;
		Uy[9] = L2x * L3y;
		Wy[9] = L2z * L3y;
		Uy[10] = L4x * L2y;
		Wy[10] = L4z * L2y;
		Uy[11] = L3x * L4y;
		Wy[11] = L3z * L4y;

		Uz[0] = L1z * L2x;
		Vz[0] = L1z * L2y;
		Uz[1] = L1z * L3x;
		Vz[1] = L1z * L3y;
		Uz[2] = L1z * L4x;
		Vz[2] = L1z * L4y;
		Uz[3] = L2z * L3x;
		Vz[3] = L2z * L3y;
		Uz[4] = L4z * L2x;
		Vz[4] = L4z * L2y;
		Uz[5] = L3z * L4x;
		Vz[5] = L3z * L4y;
		Uz[6] = L1x * L2z;
		Vz[6] = L1y * L2z;
		Uz[7] = L1x * L3z;
		Vz[7] = L1y * L3z;
		Uz[8] = L1x * L4z;
		Vz[8] = L1y * L4z;
		Uz[9] = L2x * L3z;
		Vz[9] = L2y * L3z;
		Uz[10] = L4x * L2z;
		Vz[10] = L4y * L2z;
		Uz[11] = L3x * L4z;
		Vz[11] = L3y * L4z;
		
		for (long long int i = 0; i < 12; ++i) {
			Uy[i] *= ll[i % 6];
			Uz[i] *= ll[i % 6];
			Vz[i] *= ll[i % 6];
			Vx[i] *= ll[i % 6];
			Wy[i] *= ll[i % 6];
			Wx[i] *= ll[i % 6];
		}
	}
	/* ------------------------------ */
#endif

	// 数値積分用の座標代入
	for (long long int l = 0; l < 15; ++l) {
		switch (l) {
		case 0:
			L1 = COORD_1;
			L2 = COORD_1;
			L3 = COORD_1;
			L4 = COORD_1;
			WW = WEIGHT1;
			break;
		case 1:
			L1 = COORD_2;
			L2 = COORD_2;
			L3 = COORD_2;
			L4 = COORD_4;
			WW = WEIGHT2;
			break;
		case 2:
			L1 = COORD_2;
			L2 = COORD_2;
			L3 = COORD_4;
			L4 = COORD_2;
			WW = WEIGHT2;
			break;
		case 3:
			L1 = COORD_2;
			L2 = COORD_4;
			L3 = COORD_2;
			L4 = COORD_2;
			WW = WEIGHT2;
			break;
		case 4:
			L1 = COORD_4;
			L2 = COORD_2;
			L3 = COORD_2;
			L4 = COORD_2;
			WW = WEIGHT2;
			break;
		case 5:
			L1 = COORD_3;
			L2 = COORD_3;
			L3 = COORD_3;
			L4 = COORD_5;
			WW = WEIGHT3;
			break;
		case 6:
			L1 = COORD_3;
			L2 = COORD_3;
			L3 = COORD_5;
			L4 = COORD_3;
			WW = WEIGHT3;
			break;
		case 7:
			L1 = COORD_3;
			L2 = COORD_5;
			L3 = COORD_3;
			L4 = COORD_3;
			WW = WEIGHT3;
			break;
		case 8:
			L1 = COORD_5;
			L2 = COORD_3;
			L3 = COORD_3;
			L4 = COORD_3;
			WW = WEIGHT3;
			break;
		case 9:
			L1 = COORD_6;
			L2 = COORD_6;
			L3 = COORD_7;
			L4 = COORD_7;
			WW = WEIGHT4;
			break;
		case 10:
			L1 = COORD_6;
			L2 = COORD_7;
			L3 = COORD_6;
			L4 = COORD_7;
			WW = WEIGHT4;
			break;
		case 11:
			L1 = COORD_7;
			L2 = COORD_6;
			L3 = COORD_6;
			L4 = COORD_7;
			WW = WEIGHT4;
			break;
		case 12:
			L1 = COORD_6;
			L2 = COORD_7;
			L3 = COORD_7;
			L4 = COORD_6;
			WW = WEIGHT4;
			break;
		case 13:
			L1 = COORD_7;
			L2 = COORD_6;
			L3 = COORD_7;
			L4 = COORD_6;
			WW = WEIGHT4;
			break;
		default:
			L1 = COORD_7;
			L2 = COORD_7;
			L3 = COORD_6;
			L4 = COORD_6;
			WW = WEIGHT4;
			break;
		}
		
	#ifdef CURVE
		if (data->fem.order == 1) {
			N[0] = L1;
			nx[0] = 1.0;
			ny[0] = 0.0;
			nz[0] = 0.0;
			N[1] = L2;
			nx[1] = 0.0;
			ny[1] = 1.0;
			nz[1] = 0.0;
			N[2] = L3;
			nx[2] = 0.0;
			ny[2] = 0.0;
			nz[2] = 1.0;
			N[3] = L4;
			nx[3] = -1.0;
			ny[3] = -1.0;
			nz[3] = -1.0;
		} else {
			// nxは∂N/∂L1, nyは∂N/∂L2... L4 = 1-L1-L2-L3を使用
			// N[]は直接使わないけどメモ書き?
			N[0] = L1 * (2.0 * L1 - 1.0);
			nx[0] = 4.0 * L1 - 1.0;
			ny[0] = 0.0;
			nz[0] = 0.0;
			N[1] = L2 * (2.0 * L2 - 1.0);
			nx[1] = 0.0;
			ny[1] = 4.0 * L2 - 1.0;
			nz[1] = 0.0;
			N[2] = L3 * (2.0 * L3 - 1.0);
			nx[2] = 0.0;
			ny[2] = 0.0;
			nz[2] = 4.0 * L3 - 1.0;
			N[3] = L4 * (2.0 * L4 - 1.0);
			nx[3] = -4.0 * L4 + 1.0;
			ny[3] = -4.0 * L4 + 1.0;
			nz[3] = -4.0 * L4 + 1.0;
			N[4] = 4.0 * L1 * L2;
			nx[4] = 4.0 * L2;
			ny[4] = 4.0 * L1;
			nz[4] = 0.0;
			N[5] = 4.0 * L2 * L3;
			nx[5] = 0.0;
			ny[5] = 4.0 * L3;
			nz[5] = 4.0 * L2;
			N[6] = 4.0 * L3 * L1;
			nx[6] = 4.0 * L3;
			ny[6] = 0.0;
			nz[6] = 4.0 * L1;
			N[7] = 4.0 * L1 * L4;
			nx[7] = 4.0 * (L4 - L1);
			ny[7] = -4.0 * L1;
			nz[7] = -4.0 * L1;
			N[8] = 4.0 * L2 * L4;
			nx[8] = -4.0 * L2;
			ny[8] = 4.0 * (L4 - L2);
			nz[8] = -4.0 * L2;
			N[9] = 4.0 * L3 * L4;
			nx[9] = -4.0 * L3;
			ny[9] = -4.0 * L3;
			nz[9] = 4.0 * (L4 - L3);
		}
		
		j11 = j12 = j13 = j21 = j22 = j23 = j31 = j32 = j33 = 0.0;
		for (long long int i = 0; i < data->fem.n_node_use; ++i) {
			j11 += nx[i] * x[i];
			j12 += nx[i] * y[i];
			j13 += nx[i] * z[i];
			j21 += ny[i] * x[i];
			j22 += ny[i] * y[i];
			j23 += ny[i] * z[i];
			j31 += nz[i] * x[i];
			j32 += nz[i] * y[i];
			j33 += nz[i] * z[i];
		}

		// ヤコビアン = (通常座標での体積 / 体積座標での(正規化)体積)
		J = j11 * j22 * j33 + j21 * j32 * j13 + j31 * j12 * j23 - j11 * j32
			* j23 - j21 * j12 * j33 - j31 * j22 * j13;

		// 曲辺要素に特異的な部分(要素体積の修正)
		// *question J / 6.0をする必要がある気がするが、ない方が結果は正しそう
		Ve = J;
		
		// 全体→体積の座標系変換行列で逆行列を掛けている
		// *question この式はあっているか、確認を
		j11_i = (j22 * j33 - j23 * j32);
		j12_i = (-j12 * j33 + j13 * j32);
		j13_i = (j12 * j23 - j13 * j22);
		j21_i = (-j21 * j33 + j23 * j31);
		j22_i = (j11 * j33 - j13 * j31);
		j23_i = (-j11 * j23 + j13 * j21);
		j31_i = (j21 * j32 - j22 * j31);
		j32_i = (-j11 * j32 + j12 * j31);
		j33_i = (j11 * j22 - j12 * j21);
		
		L1x = j11_i / J;
		L1y = j21_i / J;
		L1z = j31_i / J;
		L2x = j12_i / J;
		L2y = j22_i / J;
		L2z = j32_i / J;
		L3x = j13_i / J;
		L3y = j23_i / J;
		L3z = j33_i / J;
		L4x = -L1x - L2x - L3x;
		L4y = -L1y - L2y - L3y;
		L4z = -L1z - L2z - L3z;
		
		/* ------------------------------ */
		if (data->fem.order == 1) {
			Vx[0] = L1x * L2y - L2x * L1y;
			Wx[0] = L1x * L2z - L2x * L1z;
			Vx[1] = L1x * L3y - L3x * L1y;
			Wx[1] = L1x * L3z - L3x * L1z;
			Vx[2] = L1x * L4y - L4x * L1y;
			Wx[2] = L1x * L4z - L4x * L1z;
			Vx[3] = L2x * L3y - L3x * L2y;
			Wx[3] = L2x * L3z - L3x * L2z;
			Vx[4] = L4x * L2y - L2x * L4y;
			Wx[4] = L4x * L2z - L2x * L4z;
			Vx[5] = L3x * L4y - L4x * L3y;
			Wx[5] = L3x * L4z - L4x * L3z;
			
			Uy[0] = L1y * L2x - L2y * L1x;
			Wy[0] = L1y * L2z - L2y * L1z;
			Uy[1] = L1y * L3x - L3y * L1x;
			Wy[1] = L1y * L3z - L3y * L1z;
			Uy[2] = L1y * L4x - L4y * L1x;
			Wy[2] = L1y * L4z - L4y * L1z;
			Uy[3] = L2y * L3x - L3y * L2x;
			Wy[3] = L2y * L3z - L3y * L2z;
			Uy[4] = L4y * L2x - L2y * L4x;
			Wy[4] = L4y * L2z - L2y * L4z;
			Uy[5] = L3y * L4x - L4y * L3x;
			Wy[5] = L3y * L4z - L4y * L3z;
			
			Uz[0] = L1z * L2x - L2z * L1x;
			Vz[0] = L1z * L2y - L2z * L1y;
			Uz[1] = L1z * L3x - L3z * L1x;
			Vz[1] = L1z * L3y - L3z * L1y;
			Uz[2] = L1z * L4x - L4z * L1x;
			Vz[2] = L1z * L4y - L4z * L1y;
			Uz[3] = L2z * L3x - L3z * L2x;
			Vz[3] = L2z * L3y - L3z * L2y;
			Uz[4] = L4z * L2x - L2z * L4x;
			Vz[4] = L4z * L2y - L2z * L4y;
			Uz[5] = L3z * L4x - L4z * L3x;
			Vz[5] = L3z * L4y - L4z * L3y;
			
			for (long long int i = 0; i < 6; ++i) {
				Uy[i] *= ss[i] * ll[i];
				Uz[i] *= ss[i] * ll[i];
				Vz[i] *= ss[i] * ll[i];
				Vx[i] *= ss[i] * ll[i];
				Wy[i] *= ss[i] * ll[i];
				Wx[i] *= ss[i] * ll[i];
			}
		} else {
			co[12] = 4.0 * area[0] / ll[5];
			co[13] = 4.0 * area[0] / ll[4];
			co[14] = 4.0 * area[0] / ll[3];
			co[15] = 4.0 * area[1] / ll[5];
			co[16] = 4.0 * area[1] / ll[1];
			co[17] = 4.0 * area[1] / ll[2];
			co[18] = 4.0 * area[2] / ll[0];
			co[19] = 4.0 * area[2] / ll[4];
			co[20] = 4.0 * area[2] / ll[2];
			co[21] = 4.0 * area[3] / ll[0];
			co[22] = 4.0 * area[3] / ll[1];
			co[23] = 4.0 * area[3] / ll[3];
			
			Vx[0] = L1x * L2y;
			Wx[0] = L1x * L2z;
			Vx[1] = L1x * L3y;
			Wx[1] = L1x * L3z;
			Vx[2] = L1x * L4y;
			Wx[2] = L1x * L4z;
			Vx[3] = L2x * L3y;
			Wx[3] = L2x * L3z;
			Vx[4] = L4x * L2y;
			Wx[4] = L4x * L2z;
			Vx[5] = L3x * L4y;
			Wx[5] = L3x * L4z;
			Vx[6] = L1y * L2x;
			Wx[6] = L1z * L2x;
			Vx[7] = L1y * L3x;
			Wx[7] = L1z * L3x;
			Vx[8] = L1y * L4x;
			Wx[8] = L1z * L4x;
			Vx[9] = L2y * L3x;
			Wx[9] = L2z * L3x;
			Vx[10] = L4y * L2x;
			Wx[10] = L4z * L2x;
			Vx[11] = L3y * L4x;
			Wx[11] = L3z * L4x;
			
			Uy[0] = L1y * L2x;
			Wy[0] = L1y * L2z;
			Uy[1] = L1y * L3x;
			Wy[1] = L1y * L3z;
			Uy[2] = L1y * L4x;
			Wy[2] = L1y * L4z;
			Uy[3] = L2y * L3x;
			Wy[3] = L2y * L3z;
			Uy[4] = L4y * L2x;
			Wy[4] = L4y * L2z;
			Uy[5] = L3y * L4x;
			Wy[5] = L3y * L4z;
			Uy[6] = L1x * L2y;
			Wy[6] = L1z * L2y;
			Uy[7] = L1x * L3y;
			Wy[7] = L1z * L3y;
			Uy[8] = L1x * L4y;
			Wy[8] = L1z * L4y;
			Uy[9] = L2x * L3y;
			Wy[9] = L2z * L3y;
			Uy[10] = L4x * L2y;
			Wy[10] = L4z * L2y;
			Uy[11] = L3x * L4y;
			Wy[11] = L3z * L4y;
			
			Uz[0] = L1z * L2x;
			Vz[0] = L1z * L2y;
			Uz[1] = L1z * L3x;
			Vz[1] = L1z * L3y;
			Uz[2] = L1z * L4x;
			Vz[2] = L1z * L4y;
			Uz[3] = L2z * L3x;
			Vz[3] = L2z * L3y;
			Uz[4] = L4z * L2x;
			Vz[4] = L4z * L2y;
			Uz[5] = L3z * L4x;
			Vz[5] = L3z * L4y;
			Uz[6] = L1x * L2z;
			Vz[6] = L1y * L2z;
			Uz[7] = L1x * L3z;
			Vz[7] = L1y * L3z;
			Uz[8] = L1x * L4z;
			Vz[8] = L1y * L4z;
			Uz[9] = L2x * L3z;
			Vz[9] = L2y * L3z;
			Uz[10] = L4x * L2z;
			Vz[10] = L4y * L2z;
			Uz[11] = L3x * L4z;
			Vz[11] = L3y * L4z;
			
			for (long long int i = 0; i < 12; ++i) {
				Uy[i] *= ll[i % 6];
				Uz[i] *= ll[i % 6];
				Vz[i] *= ll[i % 6];
				Vx[i] *= ll[i % 6];
				Wy[i] *= ll[i % 6];
				Wx[i] *= ll[i % 6];
			}
		}
		/* ------------------------------ */
	#endif

		if (data->fem.order == 1) {
			U[0] = L1 * L2x - L2 * L1x;
			V[0] = L1 * L2y - L2 * L1y;
			W[0] = L1 * L2z - L2 * L1z;
			U[1] = L1 * L3x - L3 * L1x;
			V[1] = L1 * L3y - L3 * L1y;
			W[1] = L1 * L3z - L3 * L1z;
			U[2] = L1 * L4x - L4 * L1x;
			V[2] = L1 * L4y - L4 * L1y;
			W[2] = L1 * L4z - L4 * L1z;
			U[3] = L2 * L3x - L3 * L2x;
			V[3] = L2 * L3y - L3 * L2y;
			W[3] = L2 * L3z - L3 * L2z;
			U[4] = L4 * L2x - L2 * L4x;
			V[4] = L4 * L2y - L2 * L4y;
			W[4] = L4 * L2z - L2 * L4z;
			U[5] = L3 * L4x - L4 * L3x;
			V[5] = L3 * L4y - L4 * L3y;
			W[5] = L3 * L4z - L4 * L3z;
			
			for (long long int i = 0; i < 6; ++i) {
					U[i] *= ss[i] * ll[i];
					V[i] *= ss[i] * ll[i];
					W[i] *= ss[i] * ll[i];
			}
		} else {
			U[0] = L1 * L2x;
			V[0] = L1 * L2y;
			W[0] = L1 * L2z;
			U[1] = L1 * L3x;
			V[1] = L1 * L3y;
			W[1] = L1 * L3z;
			U[2] = L1 * L4x;
			V[2] = L1 * L4y;
			W[2] = L1 * L4z;
			U[3] = L2 * L3x;
			V[3] = L2 * L3y;
			W[3] = L2 * L3z;
			U[4] = L4 * L2x;
			V[4] = L4 * L2y;
			W[4] = L4 * L2z;
			U[5] = L3 * L4x;
			V[5] = L3 * L4y;
			W[5] = L3 * L4z;
			U[6] = L1x * L2;
			V[6] = L1y * L2;
			W[6] = L1z * L2;
			U[7] = L1x * L3;
			V[7] = L1y * L3;
			W[7] = L1z * L3;
			U[8] = L1x * L4;
			V[8] = L1y * L4;
			W[8] = L1z * L4;
			U[9] = L2x * L3;
			V[9] = L2y * L3;
			W[9] = L2z * L3;
			U[10] = L4x * L2;
			V[10] = L4y * L2;
			W[10] = L4z * L2;
			U[11] = L3x * L4;
			V[11] = L3y * L4;
			W[11] = L3z * L4;
			U[12] = L3 * L4 * L2x;
			V[12] = L3 * L4 * L2y;
			W[12] = L3 * L4 * L2z;
			U[13] = L4 * L2 * L3x;
			V[13] = L4 * L2 * L3y;
			W[13] = L4 * L2 * L3z;
			U[14] = L2 * L3 * L4x;
			V[14] = L2 * L3 * L4y;
			W[14] = L2 * L3 * L4z;
			U[15] = L4 * L3 * L1x;
			V[15] = L4 * L3 * L1y;
			W[15] = L4 * L3 * L1z;
			U[16] = L3 * L1 * L4x;
			V[16] = L3 * L1 * L4y;
			W[16] = L3 * L1 * L4z;
			U[17] = L1 * L4 * L3x;
			V[17] = L1 * L4 * L3y;
			W[17] = L1 * L4 * L3z;
			U[18] = L1 * L2 * L4x;
			V[18] = L1 * L2 * L4y;
			W[18] = L1 * L2 * L4z;
			U[19] = L2 * L4 * L1x;
			V[19] = L2 * L4 * L1y;
			W[19] = L2 * L4 * L1z;
			U[20] = L4 * L1 * L2x;
			V[20] = L4 * L1 * L2y;
			W[20] = L4 * L1 * L2z;
			U[21] = L2 * L1 * L3x;
			V[21] = L2 * L1 * L3y;
			W[21] = L2 * L1 * L3z;
			U[22] = L1 * L3 * L2x;
			V[22] = L1 * L3 * L2y;
			W[22] = L1 * L3 * L2z;
			U[23] = L3 * L2 * L1x;
			V[23] = L3 * L2 * L1y;
			W[23] = L3 * L2 * L1z;
			
			Vx[12] = (L3x * L4 + L3 * L4x) * L2y;
			Wx[12] = (L3x * L4 + L3 * L4x) * L2z;
			Vx[13] = (L4x * L2 + L4 * L2x) * L3y;
			Wx[13] = (L4x * L2 + L4 * L2x) * L3z;
			Vx[14] = (L2x * L3 + L2 * L3x) * L4y;
			Wx[14] = (L2x * L3 + L2 * L3x) * L4z;
			Vx[15] = (L4x * L3 + L4 * L3x) * L1y;
			Wx[15] = (L4x * L3 + L4 * L3x) * L1z;
			Vx[16] = (L3x * L1 + L3 * L1x) * L4y;
			Wx[16] = (L3x * L1 + L3 * L1x) * L4z;
			Vx[17] = (L1x * L4 + L1 * L4x) * L3y;
			Wx[17] = (L1x * L4 + L1 * L4x) * L3z;
			Vx[18] = (L1x * L2 + L1 * L2x) * L4y;
			Wx[18] = (L1x * L2 + L1 * L2x) * L4z;
			Vx[19] = (L2x * L4 + L2 * L4x) * L1y;
			Wx[19] = (L2x * L4 + L2 * L4x) * L1z;
			Vx[20] = (L4x * L1 + L4 * L1x) * L2y;
			Wx[20] = (L4x * L1 + L4 * L1x) * L2z;
			Vx[21] = (L2x * L1 + L2 * L1x) * L3y;
			Wx[21] = (L2x * L1 + L2 * L1x) * L3z;
			Vx[22] = (L1x * L3 + L1 * L3x) * L2y;
			Wx[22] = (L1x * L3 + L1 * L3x) * L2z;
			Vx[23] = (L3x * L2 + L3 * L2x) * L1y;
			Wx[23] = (L3x * L2 + L3 * L2x) * L1z;
			
			Uy[12] = (L3y * L4 + L3 * L4y) * L2x;
			Wy[12] = (L3y * L4 + L3 * L4y) * L2z;
			Uy[13] = (L4y * L2 + L4 * L2y) * L3x;
			Wy[13] = (L4y * L2 + L4 * L2y) * L3z;
			Uy[14] = (L2y * L3 + L2 * L3y) * L4x;
			Wy[14] = (L2y * L3 + L2 * L3y) * L4z;
			Uy[15] = (L4y * L3 + L4 * L3y) * L1x;
			Wy[15] = (L4y * L3 + L4 * L3y) * L1z;
			Uy[16] = (L3y * L1 + L3 * L1y) * L4x;
			Wy[16] = (L3y * L1 + L3 * L1y) * L4z;
			Uy[17] = (L1y * L4 + L1 * L4y) * L3x;
			Wy[17] = (L1y * L4 + L1 * L4y) * L3z;
			Uy[18] = (L1y * L2 + L1 * L2y) * L4x;
			Wy[18] = (L1y * L2 + L1 * L2y) * L4z;
			Uy[19] = (L2y * L4 + L2 * L4y) * L1x;
			Wy[19] = (L2y * L4 + L2 * L4y) * L1z;
			Uy[20] = (L4y * L1 + L4 * L1y) * L2x;
			Wy[20] = (L4y * L1 + L4 * L1y) * L2z;
			Uy[21] = (L2y * L1 + L2 * L1y) * L3x;
			Wy[21] = (L2y * L1 + L2 * L1y) * L3z;
			Uy[22] = (L1y * L3 + L1 * L3y) * L2x;
			Wy[22] = (L1y * L3 + L1 * L3y) * L2z;
			Uy[23] = (L3y * L2 + L3 * L2y) * L1x;
			Wy[23] = (L3y * L2 + L3 * L2y) * L1z;
			
			Uz[12] = (L3z * L4 + L3 * L4z) * L2x;
			Vz[12] = (L3z * L4 + L3 * L4z) * L2y;
			Uz[13] = (L4z * L2 + L4 * L2z) * L3x;
			Vz[13] = (L4z * L2 + L4 * L2z) * L3y;
			Uz[14] = (L2z * L3 + L2 * L3z) * L4x;
			Vz[14] = (L2z * L3 + L2 * L3z) * L4y;
			Uz[15] = (L4z * L3 + L4 * L3z) * L1x;
			Vz[15] = (L4z * L3 + L4 * L3z) * L1y;
			Uz[16] = (L3z * L1 + L3 * L1z) * L4x;
			Vz[16] = (L3z * L1 + L3 * L1z) * L4y;
			Uz[17] = (L1z * L4 + L1 * L4z) * L3x;
			Vz[17] = (L1z * L4 + L1 * L4z) * L3y;
			Uz[18] = (L1z * L2 + L1 * L2z) * L4x;
			Vz[18] = (L1z * L2 + L1 * L2z) * L4y;
			Uz[19] = (L2z * L4 + L2 * L4z) * L1x;
			Vz[19] = (L2z * L4 + L2 * L4z) * L1y;
			Uz[20] = (L4z * L1 + L4 * L1z) * L2x;
			Vz[20] = (L4z * L1 + L4 * L1z) * L2y;
			Uz[21] = (L2z * L1 + L2 * L1z) * L3x;
			Vz[21] = (L2z * L1 + L2 * L1z) * L3y;
			Uz[22] = (L1z * L3 + L1 * L3z) * L2x;
			Vz[22] = (L1z * L3 + L1 * L3z) * L2y;
			Uz[23] = (L3z * L2 + L3 * L2z) * L1x;
			Vz[23] = (L3z * L2 + L3 * L2z) * L1y;
			
			for (long long int i = 0; i < 12; ++i) {
				U[i] *= ll[i % 6];
				V[i] *= ll[i % 6];
				W[i] *= ll[i % 6];
			}
			for (long long int i = 12; i < data->fem.n_en; ++i) {
				U[i] *= co[i];
				V[i] *= co[i];
				W[i] *= co[i];
				Vx[i] *= co[i];
				Wx[i] *= co[i];
				Uy[i] *= co[i];
				Wy[i] *= co[i];
				Uz[i] *= co[i];
				Vz[i] *= co[i];
			}
		}
		
		/* --------PML--------------------------------------- */
		xg = L1 * x[0] + L2 * x[1] + L3 * x[2] + L4 * x[3];
		yg = L1 * y[0] + L2 * y[1] + L3 * y[2] + L4 * y[3];
		zg = L1 * z[0] + L2 * z[1] + L3 * z[2] + L4 * z[3];
		
		if (xg > pml_x2 && dx2 != 0) {
			sx = std::complex<double>(1.0, -std::pow((xg - pml_x2) / dx2, mm) * tanD);
		} else if (xg < pml_x1 && dx1 != 0) {
			sx = std::complex<double>(1.0, -std::pow((pml_x1 - xg) / dx1, mm) * tanD);
		} else {
			sx = 1.0;
		}
		if (yg > pml_y2 && dy2 != 0) {
			sy = std::complex<double>(1.0, -std::pow((yg - pml_y2) / dy2, mm) * tanD);
		} else if (yg < pml_y1 && dy1 != 0) {
			sy = std::complex<double>(1.0, -std::pow((pml_y1 - yg) / dy1, mm) * tanD);
		} else {
			sy = 1.0;
		}
		if (zg > pml_z2 && dz2 != 0) {
			sz = std::complex<double>(1.0, -std::pow((zg - pml_z2) / dz2, mm) * tanD);
		} else if (zg < pml_z1 && dz1 != 0) {
			sz = std::complex<double>(1.0, -std::pow((pml_z1 - zg) / dz1, mm) * tanD);
		} else {
			sz = 1.0;
		}

		for (long long int i = 0; i < 3; ++i) {
			for (long long int j = 0; j < 3; ++j) {
				pp_s[i][j] = pp[i][j];
				qq_s[i][j] = qq[i][j];
			}
		}

		pp_s[0][0] /= sy * sz / sx;
		pp_s[1][1] /= sz * sx / sy;
		pp_s[2][2] /= sx * sy / sz;
		qq_s[0][0] *= sy * sz / sx;
		qq_s[1][1] *= sz * sx / sy;
		qq_s[2][2] *= sx * sy / sz;
		/* --------PML END----------------------------------- */
		
		for (long long int i = 0; i < data->fem.n_en; ++i) {
			for (long long int j = 0; j < data->fem.n_en; ++j) {
				t_sk[l][i][j] 
					= Ve * WW * ((Wy[i] - Vz[i]) * (Wy[j] - Vz[j]) * pp_s[0][0] 
							+ (Uz[i] - Wx[i]) * (Uz[j] - Wx[j]) * pp_s[1][1] 
							+ (Vx[i] - Uy[i]) * (Vx[j] - Uy[j]) * pp_s[2][2]);
				t_sm[l][i][j] 
					= Ve * WW * (U[i] * U[j] * qq_s[0][0] 
							+ V[i] * V[j] * qq_s[1][1] 
							+ W[i] * W[j] * qq_s[2][2]);
			}
		}

	} // loop end of l (numerical integration)
	
	for (long long int i = 0; i < data->fem.n_en; ++i) {
		for (long long int j = 0; j < data->fem.n_en; ++j) {
			sk[i][j] = sm[i][j] = 0.0;
		}
	}
	
	for (long long int i = 0; i < data->fem.n_en; ++i) {
		for (long long int j = 0; j < data->fem.n_en; ++j) {
			for (long long int l = 0; l < 15; ++l) {
				sk[i][j] += t_sk[l][i][j];
				sm[i][j] += t_sm[l][i][j];
			}
		}
	}

	return 0;
}
