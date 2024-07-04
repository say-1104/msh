#include "fem3d.h"

/** 2022.01 変更点(By muratsubaki)
 *	出力画像の解像度は、絶対値ではなく材料依存にしたかったので修正。
 * 	ホール部は.slvではいびつな形状になるので、むしろ自作の方がよさそう
 * 		→ しかし、構造確認においてはslvファイルは有効なので、出力しておきたい。
 *　*/

// #define DIV_MAX 639 // 出力画像のサイズ
// #define DIV_MAX 319 // 出力画像のサイズ

#define NEARLY_ZERO 1.0e-6

#define NX 250
#define NY 250
#define NZ 250

/**
 * pict->element[ix][iz]に、画像上の各位置における参照要素を格納する。
 */
void detectElementForIntpol (DataTable *data) {
	FEM *fem = &(data->fem);
	Param *par = &(data->par);
	Picture *pict = &(data->pict);
	long long int i, j, k, l, m, n, ii, itmp;
	long long int ix, iz;
	long long int ecount;
	long long int *elist;
	double x0, y0, z0;
	double x[10], y[10], z[10];
	long long int **pixF;
	double a[4], b[4], c[4], d[4], alpha[4];
	double Ve;
	double L1, L2, L3, L4;
	XYZdouble emin, emax, div, ref1, ref2;

	ALLOCATION(pict->absxz, XYZdouble *, pict->nx)
	ALLOCATION(pict->realxz, XYZdouble *, pict->nx)
	for (i = 0; i < pict->nx; i++) {
		ALLOCATION(pict->absxz[i], XYZdouble, pict->nz)
		ALLOCATION(pict->realxz[i], XYZdouble, pict->nz)
		for (j = 0; j < pict->nz; j++) {
			pict->absxz[i][j].x = pict->absxz[i][j].y = pict->absxz[i][j].z = 0.0;
			pict->realxz[i][j].x = pict->realxz[i][j].y = pict->realxz[i][j].z = 0.0;
		}
	}
	ALLOCATION(pict->f_flag, long long int *, pict->nx)
	ALLOCATION(pict->element, long long int *, pict->nx)
	for (i = 0; i < pict->nx; i++) {
		ALLOCATION(pict->f_flag[i], long long int, pict->nz)
		ALLOCATION(pict->element[i], long long int, pict->nz)
	}

	div.x = (fem->max.x - fem->min.x) / pict->nx;
	div.y = (fem->max.y - fem->min.y) / pict->ny;
	div.z = (fem->max.z - fem->min.z) / pict->nz;
	y0 = par->out_slv_y_cood;

	ecount = 0;
	for (k = 0; k < fem->ne; k++) {//参照するy座標をまたぐ要素数をカウント
		emin.y = fem->max.y; emax.y = fem->min.y;
		for (i = 0; i < 4; i++) {
			itmp = fem->element[k].kn[i];
			y[i] = fem->node[itmp].y;
			if (emin.y > y[i]) emin.y = y[i];
			if (emax.y < y[i]) emax.y = y[i];
		}
		if (emin.y - 1.0e-6 < y0 && y0 < emax.y + 1.0e-6) ecount++;
	}
	pict->ecount = ecount;
	fprintf(stderr, "ecount = %lld\n", pict->ecount);
	ALLOCATION(elist, long long int, pict->ecount)
	for (i = 0; i < pict->ecount; i++) elist[i] = -1;
	ecount = 0;
	for (k = 0; k < fem->ne; k++) {//参照するy座標をまたぐ要素を抽出
		emin.y = fem->max.y; emax.y = fem->min.y;
		for (i = 0; i < 4; i++) {
			itmp = fem->element[k].kn[i];
			y[i] = fem->node[itmp].y;
			if (emin.y > y[i]) emin.y = y[i];
			if (emax.y < y[i]) emax.y = y[i];
		}
		if (emin.y - 1.0e-6 < y0 && y0 < emax.y + 1.0e-6) {
			elist[ecount] = k;//要素番号格納
			ecount++;
		}
	}
	for (i = 0; i < pict->ecount; i++) {
		if (elist[i] == -1) {
			fprintf(stderr, "Error in elist\n");
			exit(EXIT_FAILURE);
		}
	}
	fprintf(stderr, "---------- extract elements ----------\n");

	ALLOCATION(pixF, long long int *, fem->ne)
	for (k = 0; k < fem->ne; k++) {
		ALLOCATION(pixF[k], long long int, 4)
		for (j = 0; j < 4; j++) pixF[k][j] = -1;
	}
	fprintf(stderr, "---------- initialize pixF ----------\n");
	for (ix = 0; ix < pict->nx; ix++) {
		for (iz = 0; iz < pict->nz; iz++) pict->element[ix][iz] = -1;
	}
	fprintf(stderr, "---------- initialize pict->element ----------\n");
	for (j = 0; j < pict->ecount; j++) {//pixFの作成
		k = elist[j];
		emin.x = fem->max.x; emax.x = fem->min.x;
		emin.y = fem->max.y; emax.y = fem->min.y;
		emin.z = fem->max.z; emax.z = fem->min.z;
		for (i = 0; i < 4; i++) {
			itmp = fem->element[k].kn[i];
			x[i] = fem->node[itmp].x;
			y[i] = fem->node[itmp].y;
			z[i] = fem->node[itmp].z;
			if (emin.x > x[i]) emin.x = x[i];
			if (emax.x < x[i]) emax.x = x[i];
			if (emin.y > y[i]) emin.y = y[i];
			if (emax.y < y[i]) emax.y = y[i];
			if (emin.z > z[i]) emin.z = z[i];
			if (emax.z < z[i]) emax.z = z[i];
		}

		// 各要素について、座標位置を格納していく。
		for (ix = 0; ix < pict->nx; ix++) {
			ref1.x = fem->min.x + ix * div.x;
			ref2.x = ref1.x + div.x;
			if (ref1.x - div.x/10.0 < emin.x && emin.x < ref2.x + div.x/10.0) {
				if (pixF[k][0] == -1) pixF[k][0] = ix;
			}
			if (ref1.x - div.x/10.0 < emax.x && emax.x < ref2.x + div.x/10.0) {
				if (pixF[k][1] == -1) pixF[k][1] = ix;
			}
		}
		for (iz = 0; iz < pict->nz; iz++) {
			ref1.z = fem->min.z + iz * div.z;
			ref2.z = ref1.z + div.z;
			if (ref1.z - div.z/10.0 < emin.z && emin.z < ref2.z + div.z/10.0) {
				if (pixF[k][2] == -1) pixF[k][2] = iz;
			}
			if (ref1.z - div.z/10.0 < emax.z && emax.z < ref2.z + div.z/10.0) {
				if (pixF[k][3] == -1) pixF[k][3] = iz;
			}
		}
	}
	fprintf(stderr, "---------- create pixF ----------\n");
	for (j = 0; j < pict->ecount; j++) {
		k = elist[j];
		if (pixF[k][0] == -1) fprintf(stderr, "pict->element error 0.\n");
		if (pixF[k][1] == -1) fprintf(stderr, "pict->element error 1.\n");
		if (pixF[k][2] == -1) fprintf(stderr, "pict->element error 2.\n");
		if (pixF[k][3] == -1) fprintf(stderr, "pict->element error 3.\n");
	}

	alpha[0] = alpha[2] = 1.0;
	alpha[1] = alpha[3] = -1.0;
	for (j = 0; j < pict->ecount; j++) {//pict->elementの作成
		k = elist[j];
		for (i = 0; i < 4; i++) {
			itmp = fem->element[k].kn[i];
			x[i] = fem->node[itmp].x;
			y[i] = fem->node[itmp].y;
			z[i] = fem->node[itmp].z;
		}
		for (ix = pixF[k][0]; ix <= pixF[k][1]; ix++) {
			for (iz = pixF[k][2]; iz <= pixF[k][3]; iz++) {
				if (pict->element[ix][iz] != -1) continue;
				x0 = fem->min.x + (ix + 0.5) * div.x;
				z0 = fem->min.z + (iz + 0.5) * div.z;
				for (i = 0; i < 4; i++) {
					ii = i;
					l = (i + 1) % 4;
					m = (i + 2) % 4;
					n = (i + 3) % 4;
					a[ii] = alpha[ii] * (x[l] * (y[m] * z[n] - y[n] * z[m]) + x[m] * (y[n] * z[l] - y[l] * z[n]) + x[n] * (y[l] * z[m] - y[m] * z[l]));
					b[ii] = alpha[ii] * (y[l] * (z[n] - z[m]) + y[m] * (z[l] - z[n]) + y[n] * (z[m] - z[l]));
					c[ii] = alpha[ii] * (z[l] * (x[n] - x[m]) + z[m] * (x[l] - x[n]) + z[n] * (x[m] - x[l]));
					d[ii] = alpha[ii] * (x[l] * (y[n] - y[m]) + x[m] * (y[l] - y[n]) + x[n] * (y[m] - y[l]));
				}
				Ve = (a[0] + a[1] + a[2] + a[3]) / 6.0;
				L1 = (a[0] + b[0] * x0 + c[0] * y0 + d[0] * z0) / (6.0 * Ve);
				L2 = (a[1] + b[1] * x0 + c[1] * y0 + d[1] * z0) / (6.0 * Ve);
				L3 = (a[2] + b[2] * x0 + c[2] * y0 + d[2] * z0) / (6.0 * Ve);
				L4 = (a[3] + b[3] * x0 + c[3] * y0 + d[3] * z0) / (6.0 * Ve);
				if ((L1 >= -NEARLY_ZERO) && (L2 >= -NEARLY_ZERO) && (L3 >= -NEARLY_ZERO) && (L4 >= -NEARLY_ZERO)) {
					pict->element[ix][iz] = k;
				}
			}
		}
	}
	fprintf(stderr, "---------- create pict->element ----------\n");
	for (ix = 0; ix < pict->nx; ix++) {
		for (iz = 0; iz < pict->nz; iz++) {
			if (pict->element[ix][iz] == -1) {
				fprintf(stderr, "pict->element error.\n");
				exit(EXIT_FAILURE);
			}
		}
	}
	for (k = 0; k < fem->ne; k++) free(pixF[k]);
	free(pixF);
	free(elist);
	fprintf(stderr, "detectElementForIntpol END\n");
	return;
}

void intpol(DataTable *data, long long int nn, long long int mode, long long int realflag) {
	Picture *pict = &(data->pict);
	static long long int c_count = 0;
	static double xmin, xmax, ymin, ymax, zmin, zmax;
	static long long int *table[NX][NZ], table_n[NX][NZ];

	long long int ih, iv, nh, nv;
	double hmin, hmax, vmin, vmax;
	double hdiv, vdiv;
			
	long long int i, j, k, l, m, n, ii;
	long long int flag, itmp;
	XYZdouble* mat;
	CXYZdouble* phi;
	std::complex<double> phix, phiy, phiz;
	std::complex<double> phix2, phiy2, phiz2;
	double xdiv, ydiv, zdiv;
	double a[4], b[4], c[4], d[4], alpha[4];
	double U[24], V[24], W[24];
	double Vx[24], Wx[24];
	double Wy[24], Uy[24];
	double Uz[24], Vz[24];
	double ss[6], aa[6], bb[6], cc[6], ll[6];
	static long long int edge[6][2] = { { 0, 1 }, { 0, 2 }, { 0, 3 }, 
				{ 1, 2 }, { 3, 1 }, { 2, 3 } };
	double Ve;
	double L1, L2, L3, L4;
	double L1x, L2x, L3x, L4x;
	double L1y, L2y, L3y, L4y;
	double L1z, L2z, L3z, L4z;
	double x0, y0, z0;
	double x[10], y[10], z[10];
	std::complex<double> vv[24];
	double ZZ = 376.730313461;
	double k0 = 2.0 * M_PI / data->par.wavelength;
	std::complex<double> cj = std::complex<double>(0, 1);
	
	FEM *fem = &(data->fem);
	Param *par = &(data->par);
	
	double e_xmin, e_xmax, e_zmin, e_zmax;
	long long int region_x, region_z;
	long long int start_x, start_z, end_x, end_z;
	long long int iii, kk;
	long long int matID;

	int ii_pixel, jj_pixel;
	
	double area[4];
	double co[24];
	
	FILE *fp1, *fp2, *fp3;
	char string[256];
	char string2[256];

	if(c_count == 0){
		xmin = fem->min.x; xmax = fem->max.x;
		ymin = fem->min.y; ymax = fem->max.y;
		zmin = fem->min.z; zmax = fem->max.z;
		
		/**
		 *	ZX断面でテーブル作成 ← 要素検索を高速に行うため??
		 *	コアフレーム描画用にpict->element[][]も存在しているが、こちらはy = y0平面上の要素のみを格納している
		 */
		for (i = 0; i < NX; i++) {
			for (j = 0; j < NZ; j++) {
				table_n[i][j] = 0;
					
				for (ii = 0; ii < fem->ne; ii++) {
					e_xmin = xmax;
					e_xmax = xmin;
					e_zmin = zmax;
					e_zmax = zmin;
					for (i = 0; i < 4; i++) {
						itmp = fem->element[ii].kn[i];
						x[i] = fem->node[itmp].x;
						y[i] = fem->node[itmp].y;
						z[i] = fem->node[itmp].z;
						if (e_xmin > x[i]) e_xmin = x[i];
						if (e_xmax < x[i]) e_xmax = x[i];
						if (e_zmin > z[i]) e_zmin = z[i];
						if (e_zmax < z[i]) e_zmax = z[i];
					}
					start_x = (long long int) ((e_xmin - xmin) / (xmax - xmin) * NX);
					start_z = (long long int) ((e_zmin - zmin) / (zmax - zmin) * NZ);
					end_x = (long long int) ((e_xmax - xmin) / (xmax - xmin) * NX) + 1;
					end_z = (long long int) ((e_zmax - zmin) / (zmax - zmin) * NZ) + 1;
					if (end_x > NX) end_x = NX;
					if (end_z > NZ) end_z = NZ;
					for (i = start_x; i < end_x; i++) {
						for (j = start_z; j < end_z; j++) {
							table_n[i][j]++;
						}
					}
				}
				/* ------------------------------------------ */
				
				for (i = 0; i < NX; i++) {
					for (j = 0; j < NZ; j++) {
						ALLOCATION(table[i][j], long long int, (table_n[i][j]+1));
						table_n[i][j] = 0;
					}
				}
				
				for (ii = 0; ii < fem->ne; ii++) {
					e_xmin = xmax;
					e_xmax = xmin;
					e_zmin = zmax;
					e_zmax = zmin;
					for (i = 0; i < 4; i++) {
						itmp = fem->element[ii].kn[i];
						x[i] = fem->node[itmp].x;
						y[i] = fem->node[itmp].y;
						z[i] = fem->node[itmp].z;
						if (e_xmin > x[i]) e_xmin = x[i];
						if (e_xmax < x[i]) e_xmax = x[i];
						if (e_zmin > z[i]) e_zmin = z[i];
						if (e_zmax < z[i]) e_zmax = z[i];
					}
					start_x = (long long int) ((e_xmin - xmin) / (xmax - xmin) * NX);
					start_z = (long long int) ((e_zmin - zmin) / (zmax - zmin) * NZ);
					end_x = (long long int) ((e_xmax - xmin) / (xmax - xmin) * NX) + 1;
					end_z = (long long int) ((e_zmax - zmin) / (zmax - zmin) * NZ) + 1;
					if (end_x > NX) end_x = NX;
					if (end_z > NZ) end_z = NZ;
					if (start_x < 0) start_x = 0;
					if (start_z < 0) start_z = 0;
					for (i = start_x; i < end_x; i++) {
						for (j = start_z; j < end_z; j++) {
							table[i][j][table_n[i][j]] = ii;
							table_n[i][j]++;
						}
					}
				}
				
				kk = 0;
				for (i = 0; i < NX; i++) {
					for (j = 0; j < NZ; j++) {
						if(table_n[i][j] > kk)
							kk = table_n[i][j];
					}
				}
				fprintf(stderr, "max intpol-table = %lld\n", kk);
				c_count++;
			}
		}
	} // end of c_count, end of making table
					
	/* ------------------------------------------ */
	// ファイル出力開始
	if(nn < 0 || nn > 2) return;
		
	alpha[0] = alpha[2] = 1.0;
	alpha[1] = alpha[3] = -1.0;
	
	if(nn == 0){
		hmin = xmin;	vmin = ymin;
		hmax = xmax;	vmax = ymax;
		nh = pict->nx;	nv = pict->ny;
	}
	if(nn == 1){
		hmin = zmin;	vmin = ymin;
		hmax = zmax;	vmax = ymax;
		nh = pict->nz;	nv = pict->ny;
	}
	if(nn == 2){
		hmin = zmin;	vmin = xmin;
		hmax = zmax;	vmax = xmax;
		nh = pict->nz;	nv = pict->nx;
	}
	
	hdiv = (hmax - hmin) / nh;
	vdiv = (vmax - vmin) / nv;

	ALLOCATION(mat, XYZdouble,  nh*nv)
	ALLOCATION(phi, CXYZdouble, nh*nv)
	for (i = 0; i < nh*nv; i++) {
		mat[i].x = mat[i].y = mat[i].z = 0.0;
		phi[i].x = phi[i].y = phi[i].z = 0.0;
	}
	
	kk = 0;
	for (ih = 0; ih < nh; ih++) {
		for (iv = 0; iv < nv; iv++) {
			flag = 0;
			phix = phiy = phiz = 0.0;
			phix2 = phiy2 = phiz2 = 0.0;
			if(nn == 0){ // xy
				x0 = hmin + (hdiv * (ih+0.5));
				y0 = vmin + (vdiv * (iv+0.5));
				z0 = data->par.out_slv_z_cood;
			}
			if(nn == 1){ // zy
				z0 = hmin + (hdiv * (ih+0.5));
				y0 = vmin + (vdiv * (iv+0.5));
				x0 = data->par.out_slv_x_cood;
			}
			if(nn == 2){ // zx
				z0 = hmin + (hdiv * (ih+0.5));
				x0 = vmin + (vdiv * (iv+0.5));
				y0 = data->par.out_slv_y_cood;
			}
			
			region_x = (long long int) ((x0 - xmin) / (xmax - xmin) * NX);
			region_z = (long long int) ((z0 - zmin) / (zmax - zmin) * NZ);
			if (region_x >= NX) region_x = NX - 1;
			if (region_z >= NZ) region_z = NZ - 1;
			for (iii = 0; iii < table_n[region_x][region_z]; iii++) {
				ii = table[region_x][region_z][iii];
				for (i = 0; i < 4; i++) {
					itmp = fem->element[ii].kn[i];
					x[i] = fem->node[itmp].x;
					y[i] = fem->node[itmp].y;
					z[i] = fem->node[itmp].z;
				}
				for (i = 0; i < 4; i++) {
					k = i;
					l = (i + 1) % 4;
					m = (i + 2) % 4;
					n = (i + 3) % 4;
					a[k] = alpha[k] * (x[l] * (y[m] * z[n] - y[n] * z[m]) + x[m] * (y[n] * z[l] - y[l] * z[n]) + x[n] * (y[l] * z[m] - y[m] * z[l]));
					b[k] = alpha[k] * (y[l] * (z[n] - z[m]) + y[m] * (z[l] - z[n]) + y[n] * (z[m] - z[l]));
					c[k] = alpha[k] * (z[l] * (x[n] - x[m]) + z[m] * (x[l] - x[n]) + z[n] * (x[m] - x[l]));
					d[k] = alpha[k] * (x[l] * (y[n] - y[m]) + x[m] * (y[l] - y[n]) + x[n] * (y[m] - y[l]));
				}
				Ve = (a[0] + a[1] + a[2] + a[3]) / 6.0;
				for (i = 0; i < 6; i++) {
					aa[i] = x[edge[i][1]] - x[edge[i][0]];
					bb[i] = y[edge[i][1]] - y[edge[i][0]];
					cc[i] = z[edge[i][1]] - z[edge[i][0]];
					ll[i] = sqrt(aa[i] * aa[i] + bb[i] * bb[i] + cc[i] * cc[i]);
				}
				
				area[0] = sqrt(pow(bb[3] * cc[4] - bb[4] * cc[3], 2.0) 
							+ pow(cc[3] * aa[4] - cc[4] * aa[3], 2.0) + pow(aa[3] * bb[4] - aa[4] * bb[3], 2.0));
				area[1] = sqrt(pow(bb[1] * cc[2] - bb[2] * cc[1], 2.0) 
							+ pow(cc[1] * aa[2] - cc[2] * aa[1], 2.0) + pow(aa[1] * bb[2] - aa[2] * bb[1], 2.0));
				area[2] = sqrt(pow(bb[0] * cc[2] - bb[2] * cc[0], 2.0) 
							+ pow(cc[0] * aa[2] - cc[2] * aa[0], 2.0) + pow(aa[0] * bb[2] - aa[2] * bb[0], 2.0));
				area[3] = sqrt(pow(bb[0] * cc[1] - bb[1] * cc[0], 2.0) 
							+ pow(cc[0] * aa[1] - cc[1] * aa[0], 2.0) + pow(aa[0] * bb[1] - aa[1] * bb[0], 2.0));

				for (i = 0; i < 6; i++) {
					if (cc[i] > ZERO || (std::abs(cc[i]) <= ZERO && bb[i] > ZERO) 
							|| (std::abs(cc[i]) <= ZERO && std::abs(bb[i]) <= ZERO 
					&& aa[i] > ZERO)) {
						ss[i] = 1.0;
					} else {
						ss[i] = -1.0;
					}
				}
				
				L1 = (a[0] + b[0] * x0 + c[0] * y0 + d[0] * z0) / (6.0 * Ve);
				L2 = (a[1] + b[1] * x0 + c[1] * y0 + d[1] * z0) / (6.0 * Ve);
				L3 = (a[2] + b[2] * x0 + c[2] * y0 + d[2] * z0) / (6.0 * Ve);
				L4 = (a[3] + b[3] * x0 + c[3] * y0 + d[3] * z0) / (6.0 * Ve);

				if ((L1 >= -NEARLY_ZERO) && (L2 >= -NEARLY_ZERO) 
						&& (L3 >= -NEARLY_ZERO) && (L4 >= -NEARLY_ZERO)) {
					matID = fem->element[ii].matID;
					flag = 1;
								
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
					if (fem->order == 1) {

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
						
						for (i = 0; i < 6; i++) {
							U[i] *= ss[i] * ll[i];
							V[i] *= ss[i] * ll[i];
							W[i] *= ss[i] * ll[i];
						}
									
						//elementmatrix(not CURVE) より移植
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

						Vx[0] = L1x * L2y; Wx[0] = L1x * L2z;
						Vx[1] = L1x * L3y; Wx[1] = L1x * L3z;
						Vx[2] = L1x * L4y; Wx[2] = L1x * L4z;
						Vx[3] = L2x * L3y; Wx[3] = L2x * L3z;
						Vx[4] = L4x * L2y; Wx[4] = L4x * L2z;
						Vx[5] = L3x * L4y; Wx[5] = L3x * L4z;
						Vx[6] = L1y * L2x; Wx[6] = L1z * L2x;
						Vx[7] = L1y * L3x; Wx[7] = L1z * L3x;
						Vx[8] = L1y * L4x; Wx[8] = L1z * L4x;
						Vx[9] = L2y * L3x; Wx[9] = L2z * L3x;
						Vx[10] = L4y * L2x; Wx[10] = L4z * L2x;
						Vx[11] = L3y * L4x; Wx[11] = L3z * L4x;

						Uy[0] = L1y * L2x; Wy[0] = L1y * L2z;
						Uy[1] = L1y * L3x; Wy[1] = L1y * L3z;
						Uy[2] = L1y * L4x; Wy[2] = L1y * L4z;
						Uy[3] = L2y * L3x; Wy[3] = L2y * L3z;
						Uy[4] = L4y * L2x; Wy[4] = L4y * L2z;
						Uy[5] = L3y * L4x; Wy[5] = L3y * L4z;
						Uy[6] = L1x * L2y; Wy[6] = L1z * L2y;
						Uy[7] = L1x * L3y; Wy[7] = L1z * L3y;
						Uy[8] = L1x * L4y; Wy[8] = L1z * L4y;
						Uy[9] = L2x * L3y; Wy[9] = L2z * L3y;
						Uy[10] = L4x * L2y; Wy[10] = L4z * L2y;
						Uy[11] = L3x * L4y; Wy[11] = L3z * L4y;

						Uz[0] = L1z * L2x; Vz[0] = L1z * L2y;
						Uz[1] = L1z * L3x; Vz[1] = L1z * L3y;
						Uz[2] = L1z * L4x; Vz[2] = L1z * L4y;
						Uz[3] = L2z * L3x; Vz[3] = L2z * L3y;
						Uz[4] = L4z * L2x; Vz[4] = L4z * L2y;
						Uz[5] = L3z * L4x; Vz[5] = L3z * L4y;
						Uz[6] = L1x * L2z; Vz[6] = L1y * L2z;
						Uz[7] = L1x * L3z; Vz[7] = L1y * L3z;
						Uz[8] = L1x * L4z; Vz[8] = L1y * L4z;
						Uz[9] = L2x * L3z; Vz[9] = L2y * L3z;
						Uz[10] = L4x * L2z; Vz[10] = L4y * L2z;
						Uz[11] = L3x * L4z; Vz[11] = L3y * L4z;

						Vx[12] = (L3x * L4 + L3 * L4x) * L2y; Wx[12] = (L3x * L4 + L3 * L4x) * L2z;
						Vx[13] = (L4x * L2 + L4 * L2x) * L3y; Wx[13] = (L4x * L2 + L4 * L2x) * L3z;
						Vx[14] = (L2x * L3 + L2 * L3x) * L4y; Wx[14] = (L2x * L3 + L2 * L3x) * L4z;
						Vx[15] = (L4x * L3 + L4 * L3x) * L1y; Wx[15] = (L4x * L3 + L4 * L3x) * L1z;
						Vx[16] = (L3x * L1 + L3 * L1x) * L4y; Wx[16] = (L3x * L1 + L3 * L1x) * L4z;
						Vx[17] = (L1x * L4 + L1 * L4x) * L3y; Wx[17] = (L1x * L4 + L1 * L4x) * L3z;
						Vx[18] = (L1x * L2 + L1 * L2x) * L4y; Wx[18] = (L1x * L2 + L1 * L2x) * L4z;
						Vx[19] = (L2x * L4 + L2 * L4x) * L1y; Wx[19] = (L2x * L4 + L2 * L4x) * L1z;
						Vx[20] = (L4x * L1 + L4 * L1x) * L2y; Wx[20] = (L4x * L1 + L4 * L1x) * L2z;
						Vx[21] = (L2x * L1 + L2 * L1x) * L3y; Wx[21] = (L2x * L1 + L2 * L1x) * L3z;
						Vx[22] = (L1x * L3 + L1 * L3x) * L2y; Wx[22] = (L1x * L3 + L1 * L3x) * L2z;
						Vx[23] = (L3x * L2 + L3 * L2x) * L1y; Wx[23] = (L3x * L2 + L3 * L2x) * L1z;

						Uy[12] = (L3y * L4 + L3 * L4y) * L2x; Wy[12] = (L3y * L4 + L3 * L4y) * L2z;
						Uy[13] = (L4y * L2 + L4 * L2y) * L3x; Wy[13] = (L4y * L2 + L4 * L2y) * L3z;
						Uy[14] = (L2y * L3 + L2 * L3y) * L4x; Wy[14] = (L2y * L3 + L2 * L3y) * L4z;
						Uy[15] = (L4y * L3 + L4 * L3y) * L1x; Wy[15] = (L4y * L3 + L4 * L3y) * L1z;
						Uy[16] = (L3y * L1 + L3 * L1y) * L4x; Wy[16] = (L3y * L1 + L3 * L1y) * L4z;
						Uy[17] = (L1y * L4 + L1 * L4y) * L3x; Wy[17] = (L1y * L4 + L1 * L4y) * L3z;
						Uy[18] = (L1y * L2 + L1 * L2y) * L4x; Wy[18] = (L1y * L2 + L1 * L2y) * L4z;
						Uy[19] = (L2y * L4 + L2 * L4y) * L1x; Wy[19] = (L2y * L4 + L2 * L4y) * L1z;
						Uy[20] = (L4y * L1 + L4 * L1y) * L2x; Wy[20] = (L4y * L1 + L4 * L1y) * L2z;
						Uy[21] = (L2y * L1 + L2 * L1y) * L3x; Wy[21] = (L2y * L1 + L2 * L1y) * L3z;
						Uy[22] = (L1y * L3 + L1 * L3y) * L2x; Wy[22] = (L1y * L3 + L1 * L3y) * L2z;
						Uy[23] = (L3y * L2 + L3 * L2y) * L1x; Wy[23] = (L3y * L2 + L3 * L2y) * L1z;

						Uz[12] = (L3z * L4 + L3 * L4z) * L2x; Vz[12] = (L3z * L4 + L3 * L4z) * L2y;
						Uz[13] = (L4z * L2 + L4 * L2z) * L3x; Vz[13] = (L4z * L2 + L4 * L2z) * L3y;
						Uz[14] = (L2z * L3 + L2 * L3z) * L4x; Vz[14] = (L2z * L3 + L2 * L3z) * L4y;
						Uz[15] = (L4z * L3 + L4 * L3z) * L1x; Vz[15] = (L4z * L3 + L4 * L3z) * L1y;
						Uz[16] = (L3z * L1 + L3 * L1z) * L4x; Vz[16] = (L3z * L1 + L3 * L1z) * L4y;
						Uz[17] = (L1z * L4 + L1 * L4z) * L3x; Vz[17] = (L1z * L4 + L1 * L4z) * L3y;
						Uz[18] = (L1z * L2 + L1 * L2z) * L4x; Vz[18] = (L1z * L2 + L1 * L2z) * L4y;
						Uz[19] = (L2z * L4 + L2 * L4z) * L1x; Vz[19] = (L2z * L4 + L2 * L4z) * L1y;
						Uz[20] = (L4z * L1 + L4 * L1z) * L2x; Vz[20] = (L4z * L1 + L4 * L1z) * L2y;
						Uz[21] = (L2z * L1 + L2 * L1z) * L3x; Vz[21] = (L2z * L1 + L2 * L1z) * L3y;
						Uz[22] = (L1z * L3 + L1 * L3z) * L2x; Vz[22] = (L1z * L3 + L1 * L3z) * L2y;
						Uz[23] = (L3z * L2 + L3 * L2z) * L1x; Vz[23] = (L3z * L2 + L3 * L2z) * L1y;
									
						for (i = 0; i < 12; i++) {
							U[i] *= ll[i % 6];
							V[i] *= ll[i % 6];
							W[i] *= ll[i % 6];

							Uy[i] *= ll[i % 6];
							Uz[i] *= ll[i % 6];
							Vz[i] *= ll[i % 6];
							Vx[i] *= ll[i % 6];
							Wx[i] *= ll[i % 6];
							Wy[i] *= ll[i % 6];
						}
						for (i = 12; i < fem->n_en; i++) {
							U[i] *= co[i];
							V[i] *= co[i];
							W[i] *= co[i];

							Uy[i] *= co[i];
							Uz[i] *= co[i];
							Vz[i] *= co[i];
							Vx[i] *= co[i];
							Wx[i] *= co[i];
							Wy[i] *= co[i];
						}
					}
					for (i = 0; i < fem->n_en; i++) {
						itmp = fem->element[ii].kk[i];
						if (itmp >= 0) {
							vv[i] = fem->phi[mode * fem->n_edge + itmp];
						} else {
							vv[i] = 0.0;
						}
					}
					
					for (i = 0; i < fem->n_en; i++) {
						phix += U[i] * vv[i];
						phiy += V[i] * vv[i];
						phiz += W[i] * vv[i];
					}
								
					// added by satonori 20180614
					for (i = 0; i < fem->n_en; i++) {
						if (fem->unknown == 1) {
							phix2 += (Wy[i] - Vz[i]) * vv[i] / (-cj * k0 * ZZ);
							phiy2 += (Uz[i] - Wx[i]) * vv[i] / (-cj * k0 * ZZ);
							phiz2 += (Vx[i] - Uy[i]) * vv[i] / (-cj * k0 * ZZ);
						}
						else {
							phix2 += (Wy[i] - Vz[i]) * vv[i] / (cj * k0 / ZZ * data->par.er[matID]);
							phiy2 += (Uz[i] - Wx[i]) * vv[i] / (cj * k0 / ZZ * data->par.er[matID]);
							phiz2 += (Vx[i] - Uy[i]) * vv[i] / (cj * k0 / ZZ * data->par.er[matID]);
						}
					}
					// fprintf(fp, "%10.6e %10.6e %10.6e\n",
					// std::abs(phix), std::abs(phiy), std::abs(phiz));
					break;
				}
			} // loop end of iii

			if (flag == 0) {
				fprintf(stderr, "Can't find element in intpol\n");
				exit(EXIT_FAILURE);
			} else {
				if(data->mozaic.mozflag == 1){
					if(matID == data->mozaic.mozaicmatID 
						|| matID == data->mozaic.mozaicmatID+1){
						ii_pixel = data->mozaic.ii_pixel[ii];
						jj_pixel = data->mozaic.jj_pixel[ii];
						
						if(data->mozaic.BorW[kk][ii_pixel][jj_pixel] == 1){
							matID = data->mozaic.coreID;
						}else{
							matID = data->mozaic.cladID;
						}
						
					}
				}

				k = iv + ih*nv;
				
				mat[k].x = matID; 
				mat[k].y = data->par.er[matID]; 
				mat[k].z = sqrt(data->par.er[matID]);

				phi[k].x = phix2;
				phi[k].y = phiy2;
				phi[k].z = phiz2;
				
				if (nn == 2) {
					pict->absxz[iv][ih].x = std::abs(phix2);
					pict->absxz[iv][ih].y = std::abs(phiy2);
					pict->absxz[iv][ih].z = std::abs(phiz2);
					pict->realxz[iv][ih].x = std::real(phix2);
					pict->realxz[iv][ih].y = std::real(phiy2);
					pict->realxz[iv][ih].z = std::real(phiz2);
				}
			}
		}
	}
	fprintf(stderr, "finish to intpol\nNow, file output");
	
	
	if (mode == 0) {		
		if(nn == 0) sprintf(string, "%s-xy_structure.slv", data->name);
		if(nn == 1) sprintf(string, "%s-zy_structure.slv", data->name);
		if(nn == 2) sprintf(string, "%s-zx_structure.slv", data->name);
		if ((fp1 = fopen(string, "w")) == NULL) {
			fprintf(stderr, "can't open file ( %s )\n", string);
			exit(EXIT_FAILURE);
		}
	}

	if(nn == 0) sprintf(string, "%s-mode%lld-wl%1.3f-xyre.slv", data->name, mode, data->par.wavelength);
	else if(nn == 1) sprintf(string, "%s-mode%lld-wl%1.3f-zyre.slv", data->name, mode, data->par.wavelength);
	else if(nn == 2) sprintf(string, "%s-mode%lld-wl%1.3f-zxre.slv", data->name, mode, data->par.wavelength);
	if ((fp2 = fopen(string, "w")) == NULL) {
		fprintf(stderr, "can't open file ( %s )\n", string);
		exit(EXIT_FAILURE);
	}

	if(nn == 0) sprintf(string, "%s-mode%lld-wl%1.3f-xy.slv", data->name, mode, data->par.wavelength);
	else if(nn == 1) sprintf(string, "%s-mode%lld-wl%1.3f-zy.slv", data->name, mode, data->par.wavelength);
	else if(nn == 2) sprintf(string, "%s-mode%lld-wl%1.3f-zx.slv", data->name, mode, data->par.wavelength);
	if ((fp3 = fopen(string, "w")) == NULL) {
		fprintf(stderr, "can't open file ( %s )\n", string);
		exit(EXIT_FAILURE);
	}
	
	if (mode == 0) {		
		fprintf(fp1, "%lf\n", par->wavelength);
		fprintf(fp1, "%d\n", 1);
		fprintf(fp1, "%d\n", 11);
		fprintf(fp1, "%lld\n", par->n_material);
	}
	
	fprintf(fp2, "%lf\n", par->wavelength);
	fprintf(fp2, "%d\n", 1);
	fprintf(fp2, "%d\n", 11);
	fprintf(fp2, "%lld\n", par->n_material);
	
	fprintf(fp3, "%lf\n", par->wavelength);
	fprintf(fp3, "%d\n", 1);
	fprintf(fp3, "%d\n", 11);
	fprintf(fp3, "%lld\n", par->n_material);

	for (i = 0; i < par->n_material; i++) {
		if (mode == 0) fprintf(fp1, "%lf\n", par->er[i]);
		fprintf(fp2, "%lf\n", par->er[i]);
		fprintf(fp3, "%lf\n", par->er[i]);
	}

	if (mode == 0) fprintf(fp1, "%lf\n", par->er[0]);
	fprintf(fp2, "%lf\n", par->er[0]);
	fprintf(fp3, "%lf\n", par->er[0]);
	
	if (mode == 0) {
		fprintf(fp1, "%lf %lf %d\n", hmax, hmin, nh);
		fprintf(fp1, "%lf %lf %d\n", vmax, vmin, nv);
	}
	fprintf(fp2, "%lf %lf %d\n", hmax, hmin, nh);
	fprintf(fp2, "%lf %lf %d\n", vmax, vmin, nv);
	fprintf(fp3, "%lf %lf %d\n", hmax, hmin, nh);
	fprintf(fp3, "%lf %lf %d\n", vmax, vmin, nv);

	for (i = 0; i < nh*nv; i++) {
		// 材料のx成分にはmatIDが格納されている
		if (mode == 0) fprintf(fp1, "%.0lf %.6e %.6e\n", mat[i].x, mat[i].y, mat[i].z);
		fprintf(fp2, "%.6e %.6e %.6e\n", std::real(phi[i].x), std::real(phi[i].y), std::real(phi[i].z));
		fprintf(fp3, "%.6e %.6e %.6e\n",  std::abs(phi[i].x),  std::abs(phi[i].y),  std::abs(phi[i].z));
	}

	if (mode == 0) {		
		fprintf(fp1, "%d %d\n", 1, 1);
		fprintf(fp1, "%lf %lf\n", 0.0, 0.0);
	}
	
	fprintf(fp2, "%d %d\n", 1, 1);
	fprintf(fp2, "%lf %lf\n", 0.0, 0.0);
	
	fprintf(fp3, "%d %d\n", 1, 1);
	fprintf(fp3, "%lf %lf\n", 0.0, 0.0);
	
	if ((mode == 0 && fclose(fp1) == EOF) || fclose(fp2) == EOF || fclose(fp3) == EOF) {
		fprintf(stderr, "can't close file in intpol.cpp\n");
		exit(EXIT_FAILURE);
	}
	
	if (nn == 2)  WritePPMFile(data, mode);
	
	fprintf(stderr, "intpol END\n");
	free(mat);
	free(phi);
	return;
}

void addCoreframe (DataTable *data) {
	int ii, jj;
	long long int i, ix, iz, count, tmpmatID, ele;
	long long int **flag;
	double er_core = 0;
	FEM *fem = &(data->fem);
	Picture *pict = &(data->pict);

	fprintf(stderr, "addCoreframe Start\n");
	for (i = 0; i < data->par.n_material; i++) {
		if (er_core < data->par.er[i]) er_core = data->par.er[i];
	}

	
	for (ix = 0; ix < pict->nx; ix++) {
		for (iz = 0; iz < pict->nz; iz++) pict->f_flag[ix][iz] = -1;
	}
	for (ix = 0; ix < pict->nx; ix++) {
		for (iz = 0; iz < pict->nz; iz++) {
			ele = pict->element[ix][iz];
			tmpmatID = fem->element[ele].matID;

			if (data->mozaic.mozflag == 1) {
				if (tmpmatID == data->mozaic.mozaicmatID
						|| tmpmatID == data->mozaic.mozaicmatID+1) {
					ii = data->mozaic.ii_pixel[ele];
					jj = data->mozaic.jj_pixel[ele];

					if(data->mozaic.BorW[0][ii][jj] == 1){
						tmpmatID = data->mozaic.coreID;
					} else {
						tmpmatID = data->mozaic.cladID;
					}
				}
			}
			
			if (fabs(data->par.er[tmpmatID] - er_core) < 1e-6) pict->f_flag[ix][iz] = 0;
			else pict->f_flag[ix][iz] = 1;
		}
	}
	for (ix = 0; ix < pict->nx; ix++) {
		for (iz = 0; iz < pict->nz; iz++) {
			if (pict->f_flag[ix][iz] != 0 && pict->f_flag[ix][iz] != 1) {
				fprintf(stderr, "Error in addCoreframe\n");
				exit(EXIT_FAILURE);
			}
		}
	}
	
	ALLOCATION(flag, long long int *, pict->nx)
	for (ix = 0; ix < pict->nx; ix++) {
		ALLOCATION(flag[ix], long long int, pict->nz)
		for (iz = 0; iz < pict->nz; iz++) flag[ix][iz] = -1;
	}
	for (ix = 0; ix < pict->nx; ix++) {
		for (iz = 0; iz < pict->nz; iz++) {
			tmpmatID = 0;
			count = 0;
			if (pict->f_flag[ix][iz] == 1) {//クラッドのピクセルが対象
				tmpmatID += pict->f_flag[ix][iz];
				count++;
				if (ix != 0) {
					tmpmatID += pict->f_flag[ix-1][iz];
					count++;
				}
				if (ix != pict->nx-1) {
					tmpmatID += pict->f_flag[ix+1][iz];
					count++;
				}
				if (iz != 0) {
					tmpmatID += pict->f_flag[ix][iz-1];
					count++;
				}
				if (iz != pict->nz-1) {
					tmpmatID += pict->f_flag[ix][iz+1];
					count++;
				}
				if (ix != 0 && iz != 0) {
					tmpmatID += pict->f_flag[ix-1][iz-1];
					count++;
				}
				if (ix != 0 && iz != pict->nz-1) {
					tmpmatID += pict->f_flag[ix-1][iz+1];
					count++;
				}
				if (ix != pict->nx-1 && iz != 0) {
					tmpmatID += pict->f_flag[ix+1][iz-1];
					count++;
				}
				if (ix != pict->nx-1 && iz != pict->nz-1) {
					tmpmatID += pict->f_flag[ix+1][iz+1];
					count++;
				}
				if (count != tmpmatID) flag[ix][iz] = 2;//周りにコアがあればフレームにする
			}
		}
	}
	
	for (ix = 0; ix < pict->nx; ix++) {//フレームを太くしたいだけ
		for (iz = 0; iz < pict->nz; iz++) {
			if (flag[ix][iz] == 2) {
				if (ix != pict->nx-1) pict->f_flag[ix+1][iz] = 2;
				if (iz != pict->nz-1) pict->f_flag[ix][iz+1] = 2;
			}
		}
	}
	for (ix = 0; ix < pict->nx; ix++) {
		for (iz = 0; iz < pict->nz; iz++) {
			if (flag[ix][iz] == 2) pict->f_flag[ix][iz] = 2;
		}
		free(flag[ix]);
	}
	free(flag);

	fprintf(stderr, "addCoreframe END\n");
	return;
}

void WritePPMFile (DataTable *data, long long int mode) {
	long long int i, j, k, l, m, n;
	long long int num;
	double absval, absmax = 0.0;
	double realval, realmax = 0.0;
	double imagval, imagmax = 0.0;
	char tmp1[256], tmp2[256], tmp3[256], tmp4[256], cmd[256];
	char *component[3] = {"x", "y", "z"};
	FILE *fp1, *fp2, *fp3;
	Picture *pict = &(data->pict);
	RGBcolor *rgb = pict->rgb;			// Hot Color
	RGBcolor *rgb2 = pict->rgb2;		// Jet Color
	XYZdouble div;

	div.x = (data->fem.max.x - data->fem.min.x) / pict->nx;
	div.y = (data->fem.max.y - data->fem.min.y) / pict->ny;
	div.z = (data->fem.max.z - data->fem.min.z) / pict->nz;

	for (i = 0; i < pict->nx; i++) {
		for (j = 0; j < pict->nz; j++) {
			if (data->pict.modifiedScale) {
				if ((i+0.5) * div.x < data->pml.dx1*1.5 || (i+0.5) * div.x > (data->fem.max.x - data->fem.min.x) - data->pml.dx2*1.5
					|| (j+0.5) * div.z < data->pml.dz1*1.5 || (j+0.5) * div.z > (data->fem.max.z - data->fem.min.z) - data->pml.dz2*1.5)
					continue;
			}

			if (absmax < pict->absxz[i][j].x) absmax = pict->absxz[i][j].x;
			if (absmax < pict->absxz[i][j].y) absmax = pict->absxz[i][j].y;
			if (absmax < pict->absxz[i][j].z) absmax = pict->absxz[i][j].z;
			if (realmax < fabs(pict->realxz[i][j].x)) realmax = fabs(pict->realxz[i][j].x);
			if (realmax < fabs(pict->realxz[i][j].y)) realmax = fabs(pict->realxz[i][j].y);
			if (realmax < fabs(pict->realxz[i][j].z)) realmax = fabs(pict->realxz[i][j].z);
		}
	}
	for (num = 0; num < 3; num++) {
		sprintf(tmp1, "./abs%lld.ppm", data->par.myid);
		if ((fp1 = fopen(tmp1, "w")) == NULL) {
			fprintf(stderr, "cannot open file (%s)\n", tmp1);
			exit(1);
		}
		sprintf(tmp2, "./real%lld.ppm", data->par.myid);
		if ((fp2 = fopen(tmp2, "w")) == NULL) {
			fprintf(stderr, "cannot open file (%s)\n", tmp2);
			exit(1);
		}
		fprintf(fp1, "P3\n");
		fprintf(fp2, "P3\n");
		fprintf(fp1, "%lld %lld\n", pict->nz, pict->nx);
		fprintf(fp2, "%lld %lld\n", pict->nz, pict->nx);
		fprintf(fp1, "%lld\n", PPMCOLOR - 1);
		fprintf(fp2, "%lld\n", PPMCOLOR - 1);
		for (i = 0; i < pict->nx; i++) {
			for (j = 0; j < pict->nz; j++) {
				if (num == 0) {
					absval = pict->absxz[i][j].x;
					realval = pict->realxz[i][j].x;
				} else if (num == 1) {
					absval = pict->absxz[i][j].y;
					realval = pict->realxz[i][j].y;
				} else if (num == 2) {
					absval = pict->absxz[i][j].z;
					realval = pict->realxz[i][j].z;
				}
				k = (long long int)((PPMCOLOR - 1) * (absval / absmax));
				l = (long long int)((PPMCOLOR - 1) * (realval + realmax)/(2.0*realmax));
				if (k >= PPMCOLOR) k = PPMCOLOR - 1;
				if (k < 0) k = 0;
				if (l >= PPMCOLOR) l = PPMCOLOR - 1;
				if (l < 0) l = 0;
				if (pict->f_flag[i][j] == 2) {
					fprintf(fp1, "%d %d %d\n", PPMCOLOR-1, PPMCOLOR-1, PPMCOLOR-1);
					fprintf(fp2, "0 0 0\n");
				}
				else if ((i+0.5) * div.x < data->pml.dx1 || (i+0.5) * div.x > (data->fem.max.x - data->fem.min.x) - data->pml.dx2
						|| (j+0.5) * div.z < data->pml.dz1 || (j+0.5) * div.z > (data->fem.max.z - data->fem.min.z) - data->pml.dz2) {
					// PML領域は着色する。どんな色でも構わないが、軽ーく色を付けた
					fprintf(fp1, "%lld %lld %lld\n", rgb[k].b,  rgb[k].g,  rgb[k].r);
					fprintf(fp2, "%lld %lld %lld\n", rgb2[l].r, rgb2[l].g*3/5, rgb2[l].b);
				} 
				else if (data->pict.modifiedScale && ((i+0.5) * div.x < data->pml.dx1*1.5 || (i+0.5) * div.x > (data->fem.max.x - data->fem.min.x) - data->pml.dx2*1.5
						|| (j+0.5) * div.z < data->pml.dz1*1.5 || (j+0.5) * div.z > (data->fem.max.z - data->fem.min.z) - data->pml.dz2*1.5)) {
					// 修正規格化での対象外領域(PML付近)は着色する。どんな色でも構わないが、軽ーく色を付けた
					fprintf(fp1, "%lld %lld %lld\n", rgb[k].b,  rgb[k].g,  rgb[k].r);
					fprintf(fp2, "%lld %lld %lld\n", rgb2[l].r, rgb2[l].g*4/5, rgb2[l].b);
				} 
				else {
					fprintf(fp1, "%lld %lld %lld\n",  rgb[k].r,  rgb[k].g,  rgb[k].b);
					fprintf(fp2, "%lld %lld %lld\n", rgb2[l].r, rgb2[l].g, rgb2[l].b);
				}
			}
		}
		fclose(fp1);
		fclose(fp2);
		sprintf(cmd, "/usr/bin/convert ./abs%lld.ppm ./Field/abs-%.3lf/mode%lld-%s.png", 
				data->par.myid, data->par.wavelength, mode, component[num]);
		system(cmd);
		sprintf(cmd, "/usr/bin/convert ./real%lld.ppm ./Field/real-%.3lf/mode%lld-%s.png", 
				data->par.myid, data->par.wavelength, mode, component[num]);
		system(cmd);
	}
	return;
}

void setRGBTable (DataTable *data) {
	long long int i, j, k;
	double value;
	Picture *pict = &(data->pict);
	RGBcolor *rgb = pict->rgb;
	RGBcolor *rgb2 = pict->rgb2;

	for (i = 0; i < PPMCOLOR; i++) {
		value = (double)i / (PPMCOLOR-1);
		rgb[i].r = 0; rgb[i].g = 0; rgb[i].b = 0;
		if(value < 1.0/3.0){
		    rgb[i].r = (value - 0.0/3.0)*3.0*(PPMCOLOR-1);
		}else if(value < 2.0/3.0){
		    rgb[i].r = PPMCOLOR-1;
		    rgb[i].g = (value - 1.0/3.0)*3.0*(PPMCOLOR-1);
		}else{
		    rgb[i].r = PPMCOLOR-1;
			rgb[i].g = PPMCOLOR-1;
		    rgb[i].b = (value - 2.0/3.0)*3.0*(PPMCOLOR-1);
		}
		if (rgb[i].r < 0) rgb[i].r = 1;
		if (rgb[i].g < 0) rgb[i].g = 1;
		if (rgb[i].b < 0) rgb[i].b = 1;
		if (rgb[i].r > (PPMCOLOR-1)) rgb[i].r = (PPMCOLOR-1);
		if (rgb[i].g > (PPMCOLOR-1)) rgb[i].g = (PPMCOLOR-1);
		if (rgb[i].b > (PPMCOLOR-1)) rgb[i].b = (PPMCOLOR-1);
	}
	for (i = 0; i < PPMCOLOR; i++) {
		value = (double)i / (PPMCOLOR-1);
		rgb2[i].r = 0; rgb2[i].g = 0; rgb2[i].b = 0;
		if(value < 1.0/8.0){ // 0/8 ~ 1/8
		    rgb2[i].r = 0;
		    rgb2[i].g = 0;
		    rgb2[i].b = (int)((PPMCOLOR-1)*((value+1.0/8.0)*4.0)); //線形増加
		}else if(value < 3.0/8.0){ // 1/8 ~ 3/8
		    rgb2[i].r = 0;
		    rgb2[i].g = (int)((PPMCOLOR-1)*((value-1.0/8.0)*4.0));
		    rgb2[i].b = PPMCOLOR-1;
		}else if(value < 5.0/8.0){ // 3/8 ~ 5/8
		    rgb2[i].r = (int)((PPMCOLOR-1)*((value-3.0/8.0)*4.0));
		    rgb2[i].g = PPMCOLOR-1;
		    rgb2[i].b = (int)((PPMCOLOR-1)*(1.0-(value-3.0/8.0)*4.0)); //線形減少
		}else if(value < 7.0/8.0){ // 5/8 ~ 7/8
		    rgb2[i].r = PPMCOLOR-1;
		    rgb2[i].g = (int)((PPMCOLOR-1)*(1.0-(value-5.0/8.0)*4.0));
		    rgb2[i].b = 0;
		}else{ // 7/8 ~ 8/8
		    rgb2[i].r = (int)((PPMCOLOR-1)*(1.0-(value-7.0/8.0)*4.0));
		    rgb2[i].g = 0;
		    rgb2[i].b = 0;
		}
		if (rgb2[i].r < 0) rgb2[i].r = 1;
		if (rgb2[i].g < 0) rgb2[i].g = 1;
		if (rgb2[i].b < 0) rgb2[i].b = 1;
		if (rgb2[i].r > (PPMCOLOR-1)) rgb2[i].r = (PPMCOLOR-1);
		if (rgb2[i].g > (PPMCOLOR-1)) rgb2[i].g = (PPMCOLOR-1);
		if (rgb2[i].b > (PPMCOLOR-1)) rgb2[i].b = (PPMCOLOR-1);
	}
	return;
}

