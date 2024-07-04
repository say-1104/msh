#include "fem3d.h"

// added by Fujisawa
// difference with port.cpp is that this file is created for lamda sweeping
// data->calcflag == 2 is assumed (InputField-%d-%d-%1.3f is used)

#define OUT_2D_SURF_SLV_DIV_X 200
#define OUT_2D_SURF_SLV_DIV_Y 200

#define TWO_DIM_SOLVER_PATH "/home/muratsubaki/Solver/2D_VFEM-feast/edge"

#define EVEN

void createPortInfo(DataTable *data) {
	long long int i, j, k, flag;
	double z_grav;
	long long int ip;
	// long long int im;
	long long int i1, i2, i3;
	
	FEM *fem = &(data->fem);
	Param *par = &(data->par);

	long long int ii;

	if(data->calcflag != 2){
		// if data->calcflag == 2, the port data is read in inputData.cpp
		fem->n_port = 0; // par->PMLflag[i]の最大値
		for (i = 0; i < par->n_material; i++) {
			if (par->PMLflag[i] != 0) {
				if (par->PMLflag[i] > fem->n_port) {
					fem->n_port = par->PMLflag[i];
				}
			}
		}

		ALLOCATION(data->port, PortInfo, fem->n_port);
	} else if(data->calcflag == 2) {
		fprintf(stderr, "number of ports read from .param is %lld\n", fem->n_port);
	}

	// port0はPMLflag=1
	// port0が入射ポートでportangle+180してる（なぜ
	for (i = 0; i < par->n_material; i++) {
		ip = par->PMLflag[i];
		if (ip != 0) {
			data->port[ip - 1].x = par->PMLx0[i];
			data->port[ip - 1].y = 0.0;
			data->port[ip - 1].z = par->PMLz0[i];
			fprintf(stderr, "data->port[ip-1].x = %lf\n", data->port[ip-1].x);
			fprintf(stderr, "data->port[ip-1].y = %lf\n", data->port[ip-1].y);
			fprintf(stderr, "data->port[ip-1].z = %lf\n", data->port[ip-1].z);

			if (ip == 1) {
				/* kaki found error */
				data->port[ip - 1].angle = par->PMLangle[i] + 180.0;
			} else {
				data->port[ip - 1].angle = par->PMLangle[i];
			}
		}
	}
	par->INPUT_Z = data->port[0].z;
	fprintf(stderr, "number of ports is %lld\n", fem->n_port);
	
	// これらはモード数が決まってから．InputConditionで
	// ALLOCATION(fem->phi, std::complex<double>, fem->n_edge)
	// ALLOCATION(fem->phi_in, std::complex<double>, fem->n_edge)
	// for (i = 0; i < fem->n_edge; i++)
	// 	fem->phi[i] = 0.0;
	
	for (i = 0; i < fem->n_port; i++) {
		eachPort(data, &(data->port[i]), i);
		// ここでメッシュ作成後，savePortFile(2DVFEMよびだし), 
		// createS2d, normalizeInputPortPowerまで
		// port->m_func[i][j]
		// port->m_func2[i][j]
		// port->m_func_in[i][j]をアロケーションするが
		// 何に使っているか不明
	}
	
	

	for (i = 0; i < fem->ne; i++) {
		z_grav = 0.0;
		for (j = 0; j < 4; j++) {
			k = fem->element[i].kn[j];
			z_grav += fem->node[k].z;
		}
		z_grav /= 4.0;
		if (z_grav < par->INPUT_Z) {
			flag = 1;
		} else {
			flag = 2;
		}
		for (j = 0; j < fem->n_en; j++) {
			k = fem->element[i].kk[j];
			fem->inFlag[k] = flag | fem->inFlag[k];
		}
	}
	
	i1 = 0;
	i2 = 0;
	i3 = 0;
	for (i = 0; i < fem->n_edge; i++) {
		if (fem->inFlag[i] == 1)
			i1++;
		if (fem->inFlag[i] == 2)
			i2++;
		if (fem->inFlag[i] == 3)
			i3++;
	}
	
	return;
}

/*
 各ポートごとの情報を作る
 */
void eachPort(DataTable *data, PortInfo *port, long long int n_port) {
	FILE *fp;
	char string[256];
	long long int i, j, k, it, ii;
	long long int flag[4], f_count;
	long long int *e_flag = NULL, e_count;
	long long int *n1_flag = NULL, n1_count;
	long long int *n2_flag = NULL, n2_count;
	long long int n_count;
	long long int n1d_count, n2d_count;
	double buf_r, buf_i;
	std::complex<double> cj = std::complex<double>(0, 1);
	/* 各面内の節点番号 */
	/* mod by kaki- */
	static long long int tn[4][6];
	static long long int tn_a[4][6] = { { 1, 3, 2, 8, 9, 5 }, 
						{ 2, 3, 0, 9, 7, 6 }, 
						{ 3, 1, 0, 8, 4, 7 }, 
						{ 0, 1, 2, 4, 5, 6 } };
	static long long int tn_b[4][6] = { { 1, 2, 3, 5, 9, 8 }, 
						{ 2, 0, 3, 6, 7, 9 }, 
						{ 3, 0, 1, 7, 4, 8 }, 
						{ 0, 2, 1, 6, 5, 4 } };
	static long long int tn_c[4][6] = { { 1, 2, 3, 5, 9, 8 }, 
						{ 0, 3, 2, 7, 9, 6 }, 
						{ 3, 0, 1, 7, 4, 8 }, 
						{ 2, 1, 0, 5, 4, 6 } };
	/* 各面内の辺番号(CT/LN 要素) */
	static long long int te[4][3];
	static long long int te_a[4][3] = { { 4, 5, 3 }, { 5, 2, 1 }, 
						{ 4, 0, 2 }, { 0, 3, 1 } };
	static long long int te_b[4][3] = { { 3, 5, 4 }, { 1, 2, 5 }, 
						{ 2, 0, 4 }, { 1, 3, 0 } };
	static long long int te_c[4][3] = { { 3, 5, 4 }, { 2, 5, 1 }, 
						{ 2, 0, 4 }, { 3, 0, 1 } };
	/* 各面内の辺番号(LT/LN or LT/QN 要素) */
	static long long int te2[4][9];
	static long long int te2_a[4][9] = { { 10, 11, 9, 4, 5, 3, 13, 12, 14 }, 
						 { 5, 8, 1, 11, 2, 7, 15, 17, 16 }, 
						 { 4, 6, 2, 10, 0, 8, 19, 18, 20 }, 
						 { 0, 3, 7, 6, 9, 1, 21, 23, 22 } };
	static long long int te2_b[4][9] = { { 3, 5, 4, 9, 11, 10, 14, 12, 13 }, 
						 { 7, 2, 11, 1, 8, 5, 16, 17, 15 }, 
						 { 8, 0, 10, 2, 6, 4, 20, 18, 19 }, 
						 { 1, 9, 6, 7, 3, 0, 22, 23, 21 } };
	static long long int te2_c[4][9] = { { 3, 5, 4, 9, 11, 10, 14, 12, 13 }, 
						 { 2, 11, 7, 8, 5, 1, 17, 15, 16 }, 
						 { 8, 0, 10, 2, 6, 4, 20, 18, 19 }, 
						 { 9, 6, 1, 3, 0, 7, 23, 21, 22 } };
	// double zg;
	long long int imin, iele, n_tmp;
	long long int *div_flag = NULL;
	double x1, x2, z1, z2;
	long long int im;
	// long long int ip;
	long long int l_edge = 8, l_node = 6;
	
	long long int kaki;
	
	FEM *fem = &(data->fem);
	Param *par = &(data->par);
	// double max_x = 0.0, min_x = -6.0, max_z = 0.0, min_z = 1000.0;
	
	fprintf(stderr, "Port[%d]:\n", n_port);
	fprintf(stderr, "max_x = %lf min_x = %lf\n", 
		data->port_xmax[n_port], data->port_xmin[n_port]);
	fprintf(stderr, "max_z = %lf min_z = %lf\n", 
		data->port_zmax[n_port], data->port_zmin[n_port]);
	
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 6; j++) {
			if (n_port != 0) {
				tn[i][j] = tn_a[i][j];
			} else {
				tn[i][j] = tn_b[i][j];
			}
			tn[i][j] = tn_c[i][j];
		}
	}
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 3; j++) {
			if (n_port != 0) {
				te[i][j] = te_a[i][j];
			} else {
				te[i][j] = te_b[i][j];
			}
			te[i][j] = te_c[i][j];
		}
	}
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 9; j++) {
			if (n_port != 0) {
				te2[i][j] = te2_a[i][j];
			} else {
				te2[i][j] = te2_b[i][j];
			}
			te2[i][j] = te2_c[i][j];
		}
	}

	/* 要素次数に対応した要素内未知変数の数の設定 */
	switch (fem->order) {
		case 1:
			l_edge = 3;
			l_node = 3;
			break;
		case 2:
			l_edge = 6;
			l_node = 6;
			break;
		case 3:
			l_edge = 8;
			l_node = 6;
			break;
	}

	ALLOCATION(e_flag, long long int, fem->ne);
	ALLOCATION(n1_flag, long long int, fem->np);
	ALLOCATION(n2_flag, long long int, fem->n_edge);
	
	/* 面の法線ベクトル */
	x1 = sin(port->angle * M_PI / 180.0);
	z1 = cos(port->angle * M_PI / 180.0);
	
	e_count = 0;
	/* 各要素の各面がポート上にあるかを調べる */
	for (k = 0; k < fem->ne; k++) {
		f_count = 0;
		for (i = 0; i < 4; i++) {
			it = fem->element[k].kn[i];
			x2 = fem->node[it].x - port->x;
			z2 = fem->node[it].z - port->z;
			/* 導波路中心から面の中心へのベクトルが
				面の法線ベクトルと直交するかで判断 */
			if (std::abs(x1 * x2 + z1 * z2) < 1.0e-6) {
				//Port Division
				if(fem->node[it].x    > data->port_xmax[n_port]
					|| fem->node[it].x < data->port_xmin[n_port]
					|| fem->node[it].z > data->port_zmax[n_port]
					|| fem->node[it].z < data->port_zmin[n_port]){
					//f_count = -100;
					flag[i] = 0;
				}else{
					f_count++;
					flag[i] = 1;
				}
			} else {
				flag[i] = 0;
			}
		}
		im = fem->element[k].matID;
		/* 面がポート上であり PML との境界であるとき */
		if (f_count == 3 && ((par->PMLflag[im] == 0 && n_port != 0) 
				|| (par->PMLflag[im] != 0 && n_port == 0))) {
			e_count++;
			for (i = 0; i < 4; i++) {
				if (flag[i] == 0)
					e_flag[k] = i;
			}
		} else {
			e_flag[k] = -1;
		}
	}

	/* ポート上の２次元有限要素法データを構築する */
	ALLOCATION(port->element, Element2D, e_count);
	port->ne = e_count;
	e_count = 0;
	for (k = 0; k < fem->ne; k++) {
		/* 要素がポート上であるとき */
		if (e_flag[k] != -1) {
			if (fem->order == 1) {
				/* CT/LN 要素 */
				ii = e_flag[k];
				for (j = 0; j < 3; j++) {
					/* ３次元メッシュの節点番号を入れる */
					port->element[e_count].kn[j] = fem->element[k].kn[tn[ii][j]];
					port->element[e_count].ke[j] = fem->element[k].kk[te[ii][j]];
				}
				port->element[e_count].matID = fem->element[k].matID;
			} else if (fem->order == 2 || fem->order == 3) {
				/* LT/LN or LT/QN 要素 */
				ii = e_flag[k];
				for (j = 0; j < 6; j++) {
					port->element[e_count].kn[j] = fem->element[k].kn[tn[ii][j]];
					port->element[e_count].ke[j] = fem->element[k].kk[te2[ii][j]];
				}
				if (fem->order == 3) {
					if (fem->element[k].kk[te2[ii][6]] < 0) {
					kaki = 1;
					} else if (fem->element[k].kk[te2[ii][7]] < 0) {
					kaki = 2;
					} else {
					kaki = 0;
					}
					for (j = 0; j < 3; j++) {
					port->element[e_count].ke[j + 6] 
						= fem->element[k].kk[te2[ii][(j + kaki) % 3 + 6]];
					}
				}
				port->element[e_count].matID = fem->element[k].matID;
			}
			e_count++;
		}
	}

	/* ポート上の節点番号を新たに定義する */
	for (i = 0; i < fem->np; i++) n1_flag[i] = -1;
	for (i = 0; i < fem->n_edge; i++) n2_flag[i] = -1;
	for (k = 0; k < port->ne; k++) {
		for (i = 0; i < l_node; i++) {
			n1_flag[port->element[k].kn[i]]++;
		}
		for (i = 0; i < l_edge; i++) {
			n2_flag[port->element[k].ke[i]]++;
		}
	}
	
	/* for diriclet */
	for (k = 0; k < port->ne; k++) {
		for (i = 0; i < 3; i++) {
			if (n2_flag[port->element[k].ke[i]] != -1 && port->element[k].ke[i]
				>= fem->nr_edge) {
				if (fem->order == 1) {
					n2_flag[port->element[k].ke[i]] = -2;
					n1_flag[port->element[k].kn[i]] = -2;
					n1_flag[port->element[k].kn[(i + 1) % 3]] = -2;
				} else if (fem->order == 2 || fem->order == 3) {
					n2_flag[port->element[k].ke[i]] = -2;
					n2_flag[port->element[k].ke[i + 3]] = -2;
					n1_flag[port->element[k].kn[i]] = -2;
					n1_flag[port->element[k].kn[(i + 1) % 3]] = -2;
					n1_flag[port->element[k].kn[i + 3]] = -2;
				}
			}
		}
	}
	
	/* ディリクレ以外の節点数のカウント */
	n1_count = 0;
	for (i = 0; i < fem->np; i++) {
		if (n1_flag[i] != -1 && n1_flag[i] != -2) {
			n1_flag[i] = (fem->np + fem->n_edge + 1);
			n1_count++;
		}
	}
	/* ディリクレ以外の辺数のカウント */
	n2_count = 0;
	for (i = 0; i < fem->n_edge; i++) {
		if (n2_flag[i] != -1 && n2_flag[i] != -2) {
			n2_flag[i] = (fem->np + fem->n_edge + 1);
			n2_count++;
		}
	}
	
	ALLOCATION(div_flag, long long int, e_count);
	for (k = 0; k < e_count; k++)
	div_flag[k] = 0;
	
	/* 行列をバンド化するための節点番号付けバンド化は必須ではない */
	n_count = 0;
	for (k = 0; k < e_count; k++) {
		imin = (fem->n_edge + fem->np);
		iele = 0;
		for (i = 0; i < e_count; i++) {
			if (div_flag[i] == 0) {
				for (j = 0; j < l_node; j++) {
					n_tmp = n1_flag[port->element[i].kn[j]];
					if (n_tmp != -1 && n_tmp != -2 && imin > n_tmp) {
					imin = n_tmp;
					iele = i;
					}
				}
				for (j = 0; j < l_edge; j++) {
					n_tmp = n2_flag[port->element[i].ke[j]];
					if (n_tmp != -1 && n_tmp != -2 && imin > n_tmp) {
					imin = n_tmp;
					iele = i;
					}
				}
			}
		}

		div_flag[iele] = 1;
		for (i = 0; i < l_node; i++) {
			if (n1_flag[port->element[iele].kn[i]] == (fem->np + fem->n_edge + 1))
				n1_flag[port->element[iele].kn[i]] = n_count++;
		}
		for (i = 0; i < l_edge; i++) {
			if (n2_flag[port->element[iele].ke[i]] == (fem->np + fem->n_edge + 1))
				n2_flag[port->element[iele].ke[i]] = n_count++;
		}
	}
	free(div_flag);
	
	/* for diriclet */
	n1d_count = 0;
	for (i = 0; i < fem->np; i++) {
		if (n1_flag[i] == -2) {
			n1_flag[i] = n_count++;
			n1_count++;
			n1d_count++;
		}
	}
	n2d_count = 0;
	for (i = 0; i < fem->n_edge; i++) {
		if (n2_flag[i] == -2) {
			n2_flag[i] = n_count++;
			n2_count++;
			n2d_count++;
		}
	}
	/* ここまでで，全ての節点，辺の番号が決まる */
	
	port->np = n1_count;
	port->n_edge = n2_count;
	port->nr = port->np + port->n_edge - n1d_count - n2d_count;
	
	port->toLn = n1_flag;
	port->toLe = n2_flag;

	fprintf(stderr, "%d %d %d %d\n", n_port, port->np, port->n_edge, port->nr);

	savePortFile(data, port, n_port); // メッシュを出力し2DVFEM呼び出し
	createS2d(data, port, n_port); // まじでなぞ
	
	free(e_flag);

	if(data->calcflag == 0){
		for(i = 0 ; i < 10; i++){
			sprintf(string, "_input_port%lld-v%d.slv", n_port, i);
			if ((fp = fopen(string, "r")) == NULL) break;
			fclose(fp);
		}
		port->modenum = i;
		if(port->modenum == 0){
			fprintf(stderr, "input field cannot be read!!\n");
			exit(EXIT_FAILURE);
		}
		if(port->modenum == 20){
			fprintf(stderr, "too many input fields are created!!\n");
			exit(EXIT_FAILURE);
		}
		fprintf(stderr, " %d input fields are created for port %d\n", 
			port->modenum, n_port);
	}
	
	ALLOCATION(port->m_func, std::complex<double> *, port->modenum);
	ALLOCATION(port->beta, std::complex<double>, port->modenum);
	for (i = 0; i < port->modenum; i++) {
		ALLOCATION(port->m_func[i], std::complex<double>, port->n_edge+port->np);
		for (j = 0; j < port->n_edge + port->np; j++) {
			port->m_func[i][j] = 0.0;
		}
	}

	/* input & output field reading (vector field) */
	for (ii = 0; ii < port->modenum; ii++) {
		if(data->calcflag != 2){
			sprintf(string, "_input_port%lld-v%d.slv", n_port, ii);
		} else if(data->calcflag == 2) {
			if (data->par.numStructures != 1) {
				sprintf(string, "InputField_%d-%lld-%d-%1.3f", 
					data->par.strID, n_port, ii, data->par.wavelength);
			}
			else {
				sprintf(string, "InputField-%lld-%d-%1.3f", 
					n_port, ii, data->par.wavelength);
			}
		}

		if ((fp = fopen(string, "r")) == NULL) {
			fprintf(stderr, "can't open file (%s).\n", string);
			exit(EXIT_FAILURE);
		}
			
		// VFEMの仕様変更より
		fscanf(fp, "%lf %lf", &buf_r, &buf_i);
		port->beta[ii] = sqrt(buf_r + cj * buf_i);
		fprintf(stderr, "beta[%d] = complex(%lf, %lf)\n", ii, 
			real(port->beta[ii]), imag(port->beta[ii]));
		for (i = 0; i < port->nr; i++) {
			fscanf(fp, "%lf %lf", &buf_r, &buf_i);
			port->m_func[ii][i] = buf_r + cj * buf_i;
			// fprintf(stderr, "phi[%d] = complex(%lf, %lf)\n", i, buf_r, buf_i);
		}
		
		for (i = port->nr; i < port->n_edge + port->np; i++) {
			port->m_func[ii][i] = 0.0;
		}
		if (fclose(fp) == EOF) {
			fprintf(stderr, "can't close file (%s).\n", string);
			exit(EXIT_FAILURE);
		}
	}
	
	for (ii = 0; ii < port->modenum; ii++) {
		normalizeInputPortPower(data, port, n_port, ii);
	}
	
	return;
}


void savePortFile(DataTable *data, PortInfo *port, long long int n_port) {
	long long int i, j;
	long long int count;
	FILE *fp;
	char string[512] = "\0";
	double xx, zz;
	double cosA, sinA;
	
	FEM *fem = &(data->fem);
	Param *par = &(data->par);
	
	// fprintf(stderr, "angle = %lf\n", port->angle);
	
	fprintf(stderr, "angle = %lf\n", port->angle);
	
	cosA = cos(port->angle * M_PI / 180.0);
	sinA = sin(port->angle * M_PI / 180.0);
	
	sprintf(string, "_input_port%lld.msh", n_port);
	if ((fp = fopen(string, "w")) == NULL) {
		fprintf(stderr, "can't open file (%s).\n", string);
		exit(EXIT_FAILURE);
	}

	switch (fem->order) {
		case 1:
			port->l_edge = 3;
			port->l_node = 3;
			port->l_order = 1;
			break;
		case 2:
			port->l_edge = 6;
			port->l_node = 6;
			port->l_order = 2;
			break;
		case 3:
			port->l_edge = 8;
			port->l_node = 6;
			port->l_order = 3;
			break;
	}
	port->l_en = port->l_node + port->l_edge;
	

	fprintf(fp, "%lld\n", fem->unknown);
	fprintf(fp, "%lf\n", par->wavelength);
	fprintf(fp, "%d\n", 1);
	fprintf(fp, "%d\n", 11);
	fprintf(fp, "%lld\n", par->n_material);
	for (i = 0; i < par->n_material; i++) {
		fprintf(fp, "%15.10lf %15.10lf\n", sqrt(par->er[i]), 0.0);
	}
	fprintf(fp, "%d\n", 1);
	fprintf(fp, "%lf %lf\n", 0.0, 0.0);
	fprintf(fp, "%lld %lld\n", port->n_edge, port->np);//%5d
	fprintf(fp, "%lld\n", port->nr);//%5d
	count = 0;
	for (i = 0; i < fem->np; i++) {
		if (port->toLn[i] != -1) {
			xx = fem->node[i].x - port->x;
			zz = fem->node[i].z - port->z;
			//%5d
			fprintf(fp, "%lld %15.10lf %15.10lf\n", port->toLn[i] + 1, ROT1(xx,zz,cosA,sinA), fem->node[i].y);
			//fprintf(fp, "%lld %15.10lf %15.10lf\n", port->toLn[i] + 1, fem->node[i].z - port->z, fem->node[i].y);
			count++;
		}
	}
	fprintf(fp, "%lld\n", fem->order);
	fprintf(fp, "%lld\n", port->ne);
	for (i = 0; i < port->ne; i++) {
			fprintf(fp, "%lld ", port->element[i].matID + 1);
		for (j = 0; j < port->l_edge; j++) {
			fprintf(fp, " %lld", port->toLe[port->element[i].ke[j]] + 1);
		}
		for (j = 0; j < port->l_node; j++) {
			fprintf(fp, " %lld", port->toLn[port->element[i].kn[j]] + 1);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "%d %d\n", OUT_2D_SURF_SLV_DIV_X, OUT_2D_SURF_SLV_DIV_Y);

	if (fclose(fp) == EOF) {
		fprintf(stderr, "can't close file (%s).\n", string);
		exit(EXIT_FAILURE);
	}



	if (data->calcflag != 2) {
		sprintf(string,"%s _input_port%lld  -o 9 -n %lf -real -pict 320",TWO_DIM_SOLVER_PATH,n_port,par->neff);
		fprintf(stderr, ">> %s\n\n", string);
		system(string);
	}else if(data->calcflag == 2) {

	}
	return;
}

void normalizeInputPortPower(DataTable *data, PortInfo *port,
					long long int ip, 
					long long int mode){
	long long int i, j, k, it, l_it, matID;
	// PortInfo *port = &(data->port[ip]);
	double x[N_NODE], y[N_NODE], z[N_NODE], xr[N_NODE];
	std::complex<double> BB[N_EN][N_EN], BBeh[N_EN][N_EN];
	// std::complex<double> vv[N_EN];
	std::complex<double> ff[N_EN];
	// std::complex<double> wk[N_EN], wk_nEH[N_EN];
	std::complex<double> wk_n[N_EN];
	std::complex<double> amp_n;
	// std::complex<double> amp, amp_nEH;
	double rr = 1.0;
	double ZZ, k0, beta;
	// std::complex<double> cj = std::complex<double>(0, 1);;
	double cosA, sinA;
	// long long int ss[N_NODE];
	long long int l_edge = 8, l_node = 6;
	
	FEM *fem = &(data->fem);
	Param *par = &(data->par);

	//  for (i = 0; i < port->nr; i++) {
	//    fprintf(stderr, "phi[%d] = complex(%lf, %lf)\n", i, real(port->m_func[mode][i]), imag(port->m_func[mode][i]));
	//}
	
	switch (fem->order) {
		case 1:
			l_edge = 3;
			l_node = 3;
			break;
		case 2:
			l_edge = 6;
			l_node = 6;
			break;
		case 3:
			l_edge = 8;
			l_node = 6;
			break;
	}
	
	ZZ = 376.730313461;
	k0 = 2.0 * M_PI / par->wavelength;
	beta = std::real(port->beta[mode]);
	if (fem->unknown == 1) {
		rr = k0 * ZZ / beta;
	} else {
		rr = k0 / ZZ / beta;
	}
	
	cosA = cos(port->angle * M_PI / 180.0);
	sinA = sin(port->angle * M_PI / 180.0);
	
	amp_n = 0.0;
	// amp = amp_nEH = 0.0;
	for (k = 0; k < port->ne; k++) {
		for (i = 0; i < l_node; i++) {
			it = port->element[k].kn[i];
			x[i] = fem->node[it].x - port->x;
			y[i] = fem->node[it].y - port->y;
			z[i] = fem->node[it].z - port->z;
			xr[i] = (x[i] * cosA - z[i] * sinA);
		}
		
		matID = port->element[k].matID;
		elementMatrix2D(data, xr, y, matID, BB, BBeh);
		// BBehはパワーマトリックス
		// BBehは現在使用されておらず
		
		for (i = 0; i < l_edge; i++) {
			it = port->element[k].ke[i];
			l_it = port->toLe[it];
			ff[i] = port->m_func[mode][l_it];
		}
		
		for (i = 0; i < l_node; i++) {
			it = port->element[k].kn[i];
			l_it = port->toLn[it];
			ff[i + port->l_edge] = port->m_func[mode][l_it];
		}
		
		for (i = 0; i < port->l_edge; i++) {
			// wk[i] = 0.0;
			wk_n[i] = 0.0;
			for (j = 0; j < port->l_en; j++) {
				wk_n[i] += BB[i][j] * conj(ff[j]);
			}
		}
		for (i = 0; i < port->l_edge; i++) {
			amp_n += wk_n[i] * ff[i];
		}
	
	}
	
	amp_n /= rr * 2.0;
	// amp_n = sqrt(std::abs(std::real(amp_n)));
	amp_n = sqrt(std::abs(amp_n));
	
	fprintf(stderr, "Port %lld : %lld-th mode, wl = %lf\n", ip, mode, data->par.wavelength);
	fprintf(stderr, "\t\tamp_n = %lf %lf (rr = %lf)\n", std::real(amp_n), std::imag(amp_n), rr);
	fprintf(stderr, "\t\tbeta = %lf %lf\n", std::real(port->beta[mode]), std::imag(port->beta[mode]));
	
	for (i = 0; i < port->n_edge + port->np; i++) {
		port->m_func[mode][i] /= amp_n;
	}

	return;
}
