#include "fem3d.h"
#include <omp.h>
#include <sys/time.h>
#include <linux/kernel.h>
#include <sys/sysinfo.h>

// function for creating FEM matrix and solve Equation
// without memory allocation

// for iterative simulation only

void analysis2(DataTable *data, int flag) {
	struct timeval s, e;
	long long int nn;
	long long int nc;
	double time1, time2;

	int i, check = 0;
	int ii, jj;

	fprintf(stderr, "WAVELENGTH = %f\n", data->par.wavelength);
	data->par.k0 = 2.0 * M_PI / data->par.wavelength;
	data->par.k02 = data->par.k0 * data->par.k0;
	fprintf(stderr, "k0 = %f\n", data->par.k0);

	if(flag == 0 && check == 0){

		gettimeofday(&s, NULL);
		analysis(data); // 全体行列作成, Ax=bのbも用意
		gettimeofday(&e, NULL);
		time1 = (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec) * 1.0e-6;
		fprintf(stderr, "---------------- after analysis\n");

		nn = data->matrix.nn;
		nc = data->matrix.nc;
		fprintf(stderr, "nn = %lld  nc = %lld \n", nn, nc);

		fprintf(stderr, "intpol(str) start \n");
		intpol(data, -1, -1, 0); //まずはテーブルつくる
		fprintf(stderr, "intpol(str) end \n");

		check = 1;

	}else{

		//  memoryAllocation(data);
		fprintf(stderr, "before matrix\n");
		gettimeofday(&s, NULL);
		matrix(data); // A行列もとめる
		gettimeofday(&e, NULL);
		fprintf(stderr, "Proceeded Matrix: %.2lf(s) \n",
			(e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6);

		fprintf(stderr, "input condition\n");
		inputCondition2(data); // bベクトルもとめる

		solveMatrixEquation2(data);// [A]行列を2次元から1次元に，cgからia,jaに

	}

	nn = data->matrix.nn;
	nc = data->matrix.nc;
	fprintf(stderr, "nn = %lld  nc = %lld \n", nn, nc);

	// Pardiso
	fprintf(stderr, "\n[PARDISO] in analysis2\n");
	gettimeofday(&s, NULL);


	callPARDISO(nn, nc, data->matrix.ia, data->matrix.ja,
	data->matrix.avals, data->fem.phi,
	data->matrix.b,
	data->arguments, data->nrhs);
	gettimeofday(&e,NULL);
	time2 = (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec) * 1.0e-6 ;
	//  fprintf(fp, "time2:\t%lf\n", time2);

	fprintf(stderr, "calcPower start \n");
	// calculate field at each port
	outputPortField(data, flag, i);
	// overlap integral
	calcPower2(data);
	return;
}


// analysis()内,matrix()あとに

void inputCondition2(DataTable *data) {
	long long int i, j, k, it, l_it, matID;
	long long int ii, jj, nrhs, mode;
	double cosA, sinA;
	double x[N_NODE], y[N_NODE], z[N_NODE], xr[N_NODE];
	std::complex<double> phi[N_EN], phi_in[N_EN], ff[N_EN];
	std::complex<double> BB[N_EN][N_EN], BBeh[N_EN][N_EN];
	std::complex<double> cj = std::complex<double>(0, 1);
	// long long int ss[N_NODE];
	double pn;

	FEM *fem = &(data->fem);
	Param *par = &(data->par);
	PortInfo *port;

	pn = 1.0 / sqrt(par->power_adjust);

	nrhs = 0;
	for(ii = 0; ii < fem->n_port; ii++){
		port = &(data->port[ii]);
		nrhs += port->modenum;
	}
	data->nrhs = nrhs;
	fprintf(stderr, "inputCondition() foa all port (nrhs = %lld)\n", data->nrhs);

	//ALLOCATION(fem->phi   , std::complex<double>, nrhs * fem->n_edge);
	//ALLOCATION(fem->phi_in, std::complex<double>, nrhs * fem->n_edge);
	for (i = 0; i < nrhs * fem->n_edge; i++){
		fem->phi[i] = 0.0;
		fem->phi_in[i] = 0.0;
	}

	jj = 0;
	for(ii = 0; ii < fem->n_port; ii++){
		port = &(data->port[ii]);
		for (mode = 0; mode < port->modenum; mode++){
			cosA = cos(port->angle * M_PI / 180.0);
			sinA = sin(port->angle * M_PI / 180.0);
			/* setField */
			for (k = 0; k < port->ne; k++) {
	for (i = 0; i < port->l_node; i++) {
		it = port->element[k].kn[i];

		x[i] = fem->node[it].x - port->x;
		y[i] = fem->node[it].y - port->y;
		z[i] = fem->node[it].z - port->z;
		xr[i] = (x[i] * cosA - z[i] * sinA);
	}
	matID = port->element[k].matID;

	elementMatrix2DP(data, xr, y, matID, BB, BBeh);

	if(data->mozaic.mozDev == 1 && data->mozaic.DBSflag == 1
		 && data->mozaic.symflag != 1){
		// for mozaic divider
		for (i = 0; i < port->l_edge; i++) {
			it = port->element[k].ke[i];
			l_it = port->toLe[it];
			if(ii == 1 && mode == 0){
				// fprintf(stderr, "%lf %lf\n", real(port->beta[mode]), real(port->beta[mode+2]));
				phi[i] = phi_in[i]
		= (1.0/sqrt(2.0))*(port->m_func[mode][l_it]*port->beta[mode]
					+port->m_func[mode+2][l_it]*port->beta[mode+2]);
			}else{
				phi[i] = phi_in[i] = port->m_func[mode][l_it]*port->beta[mode];
			}
		}

		for (i = 0; i < port->l_node; i++) {
			it = port->element[k].kn[i];
			l_it = port->toLn[it];
			if(ii == 1 && mode == 0){
				phi[i + port->l_edge] = phi_in[i + port->l_edge]
		= (1.0/sqrt(2.0))*(port->m_func[mode][l_it]*port->beta[mode]
					+port->m_func[mode+2][l_it]*port->beta[mode+2]);
			}else{
				phi[i + port->l_edge] = phi_in[i + port->l_edge]
		= port->m_func[mode][l_it]*port->beta[mode];
			}
		}

		for (i = 0; i < port->l_edge; i++) {
			ff[i] = 0.0;
			for (j = 0; j < port->l_edge + port->l_node; j++) {
				//ff[i] += BB[i][j] * phi[j] * (-cj * 2.0) * port->beta[mode];
				ff[i] += BB[i][j] * phi[j] * (-cj * 2.0);
			}
		}

	}else{
		// for other devices
		for (i = 0; i < port->l_edge; i++) {
			it = port->element[k].ke[i];
			l_it = port->toLe[it];
			phi[i] = phi_in[i] = port->m_func[mode][l_it];
		}

		for (i = 0; i < port->l_node; i++) {
			it = port->element[k].kn[i];
			l_it = port->toLn[it];
			phi[i + port->l_edge] = phi_in[i + port->l_edge]
				= port->m_func[mode][l_it];
		}

		for (i = 0; i < port->l_edge; i++) {
			ff[i] = 0.0;
			for (j = 0; j < port->l_edge + port->l_node; j++) {
				ff[i] += BB[i][j] * phi[j] * (-cj * 2.0) * port->beta[mode];
			}
		}
	}

	for (i = 0; i < port->l_edge; i++) {
		it = port->element[k].ke[i];
		l_it = port->toLe[it];
		fem->phi[jj * fem->n_edge + it] += ff[i];
		fem->phi_in[jj * fem->n_edge + it] = phi_in[i];
	}
			} // end of loop k, ne

			for (i = 0; i < fem->nr_edge; i++) {
	if (port->toLe[i] >= 0 && port->toLe[i] < port->nr) {
		fem->phi[jj * fem->n_edge + i]
			*= fem->s3d[i] * port->s2d[port->toLe[i]];
	}
			}
			jj++;
		}
	}

	return;
}

void solveMatrixEquation2(DataTable *data) {
	long long int i, j;
	// long long int nrhs = data->port[0].modenum;
	std::complex<double> *B = NULL;

	FEM *fem = &(data->fem);

	// B作らなくてもいい気がするが？
	// ALLOCATION(B, std::complex<double>, fem->nr_edge * data->nrhs);

	/*
	for (j = 0; j < data->nrhs; j++) {
		for (i = 0; i < fem->nr_edge; i++) {
			B[j*fem->nr_edge +i] = fem->phi[j*fem->nr_edge + i];
		}
	}
	*/

	fprintf(stderr, "solveMatrixEquation() : wsmp_pre()\n");
	// [A]行列を2次元から1次元に，cgからia,jaに
	wsmp_pre2(fem->ok, fem->phi, fem->nr_edge, data, data->nrhs);

	//  free(B);

	return;
}

long long int wsmp_pre2(std::complex<double> **A, std::complex<double> *bb,
			long long int nf, DataTable *data,
			long long int nrhs) {
	long long int ne = nf;
	long long int ns = 0;
	long long int nn = ne - ns;
	long long int i, j;
	long long int count;

	CGtable *cg_table = data->fem.cg_table;

	count = 0;
	// wsmp
	/*
		for (i = ns; i < ne; i++) {
		for (j = 0; j < cg_table[i].n; j++) {
		if (cg_table[i].col[j] >= i)
		count++;
		}
		}
	*/
	// pardiso

	for (i = 0, count = 0; i < nn; i++) {
		for (j = 0; j < cg_table[i].n; j++) {
			count++;
		}
	}

	data->matrix.nn = ne - ns;
	data->matrix.nc = count;

	count = 0;
	for (i = ns; i < ne; i++) {
		data->matrix.ia[i - ns] = count + 1;
		for (j = 0; j < cg_table[i].n; j++) {
			if (cg_table[i].col[j] >= i) {
	data->matrix.ja[count] = cg_table[i].col[j] + 1;
	data->matrix.avals[count] = A[i][j];
	count++;
			}
		}
		for (j = 0; j < nrhs; j++) {
			data->matrix.b[j*nf + i - ns] = bb[j*nf + i];
		}
	}
	data->matrix.ia[nn] = count + 1;

	return 0;
}

void calcPower2(DataTable *data) {
	long long int ip;
	long long int ip_in, jj;
	long long int mode1, mode2;
	FILE *fp;
	char string[256];

	FEM *fem = &(data->fem);
	std::complex<double> cj(0.0, 1.0);

	int jjj = 0;

	if (data->par.numStructures != 1) sprintf(string, "%s-%1.3f-%d.power", data->name, data->par.wavelength, data->par.strID);
	else sprintf(string, "%s-%1.3f.power", data->name, data->par.wavelength);
	if ((fp = fopen(string, "w")) == NULL) {
		fprintf(stderr, "can't open file (%s).\n", string);
		exit(EXIT_FAILURE);
	}
	fprintf(stderr, "------------- Now, in calcPower2 ------------\n");

	fprintf(fp, "wl\tInputMode\tOutputPortID\tOutputMode\tTransmission\tAmplitude\tPhase[rad]\tRe\tIm\n");

	jj = 0;
	for (ip_in = 0; ip_in < fem->n_port; ip_in++) { //入射ポート

		for (mode1 = 0; mode1 < data->port[ip_in].modenum; mode1++) {//入射モード

			for (ip = 0; ip < fem->n_port; ip++) { //出射ポート

				for (mode2 = 0; mode2 < data->port[ip].modenum; mode2++) {
					//出射モード（理想フィールド）
					fprintf(stderr, "inputmode[%lld](%lld-th mode from port %d):",
						jj, mode1, mode2);
					fprintf(stderr, "beta_in = %lf\n",
						abs(data->port[ip_in].beta[mode1]));

					calcPortPower2(data, ip, ip_in, mode1, mode2, jj);

					fprintf(stderr, "  port[%lld] mode[%lld](beta = %6.3lf) = %15.10lf (phase = %lf)\n",
						ip, mode2, std::real(data->port[ip].beta[mode2]),
						data->port[ip].power, data->port[ip].phase);
					fprintf(fp, "%lf\t%lld\t%lld\t%lld\t", data->par.wavelength,
						jj, ip, mode2);

					fprintf(fp, "%.15lf\t%.15lf\t%.15lf\t",
						data->port[ip].power, sqrt(data->port[ip].power),
						data->port[ip].phase);
					fprintf(fp, "%.15lf\t%.15lf\t\n",
						sqrt(data->port[ip].power) * cos(data->port[ip].phase),
						sqrt(data->port[ip].power) * sin(data->port[ip].phase));

					if (data->mozaic.DBSflag != 0) {
						if(data->mozaic.mozDev == 0){
							// for SSC
							data->mozaic.TT[jj][jjj] = data->port[ip].power;
							fprintf(stderr, "mozaic.TT[%d][%d] = %lf\n",
								jj, jjj, data->mozaic.TT[jj][jjj]);
							jjj++;

						} else if(data->mozaic.mozDev == 1){
							// for mode divider
							if(jj == 0){
								data->mozaic.TT[jj][jjj] = data->port[ip].power;
								fprintf(stderr, "mozaic.TT[%d][%d] = %lf\n",
									jj, jjj, data->mozaic.TT[jj][jjj]);
								jjj++;
							}
						} else if(data->mozaic.mozDev == 2){
							// for mode MUX
							data->mozaic.TT[jj][jjj] = data->port[ip].power;
							fprintf(stderr, "mozaic.TT[%d][%d] = %lf\n",
								jj, jjj, data->mozaic.TT[jj][jjj]);
							jjj++;
						} else if(data->mozaic.mozDev == 3 || data->mozaic.mozDev == 4){
							// for 2 mode EX
							data->mozaic.TT[jj][jjj] = data->port[ip].power;
							fprintf(stderr, "mozaic.TT[%d][%d] = %lf\n",
								jj, jjj, data->mozaic.TT[jj][jjj]);
							jjj++;
						} else if(data->mozaic.mozDev == 5){
							// for 1 mode lense
							data->mozaic.TT[jj][jjj] = data->port[ip].power;
							fprintf(stderr, "mozaic.TT[%d][%d] = %lf\n",
								jj, jjj, data->mozaic.TT[jj][jjj]);
							jjj++;
						} else if (data->mozaic.mozDev == 6) {
							// for Waveguide Crossing
							data->mozaic.TT[jj][jjj] = data->port[ip].power;
							fprintf(stderr, "mozaic.TT[%d][%d] = %lf\n",
								jj, jjj, data->mozaic.TT[jj][jjj]);
							jjj++;
						} else if (data->mozaic.mozDev == 7) {
							// nxn Power Splitter
							if (data->mozaic.numModes == 2) {
								if (jj == 0) {
									data->mozaic.TT[jj][jjj] = data->port[ip].power;

									data->mozaic.c_TT[jj][jjj] 
										= std::complex<double>(sqrt(data->port[ip].power) * cos(data->port[ip].phase),
														  sqrt(data->port[ip].power) * sin(data->port[ip].phase));
									fprintf(stderr, "mozaic.TT[%d][%d] = %lf\n",
										jj, jjj, data->mozaic.TT[jj][jjj]);
									jjj++;
								}
							}
						}
					}
				} // end of mode2
			} // end of ip
			jj++;
			jjj = 0;
		} // end of mode1
	} // end of ip_in

	if (fclose(fp) == EOF) {
		fprintf(stderr, "can't close file ( %s )\n", string);
		exit(EXIT_FAILURE);
	}

	return;
}


void calcPortPower2(DataTable *data, long long int ip,
				long long int ip_in, long long int mode1,
				long long int mode2, long long int jj){
	long long int i, j, k, it, l_it, matID;
	PortInfo *port = &(data->port[ip]);
	long long int  nn2D = (port->n_edge + port->np);
	double x[N_NODE], y[N_NODE], z[N_NODE], xr[N_NODE];
	std::complex<double> BB[N_EN][N_EN], BBeh[N_EN][N_EN];
	std::complex<double> vv[N_EN], ff[N_EN], ff_in[N_EN];
	std::complex<double> wk[N_EN], wk_n[N_EN];
	std::complex<double> amp, amp_n;
	double ZZ, k0, beta;
	double cosA, sinA;
	// long long int ss[N_NODE];
	long long int l_edge = 8, l_node = 6;

	FEM *fem = &(data->fem);
	Param *par = &(data->par);

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

	// fprintf(stderr, "Port %lld: beta[%lld] = %lf\n", ip, mode, std::real(port->beta[mode]));

	cosA = cos(port->angle * M_PI / 180.0);
	sinA = sin(port->angle * M_PI / 180.0);

	amp = amp_n = 0.0;
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

		for (i = 0; i < l_edge; i++) {
			it = port->element[k].ke[i];
			l_it = port->toLe[it];
			if(ip_in == ip) ff_in[i] = port->m_func[mode1][l_it];

			ff[i]    = port->m_func[mode2][l_it];
			vv[i]    = port->vv[nn2D*jj + l_it];
		}

		for (i = 0; i < l_node; i++) {
			it = port->element[k].kn[i];
			l_it = port->toLn[it];
			if(ip_in == ip) ff_in[i + port->l_edge] = port->m_func[mode1][l_it];

			ff[i + port->l_edge]    = port->m_func[mode2][l_it]; // 理想フィールド
			// vv[i + port->l_edge]
			// = port->vv[nn2D*jj + l_it] / data->port[ip_in].beta[mode1];
			// 伝搬フィールド(betaあってる？)
			vv[i + port->l_edge]    = port->vv[nn2D*jj + l_it]; // 伝搬フィールド
		}

		//入射界の分をさっぴく
		if (ip == ip_in) {
			for (i = 0; i < l_edge + l_node; i++) {
				vv[i] -= ff_in[i];
			}
		}

		for (i = 0; i < port->l_edge; i++) {
			wk[i] = wk_n[i] = 0.0;
			for (j = 0; j < port->l_en; j++) {
				wk_n[i] += BB[i][j] * conj(ff[j]);
				wk[i] += BB[i][j] * conj(ff[j]);
			}
		}
		for (i = 0; i < port->l_edge; i++) {
			amp_n += wk_n[i] * (ff[i]);
			amp += wk[i] * (vv[i]);
		}
	}

	fprintf(stderr, "  amp, amp_n = %lf %lf\n", std::abs(amp), std::real(amp_n));
	amp = amp / amp_n;
	port->power = std::abs(amp) * std::abs(amp);
	// port->power = std::abs(amp) * std::abs(amp) * std::abs( data->port[ip].beta[0] / data->port[0].beta[0]);
	// port->power = std::abs(amp) * std::abs(amp) * std::abs(data->port[ip].beta[mode2] / data->port[ip_in].beta[mode1]);// (betaあってる？)
	port->phase = atan2(std::imag(amp), std::real(amp));

	return;
}
