#include "fem3d.h"
#include <omp.h>
#include <sys/time.h>
#include <linux/kernel.h>
#include <sys/sysinfo.h>

void analysis(DataTable *data) {
	struct timeval s, e;
		
	memoryAllocation(data);
	fprintf(stderr, "before matrix\n");
	gettimeofday(&s, NULL);
	matrix(data); // A行列もとめる
	gettimeofday(&e, NULL);
	fprintf(stderr, "Proceeded Matrix: %.2lf(s) \n", 
		(e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6);
		
	fprintf(stderr, "input condition\n");
	inputCondition(data); // bベクトルもとめる
	
	fprintf(stderr, "input condition end\n");
	solveMatrixEquation(data);// [A]行列を2次元から1次元に，cgからia,jaに

	fprintf(stderr, "solveMatrixEquation end\n");
	return;
}

void memoryAllocation(DataTable *data) {
	long long int i;
	FEM *fem = &(data->fem);
	fem->nbw = checkBandWidth(data);
	fprintf(stderr, "nr_edge = %lld \n", fem->nr_edge);
	ALLOCATION(fem->cg_table, CGtable, fem->nr_edge);
		
	fprintf(stderr, "Make cg_table \n");
	makeCGmatrixTable(fem->element, fem->ne, fem->nr_edge, 
				fem->cg_table, fem->n_en);
		
	ALLOCATION(fem->ok, std::complex<double> *, fem->nr_edge);
	for (i = 0; i < fem->nr_edge; i++) {
		ALLOCATION(fem->ok[i], std::complex<double>, fem->cg_table[i].n);
	}
		
	fprintf(stderr, "Make cg_table complete \n");
	
	fprintf(stderr, "memory was allocated\n");
		
	return;
}

long long int checkBandWidth(DataTable *data) {
	long long int i, imax;
	long long int *a = NULL;
	FEM *fem = &(data->fem);
		
	ALLOCATION(a, long long int, fem->ne);
			
	for (i = 0; i < fem->ne; i++) {
		a[i] = maxValue(fem->element[i].kk, fem->n_en, fem->nr_edge)
			- minValue(fem->element[i].kk, fem->n_en, fem->nr_edge);
	}
	
	imax = maxValue(a, fem->ne, fem->nr_edge);
	
	free(a);
	
	return imax;
}

long long int maxValue(long long int *a, long long int n, long long int nr) {
	long long int i;
	long long int imax = 0;
		
	for (i = 0; i < n; i++) {
		if (a[i] > imax && a[i] < nr)
			imax = a[i];
	}
	
	return imax;
}

long long int minValue(long long int *a, long long int n, long long int nr) {
	long long int i;
	long long int imin = nr;
		
	for (i = 0; i < n; i++) {
		if (a[i] < imin && a[i] < nr)
			imin = a[i];
	}
		
	return imin;
}

void matrix(DataTable *data) {
	long long int i, j, k, ii, jj, itmp, matID;
	double x1[10], y1[10], z1[10];
	std::complex<double> sk[24][24], sm[24][24];
	long long int flag;
	double z_grav;
	FEM *fem = &(data->fem);
	Param *par = &(data->par);
	long long int memtmp;
	double memdbl;
	long long int num_t;
	
	if(data->nm > 0){
		omp_set_num_threads(data->nm);
		num_t = data->nm;
	}else{
		omp_set_num_threads(atoi(getenv("OMP_NUM_THREADS")));
		num_t = atoi(getenv("OMP_NUM_THREADS"));
	}
	
	for (i = 0; i < fem->nr_edge; i++) {
		for (j = 0; j < fem->cg_table[i].n; j++) {
			fem->ok[i][j] = 0.0;
		}
	}
	

	//Multi-thred setting
	if(num_t > 1){
		// if(1){

		std::complex<double> ****array;
		ALLOCATION(array, std::complex<double> ***, num_t);
		memtmp = 0;
		for(i = 0; i < fem->nr_edge; i++){
			for(long long int j = 0; j < fem->cg_table[i].n; j++){
	memtmp++;
			}
		}
				
		memtmp = sizeof(std::complex<double>) * num_t * memtmp;
		memtmp /= 1024 * 1024;
		memdbl = (double)memtmp / 1024.0;
		// 各行列の実体は作らず，ポインタだけ確保しているので，
		// 代入しない要素分はメモリ消費しない？みたい
		fprintf(stderr, "Max required memory (for matrix) = %.0lf GB ( = %.2lf * %lld)\n", memdbl, memdbl/num_t, num_t);

		for(k = 0; k < num_t; k++){
			fprintf(stderr, "allocation %lld / %lld completed ... \n", k, num_t);
			ALLOCATION(array[k], std::complex<double> **, fem->nr_edge);
			for (long long int i = 0; i < fem->nr_edge; i++) {
	ALLOCATION(array[k][i], std::complex<double> *, fem->cg_table[i].n);
	for(long long int j = 0; j < fem->cg_table[i].n; j++){
		array[k][i][j] = NULL;
	}
			}
		}
		fprintf(stderr, "allocation %lld / %lld completed ... \n", num_t, num_t);	
		
		fprintf(stderr, "calculating ... \n");
#pragma omp parallel for
		for (long long int kk = 0; kk < fem->ne; kk++){
			long long int tn = omp_get_thread_num();
			mprocess(data, array[tn], kk);
		}	

		fprintf(stderr, "summarizing ... \n");
		for(k = 0; k < num_t; k++){
			fprintf(stderr, "%lld / %lld proceeded... \n", k, num_t);
			
#pragma omp parallel for
			for(long long int ii = 0; ii < fem->nr_edge; ii++){
	for(long long int i = 0; i < fem->cg_table[ii].n; i++){
		if(array[k][ii][i]){
			fem->ok[ii][i] += *array[k][ii][i];
		}
	}
			}
		}
		fprintf(stderr, "%lld / %lld proceeded... \n", num_t, num_t);
		
		fprintf(stderr, "free array ... \n");
		for(long long int k = 0; k < num_t; k++){
			for(long long int i = 0; i < fem->nr_edge; i++) {
	for(long long int j = 0; j < fem->cg_table[i].n; j++){
		free(array[k][i][j]);
	}
	free(array[k][i]);
			}
			free(array[k]);
		}
		free(array);
		
	}

	if(data->np > 0){
		omp_set_num_threads(data->np);
	}else{
		omp_set_num_threads(atoi(getenv("OMP_NUM_THREADS")));	
	}

	return;
}

void mprocess(DataTable *data, std::complex<double> ***array,
				long long int k){

	double x1[10], y1[10], z1[10];
	std::complex<double> sk[24][24], sm[24][24];
	FEM *fem = &(data->fem);
	Param *par = &(data->par);
	//std::complex<double> tmp;
	long long int itmp;
	int ii, jj, kk = 0;
	int tmp = 0;

	for (long long int i = 0; i < fem->n_node; i++) {
		itmp = fem->element[k].kn[i];
		x1[i] = fem->node[itmp].x;
		y1[i] = fem->node[itmp].y;
		z1[i] = fem->node[itmp].z;
	}
	long long int matID = fem->element[k].matID;
	
	// if(k == 0) fprintf(stderr, "k0 = %d\n", k);

	// check mozaic pattern
	if(data->mozaic.mozflag == 1){
		if(matID == data->mozaic.mozaicmatID 
			 || matID == data->mozaic.mozaicmatID+1){
			ii = data->mozaic.ii_pixel[k];
			jj = data->mozaic.jj_pixel[k];
	
			// fprintf(stderr, "myID = %d, ii,jj = %d %d\n", data->par.myid, ii, jj);

			if(data->mozaic.BorW[kk][ii][jj] == 1){
	//if(1){
	matID = data->mozaic.coreID;
			}else{
	matID = data->mozaic.cladID;
			}

			/*
			fprintf(stderr, "tmp %d k %d, ii = %d, jj = %d, matID = %d, modified ID = %d\n", 
			 tmp, k, ii, jj, 
			 fem->element[k].matID, matID);
			*/

			tmp++;

		}
	}

	//	if(k == 0) fprintf(stderr, "k0 = %d\n", k);

	// fprintf(stderr, "k %d, ii = %d, jj = %d, matID = %d\n", k, ii, jj, matID);

	elementMatrix(data, x1, y1, z1, matID, sk, sm);

	// fprintf(stderr, "k = %d\n", k);

	for (long long int j = 0; j < fem->n_en; j++) {
		for (long long int l = 0; l < fem->n_en; l++) {
			long long int jj = fem->element[k].kk[j];
			long long int ii = fem->element[k].kk[l];
			if (ii < fem->nr_edge && jj < fem->nr_edge && ii >= 0 && jj>= 0) {
	long long int flag = 0;
	if (fem->inFlag[jj] == 3) {
		if (fem->inFlag[ii] == 1) {
			flag = 0;
		} else if (fem->inFlag[ii] == 2) {
			flag = 1;
		} else {
			double z_grav = (z1[0] + z1[1] + z1[2] + z1[3]) / 4.0;
			if (z_grav < par->INPUT_Z) {
				flag = 0;
			} else {
				flag = 1;
			}
		}
	}

	itmp = columnInCGmatrix(ii, jj, fem->cg_table, fem->nr_edge);
	//fem->ok[ii][itmp] += (sk[l][j] - par->k02 * sm[l][j]);
	if(array[ii][itmp]==NULL){
		ALLOCATION(array[ii][itmp], std::complex<double>, 1);
	}
	*array[ii][itmp] += (sk[l][j] - par->k02 * sm[l][j]);

			}
		}
	}

	// fprintf(stderr, "kk = %d\n", k);

}


// analysis()内,matrix()あとに
void inputCondition(DataTable *data) {
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

	ALLOCATION(fem->phi	 , std::complex<double>, nrhs * fem->n_edge);
	ALLOCATION(fem->phi_in, std::complex<double>, nrhs * fem->n_edge);
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

			//			fprintf(stderr, "%d %d %d\n", ii, mode, port->ne);

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

				} else {
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
			}
			
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
