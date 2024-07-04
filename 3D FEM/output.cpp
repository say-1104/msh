#include "fem3d.h"

void outputPortField(DataTable *data, int flag, int ii) {
	long long int i, ip;
	long long int mode, nn, nn2D;
	// long long int j, it;
	FILE *fp;
	char string[256];
	std::complex<double> phix, phiy, phiz, phit;
	std::complex<double> phix2, phiy2, phiz2;
	PortInfo *port;
	long long int count;
	// std::complex<double> *ft;
	// std::complex<double> *fz;
	double cosA, sinA;
	// double xx[3], yy[3], zz[3];
	// double *ss;
	
	FEM *fem = &(data->fem);
	nn = fem->n_edge;
	
	// fprintf(stderr, "nn(fem->n_edge) = %lld\n", nn);
	// if ((fp = fopen(".DEBUG_field0", "w")) != NULL) {
	//     for (i = 0; i < nn; i++) {
	//         fprintf(fp, "%.2e\t%.2e\n", std::real(fem->phi[i]), std::imag(fem->phi[i]));
	//     }
	//     fclose(fp);
	// }
	// if ((fp = fopen(".DEBUG_field1", "w")) != NULL) {
	//     for (i = nn; i < nn*2; i++) {
	//         fprintf(fp, "%.2e\t%.2e\n", std::real(fem->phi[i]), std::imag(fem->phi[i]));
	//     }
	//     fclose(fp);
	// }
	// if ((fp = fopen(".DEBUG_field2", "w")) != NULL) {
	//     for (i = nn*2; i < nn*3; i++) {
	//         fprintf(fp, "%.2e\t%.2e\n", std::real(fem->phi[i]), std::imag(fem->phi[i]));
	//     }
	//     fclose(fp);
	// }
	
	long long int jj, ip_in;
	// nrhs = data->port[0].modenum;
	for (ip = 0; ip < fem->n_port; ip++) { // 出射ポート
		port = &(data->port[ip]);
		cosA = cos(port->angle * M_PI / 180.0);
		sinA = sin(port->angle * M_PI / 180.0);
		nn2D = (port->n_edge + port->np);

		if(flag == 0){
		ALLOCATION(port->vv, std::complex<double>, nn2D*data->nrhs);
		}

		for (i = 0; i < nn2D*data->nrhs; i++) {
		port->vv[i] = 0.0;
		}
		
		jj = -1;
		// for (mode = 0; mode < nrhs; mode++) { // 入射モード
		for (ip_in = 0; ip_in < fem->n_port; ip_in++) { // 入射ポート
		for (mode = 0; mode < data->port[ip_in].modenum; mode++) { // 入射モード
	jj++;
	sprintf(string, "_output_wl%1.3f_port%lld-v%lld.slv", 
		data->par.wavelength, ip, jj);
	if ((fp = fopen(string, "w")) == NULL) {
		break;
	} else {
		
		fprintf(fp, "%.10e %.10e\n", 
		std::real(data->port[ip_in].beta[mode]), 
		std::imag(data->port[ip_in].beta[mode]));
		    
		// no boundary condition
		count = 0;
		for (i = 0; i < fem->np; i++) {
		if (port->toLn[i] >= 0 && port->toLn[i] < port->nr) {
		intpolAtPoint(data, fem->node[i].x, fem->node[i].y, 
		fem->node[i].z, &phix, &phiy, &phiz, 
		&phix2, &phiy2, &phiz2, jj, nn);
		// port->vv[nn2D*jj + port->toLn[i]] 
		// = (phiz * cosA + phix * sinA) / data->port[ip_in].beta[mode]; 
		port->vv[nn2D*jj + port->toLn[i]] = (phiz * cosA + phix * sinA);
		count++;
		// if(isinf(std::real(port->vv[nn2D*jj + port->toLn[i]]))){
		//fprintf(stderr, "xx");}
		}
		}
		for (i = 0; i < fem->n_edge; i++) {
		if (port->toLe[i] >= 0 && port->toLe[i] < port->nr) {
		phit = fem->phi[nn*jj + i];
		if (fem->order == 1) {
		port->vv[nn2D*jj + port->toLe[i]] 
		= phit * std::complex<double>(fem->s3d[i], 0) 
		* std::complex<double>(port->s2d[port->toLe[i]], 0);
		} else {
		port->vv[nn2D*jj + port->toLe[i]] = phit;
		}
		// if(isinf(std::real(port->vv[nn2D*jj + port->toLe[i]]))){
		// fprintf(stderr, "tt");}
		count++;
		}
		}
		// diriclet condition
		for (i = 0; i < fem->np; i++) {
		if (port->toLn[i] >= port->nr) {
		intpolAtPoint(data, fem->node[i].x, fem->node[i].y, 
		fem->node[i].z, &phix, &phiy, &phiz, 
		&phix2, &phiy2, &phiz2, jj, nn);
		port->vv[nn2D*jj + count++] = 0.0;
		}
		}
		for (i = 0; i < fem->n_edge; i++) {
		if (port->toLe[i] >= port->nr) {
		phit = fem->phi[nn*jj + i];
		port->vv[nn2D*jj + count++] = 0.0;
		}
		}
		for (i = 0; i < port->n_edge + port->np; i++) {
		// fprintf(fp, "port->vv[%lld] = %+.10e %+.10e\n", 
		// i, std::real( port->vv[i]), std::imag(port->vv[i]));
		fprintf(fp, "%+.10e %+.10e\n", 
		std::real(port->vv[nn2D*jj + i]), 
		std::imag(port->vv[nn2D*jj + i]));
		// if(isinf(std::real(port->vv[nn2D*jj + i]))){
		// fprintf(stderr, "!!");}
		}
		fclose(fp);
	}
		}
		}
	}
	
		
	return;
}

#define NX 200
#define NZ 200

void intpolAtPoint(DataTable *data, double x0, double y0, double z0,
		 std::complex<double> *Phix, std::complex<double> *Phiy, 
		 std::complex<double> *Phiz, 
		 std::complex<double> *Phix2, std::complex<double> *Phiy2, 
		 std::complex<double> *Phiz2, 
		 long long int mode, long long int nn){

	static long long int call_count = 0;
	// long long int ix, iy;
	long long int i, j, k, l, m, n, ii, mm;
	// long long int jj;
	long long int itmp;
	// long long int flag;
	std::complex<double> phix, phiy, phiz;
	std::complex<double> phix2, phiy2, phiz2;
	double xmin, xmax, ymin, ymax, zmin, zmax;
	// double xdiv, ydiv, zdiv;
	double a[4], b[4], c[4], d[4], alpha[4];
	double U[24], V[24], W[24];
	double Vx[24], Wx[24];
	double Wy[24], Uy[24];
	double Uz[24], Vz[24];
	double ss[6], aa[6], bb[6], cc[6], ll[6];
	double area[4];
	double co[24];
	static long long int edge[6][2] = { { 0, 1 }, { 0, 2 }, { 0, 3 }, 
		{ 1, 2 }, { 3, 1 }, { 2, 3 } };
	double Ve;
	double L1, L2, L3, L4;
	double L1x, L2x, L3x, L4x;
	double L1y, L2y, L3y, L4y;
	double L1z, L2z, L3z, L4z;
	double x[10], y[10], z[10];
	std::complex<double> vv[24];
	// double xg[10], yg[10], zg[10];
	// double hx, hy, hz, ht;
	double ZZ = 376.730313461;
	double k0 = 2.0 * M_PI / data->par.wavelength;
	std::complex<double> cj = std::complex<double>(0, 1);
	
	static long long int *table[NX][NZ], table_n[NX][NZ];
	double e_xmin, e_xmax, e_zmin, e_zmax;
	long long int region_x, region_z;
	long long int start_x, start_z, end_x, end_z;
	long long int iii;
	long long int matID;
	// long long int kk;
	
	FEM *fem = &(data->fem);
	
	xmin = fem->min.x;
	xmax = fem->max.x;
	ymin = fem->min.y;
	ymax = fem->max.y;
	zmin = fem->min.z;
	zmax = fem->max.z;
	
	if (call_count == 0) {
		call_count++;
		
		for (i = 0; i < NX; i++) {
		for (j = 0; j < NZ; j++) {
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
	
	if (e_xmin > x[i])
		e_xmin = x[i];
	if (e_xmax < x[i])
		e_xmax = x[i];
	if (e_zmin > z[i])
		e_zmin = z[i];
	if (e_zmax < z[i])
		e_zmax = z[i];
		}
		
		start_x = (long long int) ((e_xmin - xmin) / (xmax - xmin) * NX);
		start_z = (long long int) ((e_zmin - zmin) / (zmax - zmin) * NZ);
		end_x = (long long int) ((e_xmax - xmin) / (xmax - xmin) * NX) + 1;
		end_z = (long long int) ((e_zmax - zmin) / (zmax - zmin) * NZ) + 1;
		
		if (end_x > NX)
	end_x = NX;
		if (end_z > NZ)
	end_z = NZ;
		
		for (i = start_x; i < end_x; i++) {
	for (j = start_z; j < end_z; j++) {
		table_n[i][j]++;
	}
		}
		}
		
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
	
	if (e_xmin > x[i])
		e_xmin = x[i];
	if (e_xmax < x[i])
		e_xmax = x[i];
	if (e_zmin > z[i])
		e_zmin = z[i];
	if (e_zmax < z[i])
		e_zmax = z[i];
		}
		
		start_x = (long long int) ((e_xmin - xmin) / (xmax - xmin) * NX);
		start_z = (long long int) ((e_zmin - zmin) / (zmax - zmin) * NZ);
		end_x = (long long int) ((e_xmax - xmin) / (xmax - xmin) * NX) + 1;
		end_z = (long long int) ((e_zmax - zmin) / (zmax - zmin) * NZ) + 1;
		
		if (end_x > NX)
	end_x = NX;
		if (end_z > NZ)
	end_z = NZ;
		
		for (i = start_x; i < end_x; i++) {
	for (j = start_z; j < end_z; j++) {
		table[i][j][table_n[i][j]] = ii;
		table_n[i][j]++;
	}
		}
		}
	}
	
	region_x = (long long int) ((x0 - xmin) / (xmax - xmin) * NX);
	region_z = (long long int) ((z0 - zmin) / (zmax - zmin) * NZ);
	
	if (region_x >= NX)
		region_x = NX - 1;
	if (region_z >= NZ)
		region_z = NZ - 1;
	
	alpha[0] = alpha[2] = 1.0;
	alpha[1] = alpha[3] = -1.0;
	
	// flag = 0;
	mm = 0;
	phix = phiy = phiz = 0.0;
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
		a[k] = alpha[k] * (x[l] * (y[m] * z[n] - y[n] * z[m]) 
			 + x[m] * (y[n] * z[l] - y[l] * z[n]) 
			 + x[n] * (y[l] * z[m] - y[m] * z[l]));
		b[k] = alpha[k] * (y[l] * (z[n] - z[m]) + y[m] * (z[l] - z[n])
			 + y[n] * (z[m] - z[l]));
		c[k] = alpha[k] * (z[l] * (x[n] - x[m]) + z[m] * (x[l] - x[n])
			 + z[n] * (x[m] - x[l]));
		d[k] = alpha[k] * (x[l] * (y[n] - y[m]) + x[m] * (y[l] - y[n])
			 + x[n] * (y[m] - y[l]));
		}
		Ve = (a[0] + a[1] + a[2] + a[3]) / 6.0;
		for (i = 0; i < 6; i++) {
		aa[i] = x[edge[i][1]] - x[edge[i][0]];
		bb[i] = y[edge[i][1]] - y[edge[i][0]];
		cc[i] = z[edge[i][1]] - z[edge[i][0]];
		ll[i] = sqrt(aa[i] * aa[i] + bb[i] * bb[i] + cc[i] * cc[i]);
		/*
	xg[i] = (x[edge[i][1]]+x[edge[i][0]])/2.0;
	yg[i] = (y[edge[i][1]]+y[edge[i][0]])/2.0;
	zg[i] = (z[edge[i][1]]+z[edge[i][0]])/2.0;
		*/
		}
		for (i = 0; i < 6; i++) {
		if (cc[i] > 0.0 || (cc[i] == 0.0 && bb[i] > 0.0) 
		|| (cc[i] == 0.0 && bb[i] == 0.0 && aa[i] > 0.0)) {
	ss[i] = 1.0;
		} else {
	ss[i] = -1.0;
		}
		}
		  
		//elementmatrix(not CURVE) より移植
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
		
		L1 = (a[0] + b[0] * x0 + c[0] * y0 + d[0] * z0) / (6.0 * Ve);
		L2 = (a[1] + b[1] * x0 + c[1] * y0 + d[1] * z0) / (6.0 * Ve);
		L3 = (a[2] + b[2] * x0 + c[2] * y0 + d[2] * z0) / (6.0 * Ve);
		L4 = (a[3] + b[3] * x0 + c[3] * y0 + d[3] * z0) / (6.0 * Ve);

		if ((L1 >= -1.0e-4) && (L2 >= -1.0e-4) 
	&& (L3 >= -1.0e-4) && (L4 >= -1.0e-4)) {
		matID = fem->element[ii].matID;
		mm = 1;
		
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
	
	for (i = 0; i < 6; ++i) {
		Uy[i] *= ss[i] * ll[i];
		Uz[i] *= ss[i] * ll[i];
		Vz[i] *= ss[i] * ll[i];
		Vx[i] *= ss[i] * ll[i];
		Wy[i] *= ss[i] * ll[i];
		Wx[i] *= ss[i] * ll[i];
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
	
	//elementmatrix(not CURVE) より移植
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
	
	for (i = 0; i < 12; i++) {
		U[i] *= ll[i % 6];
		V[i] *= ll[i % 6];
		W[i] *= ll[i % 6];
	}
	for (i = 0; i < 12; i++) {
		Uy[i] *= ll[i % 6];
		Uz[i] *= ll[i % 6];
		Vz[i] *= ll[i % 6];
		Vx[i] *= ll[i % 6];
		Wy[i] *= ll[i % 6];
		Wx[i] *= ll[i % 6];
	}
	for (i = 12; i < data->fem.n_en; ++i) {
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
		
		for (i = 0; i < fem->n_en; i++) {
	itmp = fem->element[ii].kk[i];
	if (itmp >= 0) {
		vv[i] = fem->phi[nn*mode + itmp];
	} else {
		vv[i] = 0.0;
	}
		}
		
		phix = phiy = phiz = 0.0;
		phix2 = phiy2 = phiz2 = 0.0;
		for (i = 0; i < fem->n_en; i++) {
	phix += U[i] * vv[i];
	phiy += V[i] * vv[i];
	phiz += W[i] * vv[i];
		}
		
		// added by satonori 20180614
		// magnetic fields if phi is the electric field
		for (i = 0; i < fem->n_en; i++) {
	phix2 += (Wy[i] - Vz[i]) * vv[i];
	phiy2 += (Uz[i] - Wx[i]) * vv[i];
	phiz2 += (Vx[i] - Uy[i]) * vv[i];
		}
		if (fem->unknown == 1) {
	phix2 /= -cj * k0 * ZZ;
	phiy2 /= -cj * k0 * ZZ;
	phiz2 /= -cj * k0 * ZZ;
		} else {
	phix2 /= cj * k0 / ZZ * data->par.ermatrix[matID][0][0];
	phiy2 /= cj * k0 / ZZ * data->par.ermatrix[matID][1][1];
	phiz2 /= cj * k0 / ZZ * data->par.ermatrix[matID][2][2];
		}
		break;
		}
	}
	if (mm != 1) {
		fprintf(stderr, "------------ Not Found (%lf, %lf, %lf)\n", x0, y0, z0);
	}
	
	*Phix = phix;
	*Phiy = phiy;
	*Phiz = phiz;
	*Phix2 = phix2;
	*Phiy2 = phiy2;
	*Phiz2 = phiz2;
	
	return;
}

void inputPortField(DataTable *data) {
	long long int i, ip;
	FILE *fp;
	char string[256];
	// std::complex<double> phiz, phit;
	PortInfo *port;
	// long long int count;
	// std::complex<double> *ft;
	// std::complex<double> *fz;
	std::complex<double> cj = std::complex<double>(0, 1);
	double re, im;
	
	FEM *fem = &(data->fem);
	
	for (ip = 0; ip < fem->n_port; ip++) {
		port = &(data->port[ip]);
		ALLOCATION(port->vv, std::complex<double>, port->n_edge+port->np);
		sprintf(string, "_output_port%lld.Vslv", ip);
		if ((fp = fopen(string, "r")) == NULL) {
		break;
		} else {
		for (i = 0; i < port->n_edge + port->np; i++) {
	fscanf(fp, "%lf %lf", &re, &im);
	port->vv[i] = re + cj * im;
		}
		fscanf(fp, "%lf", &re);
		port->beta[0] = re;
		
		fclose(fp);
		}
	}
	
	return;
}
