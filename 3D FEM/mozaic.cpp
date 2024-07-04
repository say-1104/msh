#include "fem3d.h"
#include <getopt.h>

// reading mozaic data, memory allocation

extern int inputMozaic(DataTable *data) {

	int i, j, k, jj;
	int tmp, nrhs;
	char filename[256];
	FILE *fp;

	sprintf(filename, "mozaic.pre");
	if ((fp = fopen(filename, "r")) == NULL) {
	fprintf(stderr, "can't open file ( %s )\n", filename);
	exit(1);
	}

	fscanf(fp, "%*s %lf %lf %lf", 
	 &data->mozaic.xs, &data->mozaic.ys, &data->mozaic.zs);
	fprintf(stderr,"Starting x, y, z (left upper position) = %lf, %lf, %lf\n", 
		data->mozaic.xs, data->mozaic.ys, data->mozaic.zs);

	fscanf(fp, "%*s %d %d %d", 
	 &data->mozaic.Nx, &data->mozaic.Ny, &data->mozaic.Nz);
	fprintf(stderr,"The number of mozaic region in x,z = %d, %d\n", 
		data->mozaic.Nx, data->mozaic.Nz);

	fscanf(fp, "%*s %lf %lf %lf", 
	 &data->mozaic.dx, &data->mozaic.dy, &data->mozaic.dz);
	fprintf(stderr,"dx, dy, dz = %lf, %lf, %lf\n", 
		data->mozaic.dx, data->mozaic.dy, data->mozaic.dz);

	fscanf(fp, "%*s %d", &data->mozaic.symflag);
	fscanf(fp, "%*s %s", &data->mozaic.filename);
	fscanf(fp, "%*s %d", &data->mozaic.mozaicmatID);
	fscanf(fp, "%*s %d", &data->mozaic.coreID);
	fscanf(fp, "%*s %d", &data->mozaic.cladID);
	fscanf(fp, "%*s %lf", &data->mozaic.gamma);

	fscanf(fp, "%*s %d", &data->mozaic.DBSitr);
	fscanf(fp, "%*s %lf", &data->mozaic.termination);
	fscanf(fp, "%*s %d", &data->mozaic.PixelReduction);

	fprintf(stderr,"symflag = %d\n", data->mozaic.symflag);

	fprintf(stderr,"Filename of mozaic = %s\n",
		data->mozaic.filename);

	fprintf(stderr,"mozaicID, coreID, cladID = %d, %d, %d\n",
		data->mozaic.mozaicmatID, data->mozaic.coreID,
		data->mozaic.cladID);

	fprintf(stderr,"DBSitr = %d\n", data->mozaic.DBSitr);
	fprintf(stderr,"Termination = %lf\n", data->mozaic.termination);
	fprintf(stderr,"PixelReduction = %d\n", data->mozaic.PixelReduction);

	data->mozaic.totalN = 1;
	ALLOCATION(data->mozaic.BorW, int **, data->mozaic.totalN);
	for(i = 0; i < data->mozaic.totalN; i++) {
		ALLOCATION(data->mozaic.BorW[i], int *, data->mozaic.Nx);
	}
	for(i = 0; i < data->mozaic.totalN; i++) {
		for(j = 0; j < data->mozaic.Nx; j++) {
			ALLOCATION(data->mozaic.BorW[i][j], int, data->mozaic.Nz);
		}
	}

	ALLOCATION(data->mozaic.x0, double *, data->mozaic.Nx);
	ALLOCATION(data->mozaic.x1, double *, data->mozaic.Nx);
	ALLOCATION(data->mozaic.y0, double *, data->mozaic.Nx);
	ALLOCATION(data->mozaic.y1, double *, data->mozaic.Nx);
	ALLOCATION(data->mozaic.z0, double *, data->mozaic.Nx);
	ALLOCATION(data->mozaic.z1, double *, data->mozaic.Nx);

	for(i = 0; i < data->mozaic.Nx; i++) {
		ALLOCATION(data->mozaic.x0[i], double, data->mozaic.Nz);
		ALLOCATION(data->mozaic.x1[i], double, data->mozaic.Nz);
		ALLOCATION(data->mozaic.y0[i], double, data->mozaic.Nz);
		ALLOCATION(data->mozaic.y1[i], double, data->mozaic.Nz);
		ALLOCATION(data->mozaic.z0[i], double, data->mozaic.Nz);
		ALLOCATION(data->mozaic.z1[i], double, data->mozaic.Nz);
	}

	fclose(fp);

	// -------------------------------------

	if(data->mozaic.DBSflag != 2){
		// Guessed matrix pattern
		
		//    sprintf(filename, "GuessedMatrix.csv");
		if ((fp = fopen(data->mozaic.filename, "r")) == NULL) {
			fprintf(stderr, "can't open file ( %s )\n", data->mozaic.filename);
			exit(1);
		}
		
		k = 0;
		for (i = 0; i < data->mozaic.Nx; i++) {
			for (j = 0; j < data->mozaic.Nz; j++) {
				if (i == data->mozaic.Nx-1 && j == data->mozaic.Nz-1){
					fscanf(fp, "%d", &data->mozaic.BorW[k][i][j]);
				} else {
					fscanf(fp, "%d,", &data->mozaic.BorW[k][i][j]);
				}
			}
		}

		for (i = 0; i < data->mozaic.Nx; i++) {
			for (j = 0; j < data->mozaic.Nz; j++) {
				fprintf(stderr, "%d ", data->mozaic.BorW[k][i][j]);
			}
			fprintf(stderr, "\n");
		}
		fprintf(stderr, "\n");

		data->mozaic.starti = 0;
		data->mozaic.startj = 0;
		data->mozaic.starttotal = 0;

	} else if(data->mozaic.DBSflag == 2){
		sprintf(filename, "GuessedMatrix-middle.csv");
		if ((fp = fopen(filename, "r")) == NULL) {
			fprintf(stderr, "can't open file ( %s )\n", filename);
			exit(1);
		}
		
		k = 0;
		for(i = 0; i < data->mozaic.Nx; i++) {
			for(j = 0; j < data->mozaic.Nz; j++) {
				if (i == data->mozaic.Nx-1 && j == data->mozaic.Nz-1){
					fscanf(fp, "%d", &data->mozaic.BorW[k][i][j]);
				} else {
					fscanf(fp, "%d,", &data->mozaic.BorW[k][i][j]);
				}
			}
		}

		for(i = 0; i < data->mozaic.Nx; i++) {
			for(j = 0; j < data->mozaic.Nz; j++) {
				fprintf(stderr, "%d ", data->mozaic.BorW[k][i][j]);
			}
			fprintf(stderr, "\n");
		}
		fprintf(stderr, "\n");

		fscanf(fp, "%*s %d %*s %d %*s %d", 
			&data->mozaic.starti, &data->mozaic.startj, 
			&data->mozaic.starttotal);
		fprintf(stderr, "GuessedMatrix-middle.csv is read\n"); 
	}

	fprintf(stderr, "starti = %d, startj = %d, starttotal = %d\n", 
		data->mozaic.starti, data->mozaic.startj, data->mozaic.starttotal);


	// checking minimum and maximum x,z coordinate in each pixel
	for(i = 0; i < data->mozaic.Nx; i++) {
		for(j = 0; j < data->mozaic.Nz; j++) {
			if (data->mozaic.BorW[k][i][j] != -1) {
				data->mozaic.z0[i][j] = (double)j*data->mozaic.dz+data->mozaic.zs;
				data->mozaic.z1[i][j] = (double)(j+1)*data->mozaic.dz+data->mozaic.zs;
				data->mozaic.x0[i][j] = data->mozaic.xs-(double)(i+1)*data->mozaic.dx;
				data->mozaic.x1[i][j] = data->mozaic.xs-(double)(i)*data->mozaic.dx;

				data->mozaic.y0[i][j] = data->mozaic.ys-data->mozaic.dy;
				data->mozaic.y1[i][j] = data->mozaic.ys;

				/*
				fprintf(stderr, "%d %d pixel: x(%lf %lf), y(%lf %lf), z(%lf %lf)\n", 
					i, j, 
					data->mozaic.x0[i][j], data->mozaic.x1[i][j],
					data->mozaic.y0[i][j], data->mozaic.y1[i][j],
					data->mozaic.z0[i][j], data->mozaic.z1[i][j]);
				*/
			}
		}
	}
	//  exit(0);

	// parameters for pixel
	ALLOCATION(data->mozaic.matID, int, data->fem.ne);
	ALLOCATION(data->mozaic.ii_pixel, int, data->fem.ne);
	ALLOCATION(data->mozaic.jj_pixel, int, data->fem.ne);

	for(k = 0; k < data->fem.ne; k++) {
		data->mozaic.matID[k] = -1;
		data->mozaic.ii_pixel[k] = -1;
		data->mozaic.jj_pixel[k] = -1;
	}
	
	ALLOCATION(data->mozaic.ne_pixel, int *, data->mozaic.Nx);
	ALLOCATION(data->mozaic.kk_pixel, int **, data->mozaic.Nx);

	for(i = 0; i < data->mozaic.Nx; i++) {
		ALLOCATION(data->mozaic.ne_pixel[i], int, data->mozaic.Nz);
		ALLOCATION(data->mozaic.kk_pixel[i], int *, data->mozaic.Nz);
	}

	for(i = 0; i < data->mozaic.Nx; i++) {
		for(j = 0; j < data->mozaic.Nz; j++) {
			data->mozaic.ne_pixel[i][j] = 0;
			ALLOCATION(data->mozaic.kk_pixel[i][j], int, 500);
		}
	}

	tmp = 0;
	for (i = 0; i < data->fem.n_port; i++) {
		tmp += data->port[i].modenum;
	}

	data->mozaic.numofT = tmp;
	fprintf(stderr, "Total transmission number = %d\n", data->mozaic.numofT);

	ALLOCATION(data->mozaic.TT, double*, data->mozaic.numModes);
	ALLOCATION(data->mozaic.c_TT, std::complex<double>*, data->mozaic.numModes);

	for (i = 0; i < data->mozaic.numofT; i++) {
		ALLOCATION(data->mozaic.TT[i], double, data->mozaic.numofT);
		ALLOCATION(data->mozaic.c_TT[i], std::complex<double>, data->mozaic.numofT);
	}
	// ----------------------
}

// Check mesh data
extern int checkMeshforMozaic(DataTable *data) {

	int i, j, k, l, m, n;
	int ii, jj, kk;
	int tmp = 0;

	int itmp;
	double x[10], y[10], z[10];
	double xg, yg, zg;
	double y0, y1;
	int matID;
	int flag1, flag2, flag3, flag4;

	FILE *fp;
	char string[256];

	FEM *fem = &(data->fem);
	
	for (k = 0; k < fem->ne; k++) {

		for (long long int i = 0; i < fem->n_node; i++) {
			itmp = fem->element[k].kn[i];
			x[i] = fem->node[itmp].x;
			y[i] = fem->node[itmp].y;
			z[i] = fem->node[itmp].z;
		}

		matID = fem->element[k].matID;

		if(matID == data->mozaic.mozaicmatID
			|| matID == data->mozaic.mozaicmatID+1){
			// search the pixel that the element belong
			xg = (x[0]+x[1]+x[2]+x[3])/4.0;
			yg = (y[0]+y[1]+y[2]+y[3])/4.0;
			zg = (z[0]+z[1]+z[2]+z[3])/4.0;

			// fprintf(stderr, "xg,yg,zg = %lf %lf %lf\n", xg, yg, zg);

			tmp++;

			for(ii = 0; ii < data->mozaic.Nx; ii++) {
				for(jj = 0; jj < data->mozaic.Nz; jj++) {

					flag1 = flag2 = flag3 = flag4 = 0;

					if((x[0] > data->mozaic.x0[ii][jj]-1e-9)
					&& (x[0] < data->mozaic.x1[ii][jj]+1e-9)
					&& (z[0] > data->mozaic.z0[ii][jj]-1e-9)
					&& (z[0] < data->mozaic.z1[ii][jj]+1e-9)
					&& (y[0] > data->mozaic.y0[ii][jj]-1e-9)
					&& (y[0] < data->mozaic.y1[ii][jj]+1e-9)){
						flag1 = 1;
					}

					if((x[1] > data->mozaic.x0[ii][jj]-1e-9)
					&& (x[1] < data->mozaic.x1[ii][jj]+1e-9)
					&& (z[1] > data->mozaic.z0[ii][jj]-1e-9)
					&& (z[1] < data->mozaic.z1[ii][jj]+1e-9)
					&& (y[1] > data->mozaic.y0[ii][jj]-1e-9)
					&& (y[1] < data->mozaic.y1[ii][jj]+1e-9)){
						flag2 = 1;
					}
					
					if((x[2] > data->mozaic.x0[ii][jj]-1e-9)
					&& (x[2] < data->mozaic.x1[ii][jj]+1e-9)
					&& (z[2] > data->mozaic.z0[ii][jj]-1e-9)
					&& (z[2] < data->mozaic.z1[ii][jj]+1e-9)
					&& (y[2] > data->mozaic.y0[ii][jj]-1e-9)
					&& (y[2] < data->mozaic.y1[ii][jj]+1e-9)){
						flag3 = 1;
					}
					
					if((x[3] > data->mozaic.x0[ii][jj]-1e-9)
					&& (x[3] < data->mozaic.x1[ii][jj]+1e-9)
					&& (z[3] > data->mozaic.z0[ii][jj]-1e-9)
					&& (z[3] < data->mozaic.z1[ii][jj]+1e-9)
					&& (y[3] > data->mozaic.y0[ii][jj]-1e-9)
					&& (y[3] < data->mozaic.y1[ii][jj]+1e-9)){
						flag4 = 1;
					}

					/*
					if((xg > data->mozaic.x0[ii][jj])
					&& (xg < data->mozaic.x1[ii][jj])
					&& (zg > data->mozaic.z0[ii][jj])
					&& (zg < data->mozaic.z1[ii][jj])
					&& (yg > data->mozaic.y0[ii][jj])
					&& (yg < data->mozaic.y1[ii][jj])){
					*/

					if(flag1 == 1 && flag2 == 1 && flag3 == 1 && flag4 == 1){
						data->mozaic.matID[k] = 1;
						data->mozaic.ii_pixel[k] = ii;
						data->mozaic.jj_pixel[k] = jj;
						
						data->mozaic.kk_pixel[ii][jj][data->mozaic.ne_pixel[ii][jj]] = k;
						data->mozaic.ne_pixel[ii][jj] += 1;
						
						if(data->mozaic.ne_pixel[ii][jj] == 500){
							fprintf(stderr, "ne of %d %d pixel exceeds 500!\n", ii, jj);
							// exit(0);
						}
					}
				}
			}

		} // if matID == mozaicID
	
	}// loop end of k

	fprintf(stderr, "The number of element in mozaic region = %d\n", tmp);

	// checking each element in mozaic region is counted properly

	tmp = 0;
	for (k = 0; k < fem->ne; k++) {

		for (long long int i = 0; i < fem->n_node; i++) {
			itmp = fem->element[k].kn[i];
			x[i] = fem->node[itmp].x;
			y[i] = fem->node[itmp].y;
			z[i] = fem->node[itmp].z;
		}

		matID = fem->element[k].matID;

		if(matID == data->mozaic.mozaicmatID 
			|| matID == data->mozaic.mozaicmatID+1){
			// search the pixel that the element belong
			xg = (x[0]+x[1]+x[2]+x[3])/4.0;
			yg = (y[0]+y[1]+y[2]+y[3])/4.0;
			zg = (z[0]+z[1]+z[2]+z[3])/4.0;

			if(data->mozaic.ii_pixel[k] == -1 || data->mozaic.jj_pixel[k] == -1){
				fprintf(stderr, "element No %d is not counted in mozaic region!\n", k);
				fprintf(stderr, "matID = %d\n", matID);

				flag1 = flag2 = flag3 = flag4 = 0;

				fprintf(stderr, "Vertex coordinate\n");
				fprintf(stderr, "falg1 = %d, (x0,y0,z0) = %lf %lf %lf\n", 
					flag1, x[0], y[0], z[0]);
				fprintf(stderr, "falg2 = %d, (x1,y1,z1) = %lf %lf %lf\n", 
					flag2, x[1], y[1], z[1]);
				fprintf(stderr, "falg3 = %d, (x2,y2,z2) = %lf %lf %lf\n", 
					flag3, x[2], y[2], z[2]);
				fprintf(stderr, "falg4 = %d, (x3,y3,z3) = %lf %lf %lf\n", 
					flag4, x[3], y[3], z[3]);

				tmp++;
			}
		}
	}

	if(tmp != 0){
		fprintf(stderr, "Non counted element No is %d\n", tmp);
		exit(0);
	}


	sprintf(string, "MozaicMesh.data");
	if ((fp = fopen(string, "w")) == NULL) {
		fprintf(stderr, "can't open file ( %s )\n", string);
		exit(1);
	}
	 
	for(ii = 0; ii < data->mozaic.Nx; ii++) {
		for(jj = 0; jj < data->mozaic.Nz; jj++) {		
			fprintf(fp, "(x,z) = (%d %d) pixel numofne %d\n", ii, jj, 
				data->mozaic.ne_pixel[ii][jj]);
			fprintf(fp, "kk: ");
			for(k = 0; k < data->mozaic.ne_pixel[ii][jj]; k++) {
				fprintf(fp, "%d ", data->mozaic.kk_pixel[ii][jj][k]);
			}
			fprintf(fp, "\n");
		}
	}

	fclose(fp);

	tmp = 0;
	// checking multiple count of element in different vocels
	for(ii = 0; ii < data->mozaic.Nx; ii++) {
		for(jj = 0; jj < data->mozaic.Nz; jj++) {	
			for(kk = 0; kk < data->mozaic.ne_pixel[ii][jj]; kk++) {

				for(l = 0; l < data->mozaic.Nx; l++) {
					for(m = 0; m < data->mozaic.Nz; m++) {	
						for(n = 0; n < data->mozaic.ne_pixel[l][m]; n++) {
							if(ii == l && jj == m && kk == n){

							} else {
								if(data->mozaic.kk_pixel[ii][jj][kk] 
									== data->mozaic.kk_pixel[l][m][n]){
									if(data->mozaic.kk_pixel[ii][jj][kk] 
										== data->mozaic.kk_pixel[l][m][n]){
										fprintf(stderr, "same element %d in different voxels!!\n",
											data->mozaic.kk_pixel[ii][jj][kk]);
										fprintf(stderr, "Voxel No xx%d zz%d\n", ii, jj);
										fprintf(stderr, "Voxel No xx%d zz%d\n", l, m);
										tmp++;
									}
								}
							}
						}
					}
				}

			}
		}
	}

	if(tmp != 0){
		fprintf(stderr, "tmp = %d, chek mesh again\n", tmp);
		exit(0);
	}
}


// Check the number of element in each voxel
// check the structured mesh

extern int checkMeshforeachVoxel(DataTable *data)
{

	int i, j, k, l, m, n, o;
	int ii, jj, kk, ll;
	int tmp = 0;

	int itmp;
	double x[10], y[10], z[10];
	double xg, yg, zg;
	double y0, y1;
	int matID;
	int flag1, flag2, flag3, flag4;

	int ***ne_pixel, ****kk_pixel;

	FEM *fem = &(data->fem);

	FILE *fp;
	char string[256];
	
	ALLOCATION(ne_pixel, int **, data->mozaic.Ny);
	ALLOCATION(kk_pixel, int ***, data->mozaic.Ny);

	for(kk = 0; kk < data->mozaic.Ny; kk++) {
	ALLOCATION(ne_pixel[kk], int *, data->mozaic.Nx);
	ALLOCATION(kk_pixel[kk], int **, data->mozaic.Nx);
	}

	for(kk = 0; kk < data->mozaic.Ny; kk++) {
	for(i = 0; i < data->mozaic.Nx; i++) {
		ALLOCATION(ne_pixel[kk][i], int, data->mozaic.Nz);
		ALLOCATION(kk_pixel[kk][i], int *, data->mozaic.Nz);
	}
	}

	for(kk = 0; kk < data->mozaic.Ny; kk++) {
	for(i = 0; i < data->mozaic.Nx; i++) {
		for(j = 0; j < data->mozaic.Nz; j++) {
	ne_pixel[kk][i][j] = 0;
	ALLOCATION(kk_pixel[kk][i][j], int, 100);
		}
	}
	}

	ll = 0;

	for (k = 0; k < fem->ne; k++) {

	for (long long int i = 0; i < fem->n_node; i++) {
		itmp = fem->element[k].kn[i];
		x[i] = fem->node[itmp].x;
		y[i] = fem->node[itmp].y;
		z[i] = fem->node[itmp].z;
	}

	matID = fem->element[k].matID;
	flag1 = flag2 = flag3 = flag4 = 0;

	if(matID == data->mozaic.mozaicmatID 
		 || matID == data->mozaic.mozaicmatID+1){
		// search the pixel that the element belong
		xg = (x[0]+x[1]+x[2]+x[3])/4.0;
		yg = (y[0]+y[1]+y[2]+y[3])/4.0;
		zg = (z[0]+z[1]+z[2]+z[3])/4.0;

		// fprintf(stderr, "xg,yg,zg = %lf %lf %lf\n", xg, yg, zg);

		tmp++;

		for(kk = 0; kk < data->mozaic.Ny; kk++) {

	y0 = data->mozaic.ys-(double)(kk+1)*data->mozaic.dy;
	y1 = data->mozaic.ys-(double)kk*data->mozaic.dy;

	for(ii = 0; ii < data->mozaic.Nx; ii++) {
		for(jj = 0; jj < data->mozaic.Nz; jj++) {

		flag1 = flag2 = flag3 = flag4 = 0;

		if((x[0] > data->mozaic.x0[ii][jj]-1e-9)
			 && (x[0] < data->mozaic.x1[ii][jj]+1e-9)
			 && (z[0] > data->mozaic.z0[ii][jj]-1e-9)
			 && (z[0] < data->mozaic.z1[ii][jj]+1e-9)
			 && (y[0] > y0-1e-9)
			 && (y[0] < y1+1e-9)){
			flag1 = 1;
		}

		if((x[1] > data->mozaic.x0[ii][jj]-1e-9)
			 && (x[1] < data->mozaic.x1[ii][jj]+1e-9)
			 && (z[1] > data->mozaic.z0[ii][jj]-1e-9)
			 && (z[1] < data->mozaic.z1[ii][jj]+1e-9)
			 && (y[1] > y0-1e-9)
			 && (y[1] < y1+1e-9)){
			flag2 = 1;
		}

		if((x[2] > data->mozaic.x0[ii][jj]-1e-9)
			 && (x[2] < data->mozaic.x1[ii][jj]+1e-9)
			 && (z[2] > data->mozaic.z0[ii][jj]-1e-9)
			 && (z[2] < data->mozaic.z1[ii][jj]+1e-9)
			 && (y[2] > y0-1e-9)
			 && (y[2] < y1+1e-9)){
			flag3 = 1;
		}

		if((x[3] > data->mozaic.x0[ii][jj]-1e-9)
			 && (x[3] < data->mozaic.x1[ii][jj]+1e-9)
			 && (z[3] > data->mozaic.z0[ii][jj]-1e-9)
			 && (z[3] < data->mozaic.z1[ii][jj]+1e-9)
			 && (y[3] > y0-1e-9)
			 && (y[3] < y1+1e-9)){
			flag4 = 1;
		}

		/*
		if((xg > data->mozaic.x0[ii][jj])
			 && (xg < data->mozaic.x1[ii][jj])
			 && (zg > data->mozaic.z0[ii][jj])
			 && (zg < data->mozaic.z1[ii][jj])
			 && (yg > y0)
			 && (yg < y1)){
		*/

		/*
		if(k == 118718 
			 && (flag1 != 0 || flag2 != 0 || flag3 != 0 || flag4 != 0)){
			fprintf(stderr, "(x0,x1) = %15.10lf %15.10lf\n", 
				data->mozaic.x0[ii][jj]-1e-9, 
				data->mozaic.x1[ii][jj]+1e-9);
			fprintf(stderr, "(y0,y1) = %15.10lf %15.10lf\n", 
				y0-1e-9, 
				y1+1e-9);

			fprintf(stderr, "(z0,z1) = %15.10lf %15.10lf\n", 
				data->mozaic.z0[ii][jj]-1e-9, 
				data->mozaic.z1[ii][jj]+1e-9);

			fprintf(stderr, "(x0,y0,z0) = %15.10lf %15.10lf %15.10lf\n", 
				x[0], y[0], z[0]);
			fprintf(stderr, "(x1,y1,z1) = %15.10lf %15.10lf %15.10lf\n", 
				x[1], y[1], z[1]);
			fprintf(stderr, "(x2,y2,z2) = %15.10lf %15.10lf %15.10lf\n", 
				x[2], y[2], z[2]);
			fprintf(stderr, "(x3,y3,z3) = %15.10lf %15.10lf %15.10lf\n", 
				x[3], y[3], z[3]);
			fprintf(stderr, "%d %d %d %d\n", flag1, flag2, flag3, flag4);

			// exit(0);
		}
		*/

		if(flag1 == 1 && flag2 == 1 && flag3 == 1 && flag4 == 1){
			kk_pixel[kk][ii][jj][ne_pixel[kk][ii][jj]] = k;
			ne_pixel[kk][ii][jj] += 1;
			
			if(ne_pixel[kk][ii][jj] == 100){
		fprintf(stderr, "ne of %d %d pixel exceeds 100!\n", ii, jj);
		// exit(0);
			}
			
		}
	
		}
	}

		}

	}
	
	}// loop end of k

	fprintf(stderr, "The number of element in mozaic region = %d\n", tmp);

	sprintf(string, "MozaicMesh.data");
	if ((fp = fopen(string, "w")) == NULL) {
	fprintf(stderr, "can't open file ( %s )\n", string);
	exit(1);
	}
	 
	for(kk = 0; kk < data->mozaic.Ny; kk++) {
	for(ii = 0; ii < data->mozaic.Nx; ii++) {
		for(jj = 0; jj < data->mozaic.Nz; jj++) {
	
	fprintf(fp, "(y,x,z) = (%d %d %d) pixel numofne %d\n", kk, ii, jj, 
		ne_pixel[kk][ii][jj]);
	fprintf(fp, "kk: ");
	for(k = 0; k < ne_pixel[kk][ii][jj]; k++) {
		fprintf(fp, "%d ", kk_pixel[kk][ii][jj][k]);
	}
	fprintf(fp, "\n");
		}
	}
	}

	fclose(fp);
	
	// checking multiple count of element in different voxels
	tmp = 0;
	for(kk = 0; kk < data->mozaic.Ny; kk++) {
	for(ii = 0; ii < data->mozaic.Nx; ii++) {
		for(jj = 0; jj < data->mozaic.Nz; jj++) {	
	for(ll = 0; ll < ne_pixel[kk][ii][jj]; ll++) {

		for(k = 0; k < data->mozaic.Ny; k++) {
		for(i = 0; i < data->mozaic.Nx; i++) {
			for(j = 0; j < data->mozaic.Nz; j++) {	
		for(l = 0; l < ne_pixel[k][i][j]; l++) {
			
			if(ii == i && jj == j && kk == k && ll == l){
			
			}else{
			if(kk_pixel[kk][ii][jj][ll] == kk_pixel[k][i][j][l]){
				fprintf(stderr, "same element %d in different voxels!!\n",
					kk_pixel[kk][ii][jj][ll]);
				fprintf(stderr, "Voxel No yy%d xx%d zz%d\n", kk, ii, jj);
				fprintf(stderr, "Voxel No yy%d xx%d zz%d\n", k, i, j);
				tmp++;
			}
			}
			
		}
			}
		}
		}
		
	}
		}
	}
	}

	if(tmp != 0){
	fprintf(stderr, "tmp = %d, chek mesh again\n", tmp);
	}

	// exit(0);

}

extern int shuffle(int *array, int size) {
	for(int i = 0; i < size; i++) {
		//int j = rand()%size;
		int j = (int)(randomperc()*(double)RAND_MAX)%size;

		int t = array[i];
		array[i] = array[j];
		array[j] = t;
	}
}
