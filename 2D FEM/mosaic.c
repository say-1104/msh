#include "optfem.h"

// reading mosaic data, memory allocation

extern int inputMosaic(DataTable *data)
{
	int i, j, k, jj;
	char filename[256];
	FILE *fp;
	
	data->mosaic.DBSitr = 5;


	sprintf(filename, "mosaic.pre");
	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "can't open file ( %s )\n", filename);
		exit(1);
	}

	fscanf(fp, "%*s %lf %lf", &data->mosaic.zs, &data->mosaic.ys);
	fprintf(stderr,"Starting z, y (left upper position) = %lf, %lf\n", 
		data->mosaic.zs, data->mosaic.ys);

	fscanf(fp, "%*s %d %d", &data->mosaic.Nz, &data->mosaic.Ny);
	fprintf(stderr,"The number of mosaic region in z,y = %d, %d\n", 
		data->mosaic.Nz, data->mosaic.Ny);

	fscanf(fp, "%*s %lf %lf", &data->mosaic.dz, &data->mosaic.dy);
	fprintf(stderr,"dz, dy = %lf, %lf\n", data->mosaic.dz, data->mosaic.dy);

	fscanf(fp, "%*s %d", &data->mosaic.symflag);
	fscanf(fp, "%*s %d", &data->mosaic.mosaicmatID);
	fscanf(fp, "%*s %d", &data->mosaic.coreID);
	fscanf(fp, "%*s %d", &data->mosaic.cladID);
	fscanf(fp, "%*s %lf", &data->mosaic.gamma);

	fscanf(fp, "%*s %d", &data->mosaic.DBSitr);
	fscanf(fp, "%*s %lf", &data->mosaic.termination);
	fscanf(fp, "%*s %d", &data->mosaic.PixelReduction);

	fprintf(stderr,"symflag = %d\n", data->mosaic.symflag);

	fprintf(stderr,"mosaicID, coreID, cladID = %d, %d, %d\n",
		data->mosaic.mosaicmatID, data->mosaic.coreID,
		data->mosaic.cladID);

	if (MAX(data->mosaic.cladID, MAX(data->mosaic.mosaicmatID, data->mosaic.coreID)) >= data->par.n_material 
	||	MIN(data->mosaic.cladID, MIN(data->mosaic.mosaicmatID, data->mosaic.coreID)) < 0) {
		fprintf(stderr, "Material No. for mosaic is improper");
		exit(1);
	}

	fprintf(stderr,"Max iteration no. = %d\n", data->mosaic.DBSitr);
	fprintf(stderr,"termination condition = %lf\n", data->mosaic.termination);
	if (data->mosaic.PixelReduction != 0) fprintf(stderr,"pixels of -1 are not considered\n");

	ALLOCATION(data->mosaic.BorW, int *, data->mosaic.Ny);
	for(j = 0; j < data->mosaic.Ny; j++) {
		ALLOCATION(data->mosaic.BorW[j], int, data->mosaic.Nz);
	}

	ALLOCATION(data->mosaic.z0, double *, data->mosaic.Ny);
	ALLOCATION(data->mosaic.z1, double *, data->mosaic.Ny);
	ALLOCATION(data->mosaic.y0, double *, data->mosaic.Ny);
	ALLOCATION(data->mosaic.y1, double *, data->mosaic.Ny);

	ALLOCATION(data->mosaic.zz, double, data->mosaic.Nz);
	ALLOCATION(data->mosaic.yy, double, data->mosaic.Ny);

	for(i = 0; i < data->mosaic.Ny; i++) {
		ALLOCATION(data->mosaic.z0[i], double, data->mosaic.Nz);
		ALLOCATION(data->mosaic.z1[i], double, data->mosaic.Nz);
		ALLOCATION(data->mosaic.y0[i], double, data->mosaic.Nz);
		ALLOCATION(data->mosaic.y1[i], double, data->mosaic.Nz);
	}	

	fclose(fp);
	// checking minimum and maximum x,y coordinate in each pixel
	for(i = 0; i < data->mosaic.Ny; i++) {
		for(j = 0; j < data->mosaic.Nz; j++) {
			data->mosaic.z0[i][j] = (double)j*data->mosaic.dz+data->mosaic.zs;
			data->mosaic.z1[i][j] = (double)(j+1)*data->mosaic.dz+data->mosaic.zs;
			data->mosaic.y0[i][j] = data->mosaic.ys-(double)(i+1)*data->mosaic.dy;
			data->mosaic.y1[i][j] = data->mosaic.ys-(double)(i)*data->mosaic.dy;
		}
	}

	// y, z coordinate of the center of each pixel
	sprintf(filename, "./Output/YY-Mosaic");
	if ((fp = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "can't open file ( %s )\n", filename);
		exit(1);
	}

	for(i = 0; i < data->mosaic.Ny; i++) {
		data->mosaic.yy[i] 
			= data->mosaic.ys-0.5*data->mosaic.dy-(double)(i)*data->mosaic.dy;
		fprintf(fp, "%lf\n", data->mosaic.yy[i]);
	}
	fclose(fp);

	sprintf(filename, "Output/ZZ-Mosaic");
	if ((fp = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "can't open file ( %s )\n", filename);
		exit(1);
	}

	for(i = 0; i < data->mosaic.Nz; i++) {
		data->mosaic.zz[i] 
			= data->mosaic.zs+0.5*data->mosaic.dz+(double)(i)*data->mosaic.dz;
		fprintf(fp, "%lf\n", data->mosaic.zz[i]);
	}
	fclose(fp);

	if (data->mosaic.mzcalcNo != -2) {
		sprintf(filename, "MosaicPattern.csv");
		if ((fp = fopen(filename, "r")) == NULL) {
			fprintf(stderr, "can't open file ( %s )\n", filename);
			exit(1);
		}
		
		for(i = 0; i < data->mosaic.Ny; i++) {
			for(j = 0; j < data->mosaic.Nz; j++) {
				if (i == data->mosaic.Ny-1 && j == data->mosaic.Nz-1) {
					fscanf(fp, "%d", &data->mosaic.BorW[i][j]);
				}
				else {
					fscanf(fp, "%d,", &data->mosaic.BorW[i][j]);
				}
			}
		}
	}
	else {
		// Guessed matrix pattern		
		sprintf(filename, "GuessedMatrix.csv");
		if ((fp = fopen(filename, "r")) == NULL) {
			fprintf(stderr, "can't open file ( %s )\n", filename);
			exit(1);
		}
		
		for(i = 0; i < data->mosaic.Ny; i++) {
			for(j = 0; j < data->mosaic.Nz; j++) {
				if (i == data->mosaic.Ny-1 && j == data->mosaic.Nz-1) {
					fscanf(fp, "%d", &data->mosaic.BorW[i][j]);
				}
				else {
					fscanf(fp, "%d,", &data->mosaic.BorW[i][j]);
				}
			}
		}

		for(i = 0; i < data->mosaic.Ny; i++) {
			for(j = 0; j < data->mosaic.Nz; j++) {
				printf("%d ", data->mosaic.BorW[i][j]);
			}
			printf("\n");
		}
		printf("\n");
	}

	// parameters for pixel
	ALLOCATION(data->mosaic.matID, int, data->fem.ne);
	ALLOCATION(data->mosaic.ii_pixel, int, data->fem.ne);
	ALLOCATION(data->mosaic.jj_pixel, int, data->fem.ne);

	for(k = 0; k < data->fem.ne; k++) {
		data->mosaic.matID[k] = -1;
		data->mosaic.ii_pixel[k] = -1;
		data->mosaic.jj_pixel[k] = -1;
	}
	
	ALLOCATION(data->mosaic.ne_pixel, int *, data->mosaic.Ny);
	ALLOCATION(data->mosaic.kk_pixel, int **, data->mosaic.Ny);

	for(i = 0; i < data->mosaic.Ny; i++) {
		ALLOCATION(data->mosaic.ne_pixel[i], int, data->mosaic.Nz);
		ALLOCATION(data->mosaic.kk_pixel[i], int *, data->mosaic.Nz);
	}

	for(i = 0; i < data->mosaic.Ny; i++) {
		for(j = 0; j < data->mosaic.Nz; j++) {
			data->mosaic.ne_pixel[i][j] = 0;
			ALLOCATION(data->mosaic.kk_pixel[i][j], int, 300);
		}
	}
}

// Check mesh data
extern int checkMeshforMosaic(DataTable *data)
{
	int		 i, j, k;
	int		 ii, jj, kk;
	int		 tmp = 0, tmp1 = 0;
	double	x[6], y[6];
	double	zg, yg;
	int		 matID;

	FEM		 *fem = &(data->fem);
	
	for (k = 0; k < fem->ne; k++) {
		for (j = 0; j < 6; j++) {
			i = fem->element[k].kk[j];
			x[j] = fem->node[i].x;
			y[j] = fem->node[i].y;
		}
		matID = fem->element[k].matID;

		if (matID == data->mosaic.mosaicmatID) {
			// search the pixel that the element belong
			zg = (x[0]+x[1]+x[2])/3.0;
			yg = (y[0]+y[1]+y[2])/3.0;

			tmp++;

			for(ii = 0; ii < data->mosaic.Ny; ii++) {
				for(jj = 0; jj < data->mosaic.Nz; jj++) {
					if ((zg > data->mosaic.z0[ii][jj]-1.0e-10) && (zg < data->mosaic.z1[ii][jj]+1.0e-10)
					&& (yg > data->mosaic.y0[ii][jj]-1.0e-10) && (yg < data->mosaic.y1[ii][jj]+1.0e-10)) {

					data->mosaic.matID[k] = 1;
					data->mosaic.ii_pixel[k] = ii;
					data->mosaic.jj_pixel[k] = jj;

					data->mosaic.kk_pixel[ii][jj][data->mosaic.ne_pixel[ii][jj]] = k;
					data->mosaic.ne_pixel[ii][jj] += 1;
					
					if (data->mosaic.ne_pixel[ii][jj] == 300) {
						fprintf(stderr, "ne of %d %d pixel exceeds 300!\n", ii, jj);
						exit(0);
					}
					tmp1++;
					}
				}
			}
		}
	}

	fprintf(stderr, "tmp = %d\n", tmp);
}

extern int shuffle(int *array, int size) {
	srand((unsigned)time(NULL));

	for(int i = 0; i < size; i++) {
		int j = rand()%size;
		int t = array[i];
		array[i] = array[j];
		array[j] = t;
	}

}
