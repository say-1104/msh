#include "optfem.h"

// This file includes fucntions for interpolation
// interpolation is used for calculating field for output file

// NX, NY
// the number of division for x and y directions
// these variables are used for fast interpolation and emerge in Preforintpol
// see note for detail
#define NX           50
#define NY           50

// function for preparing element information for fast interpolation
// see note for detail
extern int Preforintpol(DataTable *data)
{
	int i, j, k, l;
	int *element, flag = 0, count = 0, e_number, itmp;
	double xmax, xmin, ymax, ymin;
	double e_xmax, e_xmin, e_ymax, e_ymin;
	int start_x, start_y, end_x, end_y;
	double Dx, Dy, Xmax, Xmin, Ymax, Ymin;
	double x[6], y[6];

	xmin = data->fem.xmin;
	xmax = data->fem.xmax;
	ymin = data->fem.ymin;
	ymax = data->fem.ymax;

	ALLOCATION(data->fem.table_n, int *, NX)
	ALLOCATION(data->fem.table, int **, NX)

	for(i = 0; i < NX; i++) {
		ALLOCATION(data->fem.table_n[i], int, NY)
		ALLOCATION(data->fem.table[i], int *, NY)
	}

	for (i = 0; i < NX; i++)
		for (j = 0; j < NY; j++) data->fem.table_n[i][j] = 0;

	fprintf(stderr, "before pre\n");

	fprintf(stderr, "%lf %lf\n", xmin, xmax);
	fprintf(stderr, "%lf %lf\n", ymin, ymax);

	for(k = 0; k < data->fem.ne; k++) {
		e_xmin = xmax; e_xmax = xmin;
		e_ymin = ymax; e_ymax = ymin;
		for (i = 0; i < 3; i++) {
			itmp = data->fem.element[k].kk[i];
			x[i] = data->fem.node[itmp].x;
			y[i] = data->fem.node[itmp].y;
			if (e_xmin > x[i]) e_xmin = x[i];
			if (e_xmax < x[i]) e_xmax = x[i];
			if (e_ymin > y[i]) e_ymin = y[i];
			if (e_ymax < y[i]) e_ymax = y[i];
		}
		start_x = (int)((e_xmin-xmin)/(xmax-xmin)*NX);
		start_y = (int)((e_ymin-ymin)/(ymax-ymin)*NY);
		end_x = (int)((e_xmax-xmin)/(xmax-xmin)*NX)+1;
		end_y = (int)((e_ymax-ymin)/(ymax-ymin)*NY)+1;
		if (end_x > NX) end_x = NX;
		if (end_y > NY) end_y = NY;

		for (i = start_x; i < end_x; i++) {
			for (j = start_y; j < end_y; j++) {
				data->fem.table_n[i][j]++;
			}
		}
	}

	for(i = 0; i < NX; i++) {
		for(j = 0; j < NY; j++) {
			ALLOCATION(data->fem.table[i][j], int, data->fem.table_n[i][j])
			data->fem.table_n[i][j] = 0;
		}
	}

	for(k = 0; k < data->fem.ne; k++) {
		e_xmin = xmax; e_xmax = xmin;
		e_ymin = ymax; e_ymax = ymin;
		for (i = 0; i < 3; i++) {
			itmp = data->fem.element[k].kk[i];
			x[i] = data->fem.node[itmp].x;
			y[i] = data->fem.node[itmp].y;
			if (e_xmin > x[i]) e_xmin = x[i];
			if (e_xmax < x[i]) e_xmax = x[i];
			if (e_ymin > y[i]) e_ymin = y[i];
			if (e_ymax < y[i]) e_ymax = y[i];
		}
		start_x = (int)((e_xmin-xmin)/(xmax-xmin)*NX);
		start_y = (int)((e_ymin-ymin)/(ymax-ymin)*NY);
		end_x = (int)((e_xmax-xmin)/(xmax-xmin)*NX)+1;
		end_y = (int)((e_ymax-ymin)/(ymax-ymin)*NY)+1;
		if (end_x > NX) end_x = NX;
		if (end_y > NY) end_y = NY;

		for (i = start_x; i < end_x; i++) {
			for (j = start_y; j < end_y; j++) {
				data->fem.table[i][j][data->fem.table_n[i][j]] = k;
				data->fem.table_n[i][j]++;
			}
		}
	}

	fprintf(stderr, "after pre\n");
}


// function for calculating Ex, Ey, Ez, Hx, Hy, and Hz
// for designated coordinate (x0, y0)
extern int intpolFieldValue(DataTable *data, double x0, double y0,
					double xmin, double xmax, double ymin, double ymax,
					double_complex phi[6])
{
	double_complex phix, phiy, phiz;

	double x1, x2, x3, y1, y2, y3;
	double x[6], y[6];
	double a1, a2, a3, b1, b2, b3, c1, c2, c3;
	double area, L1, L2, L3, al0, al1, al2;
	double j11, j12, j21, j22;
	double L1x, L2x, L3x, L1y, L2y, L3y;
	double s12, s23, s31;
	double N[6], Nx[6], Ny[6], nx[6], ny[6];
	double_complex vr[6];

	int i, j, k, m, ii, jj, flag, iel;
	int tmp1, tmp2, tmp3;
	int region_x, region_y;

	double k0;
	k0 = data->par.k0;

	FEM *fem = &(data->fem);


	flag = 0;
	m = 0;
	phix = phiy = phiz = 0.0;

	region_x = (int)((x0-xmin)/(xmax-xmin)*NX);
	region_y = (int)((y0-ymin)/(ymax-ymin)*NY);
	if (region_x >= NX) region_x = NX-1;
	if (region_y >= NY) region_y = NY-1;

	for (ii = 0; ii < data->fem.table_n[region_x][region_y]; ii++) {
		iel = data->fem.table[region_x][region_y][ii];
		tmp1 = fem->element[iel].kk[0];
		tmp2 = fem->element[iel].kk[1];
		tmp3 = fem->element[iel].kk[2];

		x1 = fem->node[tmp1].x;
		x2 = fem->node[tmp2].x;
		x3 = fem->node[tmp3].x;
		y1 = fem->node[tmp1].y;
		y2 = fem->node[tmp2].y;
		y3 = fem->node[tmp3].y;

		a1 = x2*y3-x3*y2, a2 = x3*y1-x1*y3, a3 = x1*y2-x2*y1;
		b1 = y2-y3, b2 = y3-y1, b3 = y1-y2;
		c1 = x3-x2, c2 = x1-x3, c3 = x2-x1;
		area = x1*b1+x2*b2+x3*b3;

		for (j = 0; j < 6; j++) {
			tmp1 = fem->element[iel].kk[j];
			x[j] = fem->node[tmp1].x;
			y[j] = fem->node[tmp1].y;
		}

		L1 = (a1+b1*x0+c1*y0)/area;
		L2 = (a2+b2*x0+c2*y0)/area;
		L3 = (a3+b3*x0+c3*y0)/area;

		if ((L1 >= ZERO) && (L2 >= ZERO) && (L3 >= ZERO) && (m == 0)) {
			m = 1;
			al0 = sqrt(b3*b3 + c3*c3);
			al1 = sqrt(b1*b1 + c1*c1);
			al2 = sqrt(b2*b2 + c2*c2);

			for (jj = 0 ; jj < 6; jj++) {
				vr[jj] = data->fem.phi[fem->element[iel].kk[jj]];
			}

			// shape function for nodal element
			// second order
			N[0] = L1*(2.0*L1-1.0);  nx[0] = 4.0*L1-1.0;    ny[0] = 0.0;
			N[1] = L2*(2.0*L2-1.0);  nx[1] = 0.0;           ny[1] = 4.0*L2-1.0;
			N[2] = L3*(2.0*L3-1.0);  nx[2] = 1.0-4.0*L3;    ny[2] = 1.0-4.0*L3;
			N[3] = 4.0*L1*L2;        nx[3] = 4.0*L2;        ny[3] = 4.0*L1;
			N[4] = 4.0*L2*L3;        nx[4] = -4.0*L2;       ny[4] = 4.0*(L3-L2);
			N[5] = 4.0*L3*L1;        nx[5] = 4.0*(L3-L1);   ny[5] = -4.0*L1;

			for(j11 = 0.0, j = 0; j < 6; j++) j11 += nx[j]*x[j];
			for(j12 = 0.0, j = 0; j < 6; j++) j12 += nx[j]*y[j];
			for(j21 = 0.0, j = 0; j < 6; j++) j21 += ny[j]*x[j];
			for(j22 = 0.0, j = 0; j < 6; j++) j22 += ny[j]*y[j];

			area = j11*j22-j12*j21;
			L1x = j22/area;
			L2x = -j12/area;
			L3x = (-j22+j12)/area;
			L1y = -j21/area;
			L2y = j11/area;
			L3y = (j21-j11)/area;

			for (j = 0; j < 6; j++) {
				Nx[j] = ( j22*nx[j]-j12*ny[j])/area;
				Ny[j] = (-j21*nx[j]+j11*ny[j])/area;
			}

			// -------------------------------------------------
			// Electric field

			phix = phiy = phiz = 0.0;
			for (jj = 0 ; jj < 6; ++jj) {
				phix += N[jj]*vr[jj];
			}

			phi[0] = (phix);

			// Magnetic field
			phix = phiy = phiz = 0.0;
			for (jj = 0 ; jj < 6; ++jj) {
			//	phix += (data->par.beta/(k0*Z0))*(Ny[jj]*vr[jj]);
			//	phiy += (data->par.beta/(k0*Z0))*(Nx[jj]*vr[jj]);
			}

			phi[1] = (phix);
			phi[2] = (phiy);

			flag = 1;
			break;
		}
	}

	if (flag == 0) {
		phi[0] = 0.0;
		phi[1] = 0.0;
		phi[2] = 0.0;
		phi[3] = 0.0;
		phi[4] = 0.0;
		phi[5] = 0.0;
	}
}


// function for interpolating output field in xy plane
extern int intpol(DataTable *data, double w_c, int mm)
{
	int DIV_X, DIV_Y;

	DIV_X = data->par.div_x;
	DIV_Y = data->par.div_y;

	int i, j, tmp, tmp2;

	double xmin, xmax, ymin, ymax;
	double xdiv, ydiv;
	double x0, y0;
	double_complex phi[6];
	double **Ex, **Ey, **Ez;
	double **Hx, **Hy, **Hz;
	ALLOCATION(Ex, double *, DIV_X)
	ALLOCATION(Ey, double *, DIV_X)
	ALLOCATION(Ez, double *, DIV_X)
	ALLOCATION(Hx, double *, DIV_X)
	ALLOCATION(Hy, double *, DIV_X)
	ALLOCATION(Hz, double *, DIV_X)
	for (i = 0; i < DIV_X; i++) {
		ALLOCATION(Ex[i], double, DIV_Y)
		ALLOCATION(Ey[i], double, DIV_Y)
		ALLOCATION(Ez[i], double, DIV_Y)
		ALLOCATION(Hx[i], double, DIV_Y)
		ALLOCATION(Hy[i], double, DIV_Y)
		ALLOCATION(Hz[i], double, DIV_Y)
	}

	double *xx, *yy;
	ALLOCATION(xx, double, DIV_X)
	ALLOCATION(yy, double, DIV_Y)

	FEM *fem = &(data->fem);
	FILE *fp, *fq, *fr, *fs, *ft, *fu, *fv, *fw;
	char file[256];

	xmin = data->fem.xmin;
	xmax = data->fem.xmax;
	ymin = data->fem.ymin;
	ymax = data->fem.ymax;

	xdiv = (xmax-xmin)/DIV_X;
	ydiv = (ymax-ymin)/DIV_Y;
	sprintf(file, "Output/Field_%d-%1.3f.slv\0", mm, w_c);
	fp = fopen(file, "w");

	/*
	sprintf(file, "Ey-%1.3f.slv\0", w_c);
	fq = fopen(file, "w");

	sprintf(file, "Ez-%1.3f.slv\0", w_c);
	fr = fopen(file, "w");

	sprintf(file, "Hx-%1.3f.slv\0", w_c);
	fs = fopen(file, "w");

	sprintf(file, "Hy-%1.3f.slv\0", w_c);
	ft = fopen(file, "w");

	sprintf(file, "Hz-%1.3f.slv\0", w_c);
	fu = fopen(file, "w");
	*/

	sprintf(file, "Output/XX");
	fv = fopen(file, "w");

	sprintf(file, "Output/YY");
	fw = fopen(file, "w");

	// set coordinate for interpolation
	for (i = 0; i < DIV_X; i++) {
		xx[i] = (xmin+0.5*xdiv)+(xdiv*(double)i);
		//    fprintf(fs, "%lf\n", x0);
	}

	for (j = 0; j < DIV_Y; j++) {
		yy[j] = (ymin+0.5*ydiv)+(ydiv*(double)j);
		//    fprintf(ft, "%lf\n", y0);
	}

	// interpolation of field value
	for (i = 0; i < DIV_X; i++) {
		for (j = 0; j < DIV_Y; j++) {

		x0 = xx[i];
		y0 = yy[j];

		intpolFieldValue(data, x0, y0, xmin, xmax, ymin, ymax, phi);

			if (data->realflag == 0) {
				Ex[i][j] = abs(phi[0]);
			}
			else {
				Ex[i][j] = real(phi[0]);
			}
			// *added field interpolation(intpolFieldValue) other than Ex has not been implemented yet.
		}
	}

	for (i = 0; i < DIV_X; i++) {
		fprintf(fv, "%15.10lf\n", xx[i]);
	}

	for (j = 0; j < DIV_Y; j++) {
		fprintf(fw, "%15.10lf\n", yy[j]);
	}
	fclose(fv);
	fclose(fw);

	for (i = 0; i < DIV_X; i++) {
		for (j = 0; j < DIV_Y; j++) {
			fprintf(fp, "%15.10lf ", Ex[i][j]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);


	for (i = 0; i < DIV_X; i++) {
		free(Ex[i]);
		free(Ey[i]);
		free(Ez[i]);
		free(Hx[i]);
		free(Hy[i]);
		free(Hz[i]);
	}
	free(Ex);
	free(Ey);
	free(Ez);
	free(Hx);
	free(Hy);
	free(Hz);

	free(xx);
	free(yy);
}

// function for interpolating x cross section of the field
extern int intpolXX(DataTable *data, double wl, int mm)
{
	int DIVISION_X = data->par.division_x;

	int i, j;
	double xmin, xmax, ymin, ymax;
	double xdiv, ydiv, xcenter, ycenter;
	double x0, y0;
	double_complex phi[6];

	FEM *fem = &(data->fem);
	char file[256];
	FILE *fp;

	for (j = 0; j < data->par.numx; j++) {
		sprintf(file, "Output/fieldX_%d-%d-%1.3f\0", mm, j, wl);
		fp = fopen(file, "w");
		fclose(fp);
	}

	xmin = data->fem.xmin;
	xmax = data->fem.xmax;
	ymin = data->fem.ymin;
	ymax = data->fem.ymax;

	xdiv = (xmax-xmin)/DIVISION_X;

	for (j = 0; j < data->par.numx; j++) {
		sprintf(file, "Output/fieldX_%d-%d-%1.3f\0", mm, j, wl);
		fp = fopen(file, "a+");
		for (i = 0; i <= DIVISION_X; i++) {
			x0 = xmin+(xdiv*(double)i);
			y0 = data->par.center_y[j];

			intpolFieldValue(data, x0, y0, xmin, xmax, ymin, ymax, phi);

			if (data->realflag == 0) {
				fprintf(fp, "%10.6e %10.6e\n", x0, abs(phi[0]));
			}
			else {
				fprintf(fp, "%10.6e %10.6e\n", x0, real(phi[0]));
			}
		}
		fclose(fp);
	}
}

// function for interpolating y cross section of the field
extern int intpolYY(DataTable *data, double wl, int mm)
{
	int DIVISION_Y = data->par.division_y;

	int i, j;
	double xmin, xmax, ymin, ymax;
	double xdiv, ydiv, xcenter, ycenter;
	double x0, y0;
	double_complex phi[6];

	FEM *fem = &(data->fem);
	char file[256];
	FILE *fp;

	for (j = 0; j < data->par.numy; j++) {
		sprintf(file, "Output/fieldY_%d-%d-%1.3f\0", mm, j, wl);
		fp = fopen(file, "w");
		fclose(fp);
	}

	xmin = data->fem.xmin;
	xmax = data->fem.xmax;
	ymin = data->fem.ymin;
	ymax = data->fem.ymax;

	ydiv = (ymax-ymin)/DIVISION_Y;

	for (j = 0; j < data->par.numy; j++) {
		sprintf(file, "Output/fieldY_%d-%d-%1.3f\0", mm, j, wl);
		fp = fopen(file, "a+");
		for (i = 0; i < DIVISION_Y; i++) {

			x0 = data->par.center_x[j];
			//      y0 = ymin+(ydiv*(double)i);
			y0 = ymin+0.5*ydiv+(ydiv*(double)i);

			intpolFieldValue(data, x0, y0, xmin, xmax, ymin, ymax, phi);

			if (data->realflag == 0) {
				fprintf(fp, "%10.6e %10.6e\n", y0, abs(phi[0]));
			}
			else {
				fprintf(fp, "%10.6e %10.6e\n", y0, real(phi[0]));
			}
		}
		fclose(fp);
	}
}