#include "optfem.h"

extern int checkArgument1(DataTable *data, int argc, char **argv)
{
	int i, j;
	double R = 1.0e-4;
	
	data->port.inputDirection = BOTH_DIRECTION;
	data->port.forward = 1.0;
	data->port.backward = 1.0;
	data->par.reflect = -log(R);
	data->par.mo_effect = OFF;
	data->par.PMLtype = 1;
	data->par.inputPort = 0;
	data->par.modEIMflag = 0;

	data->mosaic.Device = DEVICE_NONE;
	data->mosaic.DBSflag = -1;
	data->mosaic.mzcalcNo = -1;
	data->mosaic.evaluation = 0;

	for (i = 0; i < argc; i++) {
		if (argv[i][0] == '-') {
			if (argv[i][1] == 'm' && argv[i][2] == 'z') {	
				data->mosaic.mosflag = 1;
				fprintf(stderr, "Mosaic pattern simulation\n");
			}
			else if (argv[i][1] == 'f' && argv[i][2] == 'm' && argv[i][3] == 'z') {
				sscanf(argv[i]+4, "%d", &data->mosaic.mzcalcNo);
				fprintf(stderr, "%d th mosaic structure is calculated\n", 
					data->mosaic.mzcalcNo);
			}
			else if (argv[i][1] == 'D' && argv[i][2] == 'B' && argv[i][3] == 'S') {
				sscanf(argv[i]+4, "%d", &data->mosaic.DBSflag);
				fprintf(stderr, "Direct Binary Serch optimization\n");
				fprintf(stderr, "%d th mosaic structure is used\n", 
					data->mosaic.DBSflag);
			}
			else if (argv[i][1] == 'D' && argv[i][2] == 'I' && argv[i][3] == 'V') {
				fprintf(stderr, "modeDivider is optimized\n");
				data->mosaic.Device = MODE_DIVIDER;
			}
			else if (argv[i][1] == 'M' && argv[i][2] == 'U' && argv[i][3] == 'X') {
				fprintf(stderr, "modeMUX is optimized\n");
				data->mosaic.Device = MODE_MUX;
			}
			else if (argv[i][1] == 'P' && argv[i][2] == 'S') {
				fprintf(stderr, "Power splitter is optimized\n");
				data->mosaic.Device = POWER_SPLITTER;
				if (argv[i][3] >= '0' && argv[i][3] <= '9') {
					data->mosaic.evaluation = argv[i][3]-'0';
				}
				
				if (data->mosaic.evaluation == 1) {
					fprintf(stderr, "Euclidean distance\n");
				}
				else if (data->mosaic.evaluation == 2) {
					fprintf(stderr, "City block distance\n");
				}
				else if (data->mosaic.evaluation == 3) {
					fprintf(stderr, "Chess-board distance\n");
				}
				else {
					fprintf(stderr, "Manual set\n");
				}
			}
			else if (argv[i][1] == 'W' && argv[i][2] == 'C') {
				fprintf(stderr, "Waveguide Crossing is optimized\n");
				data->mosaic.Device = WAVEGUIDE_CROSSING;
			}		
		} 
	}
}

extern int checkArgument2(DataTable *data, int argc, char **argv)
{
	int i, j, ll, mm, n_mode, numoflambda;
	double R = 1.0e-4;
	
	n_mode = data->par.num_mode;
	numoflambda = data->par.numoflambda;

	data->port.solver = SOLVER_NONE;
	for (ll = 0; ll < numoflambda; ll++) {
		for (mm = 0; mm < n_mode; mm++) {
			for (i = 0; i < data->port.number; i++) {
				sprintf(data->port.dataMemory[ll][mm][i].file, "port%d_%d-%d", i, mm, ll);
			}
		}
	}

	data->realflag = 0;
	
	for (i = 0; i < argc; i++) {
		if (argv[i][0] == '-') {
			if (argv[i][1] == 'l') {
				data->par.w_start = atof(argv[i]+2);
				fprintf(stderr, "modified wavelength = %lf\n", data->par.w_start);
			}
			if (argv[i][1] == 'f') {
				// file input
				data->port.solver = INPUTEIGEN;
			}
			
			if (argv[i][1] == 'r') {
				data->realflag = 1;
			}

			if (argv[i][1] == 'p' && argv[i][2] == 't') {
				// pointsource
				sscanf(argv[i]+3, "%lf,%lf", &data->port.x0, &data->port.y0);
				fprintf(stderr, "Point source center x,y = %lf, %lf\n", data->port.x0, data->port.y0);
				data->port.solver = POINTSOURCE;
				data->realflag = 1;
				data->par.num_mode = 1;
			}

			if (argv[i][1] == 's' && argv[i][2] == 'x') {
				// sheetsource, sheet is normal to x direction
				sscanf(argv[i]+3, "%lf", &data->port.x0);
				fprintf(stderr, "Sheet source for x = %lf\n", data->port.x0);
				data->port.y0 = 1e50;
				data->port.solver = SHEETSOURCE;
				data->realflag = 1;
				data->par.num_mode = 1;
			}

			if (argv[i][1] == 's' && argv[i][2] == 'y') {
				// sheetsource, sheet is normal to y direction
				sscanf(argv[i]+3, "%lf", &data->port.x0);
				fprintf(stderr, "Sheet source for y = %lf\n", data->port.y0);
				data->port.x0 = 1e50;
				data->port.solver = SHEETSOURCE;
				data->realflag = 1;
				data->par.num_mode = 1;
			}

			if (argv[i][1] == 's' && argv[i][2] == 'd') {
				fprintf(stderr, "Single direction input\n");
				data->port.inputDirection = SINGLE_DIRECTION;
				data->port.forward = 1.0;
				data->port.backward = 0.0;
			}
		}
	}
}

extern int checkBandWidth(FEM *fem, DataTable *data)
{
	int i, j, k;
	int min, max, bandwidth;
	
	bandwidth = 0;
	for (k = 0; k < fem->ne; k++) {
		min = fem->np;
		max = 0;
		for (i = 0; i < 6; i++) {
			if (fem->element[k].kk[i] < fem->nr) {
				if (fem->element[k].kk[i] < min) min = fem->element[k].kk[i];
				if (fem->element[k].kk[i] > max) max = fem->element[k].kk[i];
			}
		}
		if (bandwidth < max-min) bandwidth = max-min;
	}	
	return bandwidth;
}

extern int inputData(DataTable *data, char *file)
{
	int		i, j,k, ii, jj, ll, mm;
	char 	string[256];
	FILE	*fp, *fpM;
	Param	*par = &(data->par);
	FEM *fem = &(data->fem);
	Port *port = &(data->port);
	int count, n_mode, numoflambda;
	double re, im;
	double re_tmp, im_tmp;
	double_complex cj(0.0, 1.0);
	double xr, yr, xin, yin, angle;
	double temp1, temp2;
	
	strcpy(data->file, file);
	sprintf(string, "%s.msh", data->file);
	fp = fopen(string, "r");
	
	fscanf(fp, "%lf", &(data->par.wavelength));
	fscanf(fp, "%d", &(data->par.modeID));
	fscanf(fp, "%d", &(data->par.modeNo));
	fscanf(fp, "%d", &data->par.n_material);
	for(i = 0; i < data->par.n_material; i++) {
		fscanf(fp, "%lf %lf", &re_tmp, &im_tmp);
	}
	
	fscanf(fp, "%d", &(data->par.symFlag));
	fscanf(fp, "%lf %lf", &(data->par.symAxisX), &(data->par.symAxisY));
	
	fscanf(fp, "%d", &(data->fem.np));
	fscanf(fp, "%d", &(data->fem.nr));

	data->fem.nr = data->fem.np;

	ALLOCATION(fem->node, XYdouble, fem->np);
	fem->xmax = fem->ymax = -1.0e50;
	fem->xmin = fem->ymin =	1.0e50;
	for (i = 0; i < fem->np; i++) {
		fscanf(fp, "%lf %lf", &(fem->node[i].x), &(fem->node[i].y));
		if (i == 0 || fem->node[i].x > fem->xmax) fem->xmax = fem->node[i].x;
		if (i == 0 || fem->node[i].y > fem->ymax) fem->ymax = fem->node[i].y;
		if (i == 0 || fem->node[i].x < fem->xmin) fem->xmin = fem->node[i].x;
		if (i == 0 || fem->node[i].y < fem->ymin) fem->ymin = fem->node[i].y;
	}
	
	fscanf(fp, "%*d");
	fscanf(fp, "%d", &(data->fem.ne));
	ALLOCATION(fem->element, Element3D, fem->ne);
	for (i = 0; i < fem->ne; i++) {
		fscanf(fp, "%d", &(fem->element[i].matID));
		fem->element[i].matID--;
		for (j = 0; j < 6; j++) {
			fscanf(fp, "%d", &(fem->element[i].kk[j]));
			fem->element[i].kk[j]--;
		}
	}

	fclose(fp);
	
	fprintf(stderr, "x = (%lf -- %lf), y = (%lf -- %lf)\n", fem->xmin, fem->xmax, fem->ymin, fem->ymax);
	
	ALLOCATION(fem->inFlag, int, fem->np);
	ALLOCATION(fem->FinR, double_complex, fem->np);
	ALLOCATION(fem->FinL, double_complex, fem->np);
	
	//	angle = port->data[data->par.inputPort].angle*PI/180.0;
	angle = 0.0;

	fprintf(stderr, "angle = %lf\n", angle);
	
	// read .pre file
	//------------------------------------------------------------------
	sprintf(string, "%s.pre", data->file);
	if ( (fp = fopen(string, "r")) == NULL ) {
		fprintf(stderr, "can't open input file ( %s ).\n", string);
		exit(1);
	}
	
	fscanf(fp, "%*s %*s");
	// 1st initial wavelength 2nd final wavelength 3rd wavelength step
	fscanf(fp, "%lf %lf %lf", &(data->par.w_start), &(data->par.w_end), &(data->par.w_step));

	data->par.numoflambda = (int)((data->par.w_end-data->par.w_start)/data->par.w_step)+1;

	data->par.wavelength = data->par.w_start;
	data->par.k0 = 2.0*PI/data->par.wavelength;
	data->par.k02 = data->par.k0*data->par.k0;
	fprintf(stderr, "start,end,step, numof	wavelength = %lf, %lf, %lf, %d\n", 
		data->par.w_start, data->par.w_end, data->par.w_step,
		data->par.numoflambda);

	fscanf(fp, "%*s %d", &(data->par.modeID));

	// input power
	fscanf(fp, "%*s = %lf", &data->par.power);

	// the number of material
	fscanf(fp, "%*s %*s %*s %*s = %d", &data->par.n_material);
 
	// modifiedEIM flag
	fscanf(fp, "%*s %*s = %d", &data->par.modEIMflag);
	
	fscanf(fp, "%*s %*s %*s %*s");
	// 1st real part of refractive index
	// 2nd real part of magnetic permitivity
	// 3rd nonlinear index
	// 4th saturable nonlinear index 

	for(i = 0; i < data->par.n_material; i++) {
		fscanf(fp, "%lf %lf %lf %lf", &(par->er[i]), &(par->mr[i]),	
			&(data->nonlinear.kerr[i]), &(data->nonlinear.sat[i]));

		par->er[i] *= par->er[i];
		par->mr[i] *= par->mr[i];	
	
		par->muT[i].xx = par->muT[i].yy = par->muT[i].zz = par->mr[i];
		par->muT[i].xy = par->muT[i].xz = 0.0;
		par->muT[i].yx = par->muT[i].yz = 0.0;
		par->muT[i].zx = par->muT[i].zy = 0.0;
		par->epT[i].xx = par->epT[i].yy = par->epT[i].zz = par->er[i];
		par->epT[i].xy = par->epT[i].xz = 0.0;
		par->epT[i].yx = par->epT[i].yz = 0.0;
		par->epT[i].zx = par->epT[i].zy = 0.0;

		fprintf(stderr, "%d th index = %lf\n", i, data->par.er[i]);
	}
	
	// field related parameters
	fscanf(fp, "%*s %*s %*s");
	fscanf(fp, "%*s = %d", &(data->par.div_x));
	fscanf(fp, "%*s = %d", &(data->par.div_y));

	// X cross section
	fscanf(fp, "%*s %*s %*s %d", &(data->par.numx));
	fscanf(fp, "%*s = %d", &(data->par.division_x));

	fprintf(stderr, "%d %d\n", data->par.numx, data->par.division_x);

	ALLOCATION(data->par.center_y, double, data->par.numx);
	fscanf(fp, "%*s");
	for(i = 0; i < data->par.numx; i++) {
		fscanf(fp, "%lf %d", &(data->par.center_y[i]), 
			&(data->par.inflag[i])); 
	}

	// Y cross section
	fscanf(fp, "%*s %*s %*s %d", &(data->par.numy));
	fscanf(fp, "%*s = %d", &(data->par.division_y));

	fprintf(stderr, "%d %d\n", data->par.numy, data->par.division_y);

	ALLOCATION(data->par.center_x, double, data->par.numy);
	fscanf(fp, "%*s");
	for(i = 0; i < data->par.numy; i++) {
		fscanf(fp, "%lf %d", &(data->par.center_x[i]), 
			&(data->par.inflag[data->par.numx+i])); 

		// check the input port number
		if (data->par.inflag[data->par.numx+i] == 1) {
			data->par.inputPort = data->par.numx+i;
		}
	}

	// number of input modes
	fscanf(fp, "%*s %*s %*s %*s = %d", &(data->par.num_mode));

	fprintf(stderr, "Input port No = %d\n", data->par.inputPort);
	fprintf(stderr, "Number of Light Input = %d\n", data->par.num_mode);

	// PML information
	// order, tandelta
	fscanf(fp, "%*s = %d, %lf", &(data->pml.m), &(data->pml.tanD));

	fscanf(fp, "%*s = %d, %d", &(data->pml.nx), &(data->pml.ny));
	for (i = 0; i < data->pml.nx; i++) {
		fscanf(fp, "%lf %lf %lf %lf",
			&(data->pml.x_data[i].x0),
			&(data->pml.x_data[i].y0), &(data->pml.x_data[i].y1),
			&(data->pml.x_data[i].dx));
		data->pml.x_data[i].x1 = data->pml.x_data[i].x0+data->pml.x_data[i].dx;
	}
	for (i = 0; i < data->pml.ny; i++) {
		fscanf(fp, "%lf %lf %lf %lf",
			&(data->pml.y_data[i].y0),
			&(data->pml.y_data[i].x0), &(data->pml.y_data[i].x1),
			&(data->pml.y_data[i].dy));
		data->pml.y_data[i].y1 = data->pml.y_data[i].y0+data->pml.y_data[i].dy;
	}
	//------------------------------------------------------------------
	// port info check
	port->number = data->par.numx + data->par.numy;
	ii = 0;
	n_mode = data->par.num_mode;
	numoflambda = data->par.numoflambda;

	ALLOCATION(port->data, PortData*, data->port.number);
	ALLOCATION(port->dataMemory, PortData**, numoflambda);

	for (ll = 0; ll < numoflambda; ll++) {
		ALLOCATION(port->dataMemory[ll], PortData*, n_mode);
		for (mm = 0; mm < n_mode; mm++) {
			ALLOCATION(port->dataMemory[ll][mm], PortData, port->number);
		}
	}
	fprintf(stderr, "Alloc %d, %d, %d\n", numoflambda, n_mode, port->number);


	for (i = 0; i < data->par.numx; i++) {
		for (ll = 0; ll < numoflambda; ll++) {
			for (mm = 0; mm < n_mode; mm++) {
				port->dataMemory[ll][mm][i].normal.x = 0.0;
				port->dataMemory[ll][mm][i].normal.y = 1.0;

				// the number of point on the port
				port->dataMemory[ll][mm][i].np = 0;
				for (j = 0; j < fem->np; j++) {
					if (fabs(fem->node[j].y-data->par.center_y[i]) < 1e-9) {
					port->dataMemory[ll][mm][i].np++;
					}
				}

				port->dataMemory[ll][mm][i].ne = (port->dataMemory[ll][mm][i].np-1)/2;


				// the point number of each point on the port
				ALLOCATION(port->dataMemory[ll][mm][i].kp, int, port->dataMemory[ll][mm][i].np);
				ALLOCATION(port->dataMemory[ll][mm][i].xp, double, port->dataMemory[ll][mm][i].np);
				
				jj = 0;
				for (j = 0; j < fem->np; j++) {
					if (fabs(fem->node[j].y-data->par.center_y[i]) < 1e-9) {
					port->dataMemory[ll][mm][i].xp[jj] = fem->node[j].x;
					port->dataMemory[ll][mm][i].kp[jj] = j;
					
					/* printf("# xp, yp, kp = %lf %lf %d\n",
							port->dataMemory[ll][mm][i].xp[jj],
							data->fem.node[port->dataMemory[ll][mm][i].kp[jj]].y, 
							port->dataMemory[ll][mm][i].kp[jj]); */
							
					jj++;
					}
				}

				// sorting x coordinate		 
				for (j = 0; j < port->dataMemory[ll][mm][i].np; j++) {
					temp1 = port->dataMemory[ll][mm][i].xp[j];
					temp2 = port->dataMemory[ll][mm][i].kp[j];
					k = j-1;
					while(k > -1 && port->dataMemory[ll][mm][i].xp[k] > temp1) {
					port->dataMemory[ll][mm][i].xp[k+1] = port->dataMemory[ll][mm][i].xp[k];
					port->dataMemory[ll][mm][i].kp[k+1] = port->dataMemory[ll][mm][i].kp[k];
					k--;
					}
					port->dataMemory[ll][mm][i].xp[k+1] = temp1;
					port->dataMemory[ll][mm][i].kp[k+1] = temp2;
				}
			}
		}
	}

	for(i = 0; i < data->par.numy; i++) {
		for (ll = 0; ll < numoflambda; ll++) {
			for (mm = 0; mm < n_mode; mm++) {
				port->dataMemory[ll][mm][data->par.numx+i].normal.x = 1.0;
				port->dataMemory[ll][mm][data->par.numx+i].normal.y = 0.0;

				// the number of point on the port
				port->dataMemory[ll][mm][data->par.numx+i].np = 0;
				for (j = 0; j < fem->np; j++) {
					if (fabs(fem->node[j].x-data->par.center_x[i]) < 1e-9) {
					port->dataMemory[ll][mm][data->par.numx+i].np += 1;
					}
				}

				port->dataMemory[ll][mm][data->par.numx+i].ne = 
					(port->dataMemory[ll][mm][data->par.numx+i].np-1)/2;


				// the point number of each point on the port
				ALLOCATION(port->dataMemory[ll][mm][data->par.numx+i].kp,
							int, port->dataMemory[ll][mm][data->par.numx+i].np);
				ALLOCATION(port->dataMemory[ll][mm][data->par.numx+i].yp,
							double, port->dataMemory[ll][mm][data->par.numx+i].np);
				
				jj = 0;
				for (j = 0; j < fem->np; j++) {
					if (fabs(fem->node[j].x-data->par.center_x[i]) < 1e-9) {
					port->dataMemory[ll][mm][data->par.numx+i].yp[jj] = fem->node[j].y;
					port->dataMemory[ll][mm][data->par.numx+i].kp[jj] = j;
					
					/* printf("# xp, yp, kp = %lf %lf %d\n", 
							data->fem.node[port->dataMemory[ll][mm][data->par.numx+i].kp[jj]].x,
							port->dataMemory[ll][mm][data->par.numx+i].yp[jj], 
							port->dataMemory[ll][mm][data->par.numx+i].kp[jj]); */
					jj++;
					}
				}
				
				// sorting y coordinate		 
				for (j = 0; j < port->dataMemory[ll][mm][data->par.numx+i].np; j++) {
					temp1 = port->dataMemory[ll][mm][data->par.numx+i].yp[j];
					temp2 = port->dataMemory[ll][mm][data->par.numx+i].kp[j];
					k = j-1;
					while(k > -1 && port->dataMemory[ll][mm][data->par.numx+i].yp[k] > temp1) {
					port->dataMemory[ll][mm][data->par.numx+i].yp[k+1] = port->dataMemory[ll][mm][data->par.numx+i].yp[k];
					port->dataMemory[ll][mm][data->par.numx+i].kp[k+1] = port->dataMemory[ll][mm][data->par.numx+i].kp[k];
					k--;
					}
					port->dataMemory[ll][mm][data->par.numx+i].yp[k+1] = temp1;
					port->dataMemory[ll][mm][data->par.numx+i].kp[k+1] = temp2;
				}
			}
		} 
	}

	fprintf(stderr, "The number of ports\n\
			X = %d, Y = %d, Total = %d\n",
			data->par.numx, data->par.numy, port->number);

	for(i = 0; i < port->number; i++) {
		for (ll = 0; ll < numoflambda; ll++) {
			for (mm = 0; mm < n_mode; mm++) {
				fprintf(stderr, "The number of point and ne on port%d (orthogonal to %s axis, analyzing %d th wavelength and TE%d mode) = %d %d\n",
				i, (port->dataMemory[ll][mm][i].normal.x==1.0 ? "X" : "Y"), ll, mm, port->dataMemory[ll][mm][i].np, port->dataMemory[ll][mm][i].ne);
					
				/*
				for (j = 0; j < port->dataMemory[ll][mm][i].np; j++) {
					printf("xp, yp, kp = %lf %lf %d\n", 
					dataMemory->fem.node[port->dataMemory[ll][mm][i].kp[j]].x,
					port->dataMemory[ll][mm][i].yp[j], 
					port->dataMemory[ll][mm][i].kp[j]);
				}

				printf("\n");
				*/
			}
		}
	}
	
	//------------------------------------------------------------------
	data->fem.nbw = checkBandWidth(fem, data);
	fprintf(stderr, "half band width = %d\n", data->fem.nbw);
}

extern int rotationCoordinate(double x0, double y0,
					double *xr0, double *yr0, double angle)
{
	double cosA, sinA;
	
	cosA = cos(angle); sinA = sin(angle);
	*xr0 = cosA*x0-sinA*y0;
	*yr0 = sinA*x0+cosA*y0;
	// fprintf(stderr, "(%lf, %lf) --> (%lf, %lf)\n", x0, y0, *xr0, *yr0);
}
