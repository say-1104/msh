#include "optfem.h"

#define NX		50
#define NY		50


double calcFitness(DataTable *data, int ll) {
	int mm;
	int n_mode = data->par.num_mode;
	double gamma; 
	double_complex cj(0.0, 1.0);

	double temp = 0.0;
	if (data->mosaic.Device == POWER_SPLITTER) {
		if (n_mode != 1) {
			fprintf(stderr, "マルチモードには対応していません。");
			exit(1);
		}
		gamma = data->mosaic.gamma;
		if (data->mosaic.evaluation == 1) {
			// ユークリッド距離用
			double norm = sqrt(pow(gamma, 2) + pow(1-gamma, 2));
			temp +=   1 - sqrt(pow(data->mosaic.TT[1][ll]-gamma, 2) + pow(data->mosaic.TT[2][ll]-(1.0-gamma), 2))/norm;
		}
		else if (data->mosaic.evaluation == 2) {
			// 市街地距離用
			temp +=   1 - (fabs(data->mosaic.TT[1][ll]-gamma) + fabs(data->mosaic.TT[2][ll]-(1.0-gamma)));
		}
		else if (data->mosaic.evaluation == 3) {
			// チェス版距離用
			temp +=   1 - 2*MAX(fabs(data->mosaic.TT[1][ll]-gamma), fabs(data->mosaic.TT[2][ll]-(1.0-gamma)));
		}
		else if (data->mosaic.evaluation == 4) {
			temp +=   1 - (fabs(data->mosaic.TT[1][ll]-gamma) + fabs(data->mosaic.TT[2][ll]-(1.0-gamma)))
						- fabs(data->mosaic.TT[1][ll] - data->mosaic.TT[2][ll]);
		}
		else if (data->mosaic.evaluation == 5) {
			double p_diff = 180*atan2(imag(data->mosaic.amp[1][ll]/data->mosaic.amp[2][ll]), real(data->mosaic.amp[1][ll]/data->mosaic.amp[2][ll]))/M_PI;
			double arg_x = data->mosaic.TT[1][ll]+data->mosaic.TT[2][ll]-0.8;
			p_diff = fabs(fabs(p_diff) - 90.0);
			temp +=   1 - (fabs(data->mosaic.TT[1][ll]-gamma) + fabs(data->mosaic.TT[2][ll]-(1.0-gamma)))
						- fabs(data->mosaic.TT[1][ll] - data->mosaic.TT[2][ll]) - 1/(1+exp(-20*arg_x))*p_diff/100.0;
		}
		else if (data->mosaic.evaluation == 6) {
			double p_diff = 180*atan2(imag(data->mosaic.amp[1][ll]/data->mosaic.amp[2][ll]), real(data->mosaic.amp[1][ll]/data->mosaic.amp[2][ll]))/M_PI;
			temp +=   1 - (fabs(data->mosaic.TT[1][ll]-gamma) + fabs(data->mosaic.TT[2][ll]-(1.0-gamma)))
						- fabs(data->mosaic.TT[1][ll] - data->mosaic.TT[2][ll]) - 0.01*fabs(fabs(p_diff) - 90.0);
		}
		else {
			double_complex d1 = data->mosaic.amp[1][ll] - cj*data->mosaic.amp[2][ll];
			double_complex d2 = cj*data->mosaic.amp[1][ll] - data->mosaic.amp[2][ll];
			double tot = data->mosaic.TT[1][ll] + data->mosaic.TT[2][ll];
			temp += tot - (abs(d1) < abs(d2) ? abs(d1) : abs(d2));

			/*
			double_complex d = data->mosaic.amp[1][ll] - cj*data->mosaic.amp[2][ll];
			double tot = data->mosaic.TT[1][ll] + data->mosaic.TT[2][ll];
			temp += tot - sqrt(real(d)*real(d)+imag(d)*imag(d));
			*/
		}
	}
	else if (data->mosaic.Device == WAVEGUIDE_CROSSING) {
		for (mm = 0; mm < n_mode; mm++) {
			if (data->par.modEIMflag == 1) temp += data->mosaic.TT[mm][ll*n_mode+mm];
			else temp += data->mosaic.TT[mm][0];
		}
		temp /= n_mode;
	}
	else if (data->mosaic.Device == MODE_MUX) {
		for (mm = 0; mm < n_mode; mm++) {
			if (data->par.modEIMflag == 1) temp += data->mosaic.TT[mm][ll*n_mode+mm];
			else temp += data->mosaic.TT[mm][0];
		}
		temp /= n_mode;
	}
	else {
		temp += (1.0-data->mosaic.TT[1][ll]);
	}
	return temp;
}

// Direct Binary Search for mosaic pattern
extern int DirectBinarySearch(DataTable *data)
{
	int i, j, k, l;
	int ii, jj, ll;
	int state, change, sum_change, Ny, Nz, tmpi, tmpj;
	Param	*par = &(data->par);
	Port	*port = &(data->port);
	double temp;
	char filename[256];
	FILE *fp;
	int totalNN;
	int *sequence;
	double* fitness;
	double_complex cj(0.0, 1.0);
	int ysym = (data->mosaic.symflag & Y_SYMMETRY), zsym = (data->mosaic.symflag & Z_SYMMETRY), yeqzsym = (data->mosaic.symflag & YeqZ_SYMMETRY);

	fprintf(stderr, "Execute Direct Binary Search\n");

	// Prepare output files
	sprintf(filename, "Output/fitness.data");
	if ((fp = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "can't open file ( %s )\n", filename);
		exit(1);
	}
	fclose(fp);

	sprintf(filename, "Output/pattern.data");
	if ((fp = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "can't open file ( %s )\n", filename);
		exit(1);
	}
	fclose(fp);
	
	sprintf(filename, "Output/sequence.data");
	if ((fp = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "can't open file ( %s )\n", filename);
		exit(1);
	}
	fclose(fp);
	
	if (data->mosaic.Device == POWER_SPLITTER) {
		sprintf(filename, "Output/extended.data");
		if ((fp = fopen(filename, "w")) == NULL) {
			fprintf(stderr, "can't open file ( %s )\n", filename);
			exit(1);
		}
		fclose(fp);
	}

	
	// Check the pixels which need to cosider
	if (ysym) {
		Ny = (data->mosaic.Ny+1)/2;
	}
	else Ny = data->mosaic.Ny;

	if (zsym) {
		Nz = (data->mosaic.Nz+1)/2;
	}
	else Nz = data->mosaic.Nz;
	
	if (yeqzsym) {
		if (Ny != Nz) {
			fprintf(stderr, "this structure can't implement Y=Z symmetry, because Ny isn't equal to Nz");
			exit(1);
		}

		totalNN = Ny*(Ny+1)/2;
	}
	else {
		totalNN = Ny*Nz;
	}

	if (data->mosaic.PixelReduction != 0) {
		totalNN = 0;
		for (i=0; i<Ny; i++) {
			for (j=0; j<Nz; j++) {
				if (!yeqzsym || j <= i) {
					if (data->mosaic.BorW[i][j] != -1) {
						totalNN++;
					}
				}
			}
		}
	}

	// Calculate the first FOM
	fprintf(stderr, "Performance of initial pattern\n");
	ALLOCATION(sequence, int, totalNN);
	ALLOCATION(fitness, double, data->par.numoflambda);
	ALLOCATION(data->mosaic.DBSseq, int, totalNN);
	ALLOCATION(data->mosaic.FOM_history, double, data->mosaic.DBSitr+1);

	k=0;
	for(i = 0; i < Ny*Nz; i++) {
		jj = i%Nz;
		ii = i/Nz;
		if (data->mosaic.PixelReduction != 0 && data->mosaic.BorW[ii][jj] == -1) continue;
		if (!yeqzsym || jj <= ii) {
			sequence[k] = i;
			k++;
		}
	}

	data->mosaic.FOM = 0.0;
	calcFEM(data);
	//outputOverlap(data, k);
	for(i = 0; i < data->par.numoflambda; i++) {
		fitness[i] = calcFitness(data, i);
		data->mosaic.FOM += fitness[i];
	}
	data->mosaic.FOM /= data->par.numoflambda;
	data->mosaic.FOM_history[0] = data->mosaic.FOM;

	
	sprintf(filename, "Output/fitness.data");
	if ((fp = fopen(filename, "a+")) == NULL) {
		fprintf(stderr, "can't open file ( %s )\n", filename);
		exit(1);
	}
	
	fprintf(fp, "0: ");
	fprintf(stderr, "0: ");
	fprintf(fp, "%lf, %d pixels, eval:%d", data->mosaic.FOM, totalNN, data->mosaic.evaluation);
	fprintf(stderr, "%lf, %d pixels, eval:%d", data->mosaic.FOM, totalNN, data->mosaic.evaluation);
	if (data->par.numoflambda != 1) {
		fprintf(fp, " / ");
		fprintf(stderr, " / ");
		
		for(i = 0; i < data->par.numoflambda; i++) {
			fprintf(fp, "%lf,", fitness[i]);
			fprintf(stderr, "%lf,", fitness[i]);
		}
	}
	fprintf(fp, "\n");
	fprintf(stderr, "\n");
	fclose(fp);

	
	if (data->mosaic.Device == POWER_SPLITTER) {
		sprintf(filename, "Output/extended.data");
		if ((fp = fopen(filename, "a+")) == NULL) {
			fprintf(stderr, "can't open file ( %s )\n", filename);
			exit(1);
		}
		fprintf(fp, "0: ");
		i=0;
		fprintf(fp, "%lf,%lf,%lf,%lf", real(data->mosaic.amp[1][i]), imag(data->mosaic.amp[1][i]), real(data->mosaic.amp[2][i]), imag(data->mosaic.amp[2][i]));
		for(i = 1; i < data->par.numoflambda; i++) {
			fprintf(fp, " / %lf,%lf,%lf,%lf", real(data->mosaic.amp[1][i]), imag(data->mosaic.amp[1][i]), real(data->mosaic.amp[2][i]), imag(data->mosaic.amp[2][i]));
		}
		fprintf(fp, "\n");
		fclose(fp);
	}

	change = 0;
	for (k=0; k<data->mosaic.DBSitr; k++) {
		// Unable to resume with prepared sequence.
		shuffle(sequence, totalNN);
		sprintf(filename, "Output/sequence.data");

		if ((fp = fopen(filename, "a+")) == NULL) {
			fprintf(stderr, "can't open file ( %s )\n", filename);
			exit(1);
		}
		fprintf(fp, "%d-th iteration\n", k);
		for (i=0; i<totalNN; i++) {
			jj = sequence[i]%Nz;
			ii = sequence[i]/Nz;
			fprintf(fp, "%d jj: %d ii: %d\n", sequence[i], jj, ii);
		}
		fclose(fp);

		for (l=0; l<totalNN; l++) {

			// toggle one pixel state
			jj = sequence[l]%Nz;
			ii = sequence[l]/Nz;
			tmpi = (ysym ? data->mosaic.Ny-1-ii : ii);
			tmpj = (zsym ? data->mosaic.Nz-1-jj : jj);
			state = data->mosaic.BorW[ii][jj];
			data->mosaic.BorW[ii][jj] = 1-state;
			data->mosaic.BorW[ii][tmpj] = 1-state;
			data->mosaic.BorW[tmpi][jj] = 1-state;
			data->mosaic.BorW[tmpi][tmpj] = 1-state;
			if (yeqzsym) {
				data->mosaic.BorW[jj][ii] = 1-state;
				data->mosaic.BorW[tmpj][ii] = 1-state;
				data->mosaic.BorW[jj][tmpi] = 1-state;
				data->mosaic.BorW[tmpj][tmpi] = 1-state;
			}

			temp = 0.0;
			// initialize fitness
			calcFEM(data);
			//outputOverlap(data, k);
			for (i = 0; i < data->par.numoflambda; i++) {
				fitness[i] = calcFitness(data, i);
				temp += fitness[i];
			}
			temp /= data->par.numoflambda;

			change = (temp>data->mosaic.FOM);

			sprintf(filename, "Output/fitness.data");
			if ((fp = fopen(filename, "a+")) == NULL) {
				fprintf(stderr, "can't open file ( %s )\n", filename);
				exit(1);
			}
			
			fprintf(fp, "%d: ", k*totalNN+l+1);
			fprintf(stderr, "%d: ", k*totalNN+l+1);
			fprintf(fp, "%lf,%lf,%d", temp, data->mosaic.FOM, change);
			fprintf(stderr, "%lf,%lf,%d", temp, data->mosaic.FOM, change);
			if (data->par.numoflambda != 1) {
				fprintf(fp, " / ");
				fprintf(stderr, " / ");
				for(i = 0; i < data->par.numoflambda; i++) {
					fprintf(fp, "%lf,", fitness[i]);
					fprintf(stderr, "%lf,", fitness[i]);
				}
			}
			fprintf(fp, "\n");
			fprintf(stderr, "\n");
			fclose(fp);
			
			if (data->mosaic.Device == POWER_SPLITTER) {
				sprintf(filename, "Output/extended.data");
				if ((fp = fopen(filename, "a+")) == NULL) {
					fprintf(stderr, "can't open file ( %s )\n", filename);
					exit(1);
				}
				fprintf(fp, "%d-%dth pixel: ", k, l);
				i=0;
				fprintf(fp, "%lf,%lf,%lf,%lf", real(data->mosaic.amp[1][i]), imag(data->mosaic.amp[1][i]), real(data->mosaic.amp[2][i]), imag(data->mosaic.amp[2][i]));
				for(i = 1; i < data->par.numoflambda; i++) {
					fprintf(fp, " / %lf,%lf,%lf,%lf", real(data->mosaic.amp[1][i]), imag(data->mosaic.amp[1][i]), real(data->mosaic.amp[2][i]), imag(data->mosaic.amp[2][i]));
				}
				for(i = 0; i < data->par.numoflambda; i++) {
					fprintf(fp, "::%lf,%lf", real(data->mosaic.amp[1][i] - cj*data->mosaic.amp[2][i]), imag(data->mosaic.amp[1][i] - cj*data->mosaic.amp[2][i]));
				}
				fprintf(fp, "\n");
				fclose(fp);
			}
			
			if (change) {
				data->mosaic.FOM = temp;
				
				sprintf(filename, "Output/pattern.data");
				if ((fp = fopen(filename, "a+")) == NULL) {
					fprintf(stderr, "can't open file ( %s )\n", filename);
					exit(1);
				}
					
				fprintf(fp, "%d, %lf:", k*totalNN+l+1, data->mosaic.FOM);
				fprintf(stderr, "%d, %lf:", k*totalNN+l+1, data->mosaic.FOM);
				for(i = 0; i < data->mosaic.Ny; i++) {
					for(j = 0; j < data->mosaic.Nz; j++) {
						if (i == data->mosaic.Ny-1 && j == data->mosaic.Nz-1) {
							fprintf(fp, "%d", data->mosaic.BorW[i][j]);
						}
						else {
							fprintf(fp, "%d,", data->mosaic.BorW[i][j]);
						}
					}
				}
				fprintf(fp, "\n");
				fclose(fp);
			}
			else {
				// go back to the original state
				jj = sequence[l]%Nz;
				ii = sequence[l]/Nz;
				tmpi = (ysym ? data->mosaic.Ny-1-ii : ii);
				tmpj = (zsym ? data->mosaic.Nz-1-jj : jj);
				state = data->mosaic.BorW[ii][jj];
				data->mosaic.BorW[ii][jj] = 1-state;
				data->mosaic.BorW[ii][tmpj] = 1-state;
				data->mosaic.BorW[tmpi][jj] = 1-state;
				data->mosaic.BorW[tmpi][tmpj] = 1-state;
				if (yeqzsym) {
					data->mosaic.BorW[jj][ii] = 1-state;
					data->mosaic.BorW[tmpj][ii] = 1-state;
					data->mosaic.BorW[jj][tmpi] = 1-state;
					data->mosaic.BorW[tmpj][tmpi] = 1-state;
				}
			}
		}

		sprintf(filename, "Output/GuessedMatrix-%d.csv", k);
		if ((fp = fopen(filename, "w")) == NULL) {
			fprintf(stderr, "can't open file ( %s )\n", filename);
			exit(1);
		}
			
		for(i = 0; i < data->mosaic.Ny; i++) {
			for(j = 0; j < data->mosaic.Nz; j++) {
				if (i == data->mosaic.Ny-1 && j == data->mosaic.Nz-1) {
					fprintf(fp, "%d", data->mosaic.BorW[i][j]);
				}
				else {
					fprintf(fp, "%d,", data->mosaic.BorW[i][j]);
				}
			}
			fprintf(fp, "\n");
		}
		fclose(fp);

		data->mosaic.FOM_history[k+1] = data->mosaic.FOM;
		if (data->mosaic.FOM_history[k+1]-data->mosaic.FOM_history[k] < data->mosaic.termination + 1e-6) {
			fprintf(stderr, "FOM improvement is %lf, which is no greater than %lf", data->mosaic.FOM_history[k+1]-data->mosaic.FOM_history[k], data->mosaic.termination);
			break;
		}
	}

	sprintf(filename, "Output/fitness.data");
	if ((fp = fopen(filename, "a+")) == NULL) {
		fprintf(stderr, "can't open file ( %s )\n", filename);
		exit(1);
	}

	for(i=0; i<=k+1; i++) {
		fprintf(fp, "%d-th itr: %lf, ", i, data->mosaic.FOM_history[i]);
	}
	fprintf(fp, "\n");
	fclose(fp);

	sprintf(filename, "Output/GuessedMatrix.csv");
	if ((fp = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "can't open file ( %s )\n", filename);
		exit(1);
	}
		
	for(i = 0; i < data->mosaic.Ny; i++) {
		for(j = 0; j < data->mosaic.Nz; j++) {
			if (i == data->mosaic.Ny-1 && j == data->mosaic.Nz-1) {
				fprintf(fp, "%d", data->mosaic.BorW[i][j]);
			}
			else {
				fprintf(fp, "%d,", data->mosaic.BorW[i][j]);
			}
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

extern int calcFEM(DataTable *data)
{
	int ll, mm;
	Param	*par = &(data->par);
	Port	*port = &(data->port);
	char filename[256];

	for (ll = 0; ll < par->numoflambda; ll++) {
		par->wavelength = par->w_start+(double)ll*par->w_step;
		fprintf(stderr, "wavelength = %lf\n", par->wavelength);
			
		for (mm = 0; mm < par->num_mode; mm++) {
			SetModifiedEI(data, par->wavelength, ll, mm);
			if (port->solver == INPUTEIGEN) {
				setPort(data, par->wavelength, ll, mm);
			}
			
			FEManalysis(data, par->wavelength, ll);
			
			if (port->solver != POINTSOURCE) {
				portField(data, ll);
				probability(data, par->wavelength, ll, mm);
			}
			
			//intpol(data, par->wavelength, mm);
		}
	}
}

