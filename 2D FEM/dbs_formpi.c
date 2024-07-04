#include "optfem.h"

#define NX		50
#define NY		50


double calcFitness(DataTable *data){
	int ll, mm;
	int n_mode = data->par.num_mode;
	double gamma; 

	double temp = 0.0;
	if (data->mosaic.Device == POWER_SPLITTER) {
		gamma = data->mosaic.gamma;
		for (ll = 0; ll < data->par.numoflambda; ll++) {
			temp += 1 	- (fabs(data->mosaic.TT[1][ll]-gamma) + fabs(data->mosaic.TT[2][ll]-(1.0-gamma)))
						- fabs(data->mosaic.TT[1][ll] - data->mosaic.TT[2][ll]);
		}
		temp /= (double)data->par.numoflambda;
	}
	else if (data->mosaic.Device == WAVEGUIDE_CROSSING) {
		for (ll = 0; ll < data->par.numoflambda; ll++) {
			for (mm = 0; mm < n_mode; mm++) {
				temp += data->mosaic.TT[mm][ll*n_mode+mm];
			}
		}
		temp /= n_mode;
		temp /= (double)data->par.numoflambda;  
	}
	else {
		for (ll = 0; ll < data->par.numoflambda; ll++) {
			temp += (1.0-data->mosaic.TT[1][ll]);
		}
		temp /= (double)data->par.numoflambda;
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
	int ysym = (data->mosaic.symflag & Y_SYMMETRY), zsym = (data->mosaic.symflag & Z_SYMMETRY), yeqzsym = (data->mosaic.symflag & YeqZ_SYMMETRY);

	fprintf(stderr, "Direct Binary Search\n");

	if (data->mpi.my_rank == 0) {
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
	}

	if (ysym != 0) {
		Ny = (data->mosaic.Ny+1)/2;
	}
	else Ny = data->mosaic.Ny;

	if (zsym != 0) {
		Nz = (data->mosaic.Nz+1)/2;
	}
	else Nz = data->mosaic.Nz;
	
	if (yeqzsym != 0) {
		if (Ny != Nz) {
			fprintf(stderr, "this structure can't implement Y=Z symmetry, because Ny isn't equal to Nz");
			exit(1);
		}

		totalNN = Ny*(Ny+1)/2;
		if (data->mosaic.PixelReduction != 0) {
			totalNN = 0;
			for (i=0; i<Ny; i++) {
				for (j=0; j<Nz; j++) {
					if (j <= i) {
						if (data->mosaic.BorW[i][j] != -1) {
							totalNN++;
						}
					}
				}
			}
		}
	}
	else {
		totalNN = Ny*Nz;
	}

	// Calculate the first FOM
	fprintf(stderr, "Performance of initial pattern\n");
	ALLOCATION(fitness, double, data->mpi.num_proc);
	ALLOCATION(sequence, int, totalNN);
	ALLOCATION(data->mosaic.fitness, double, data->mpi.num_proc);
	ALLOCATION(data->mosaic.DBSseq, int, totalNN);
	ALLOCATION(data->mosaic.FOM_history, double, data->mosaic.DBSitr+1);

	k=0;
	if (data->mpi.my_rank == 0) {
		for(i = 0; i < Ny*Nz; i++) {
			jj = i%Nz;
			ii = i/Nz;
			if (data->mosaic.PixelReduction != 0 && data->mosaic.BorW[ii][jj] == -1) continue;
			if (yeqzsym == 0 || jj <= ii) {
				sequence[k] = i;
				k++;
			}
		}
	}
	else {
		for(i = 0; i < totalNN; i++) {
			sequence[i] = 0;
		}
	}

	// calc fitness
	for (i = 0; i < data->mpi.num_proc; i++) {
		fitness[i] = 0.0;
	}
	calcFEM(data);
	//outputOverlap(data, k);
	fitness[data->mpi.my_rank] = calcFitness(data);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(fitness, data->mosaic.fitness, data->mpi.num_proc, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	data->mosaic.FOM = 0.0;
	for(i = 0; i < data->mpi.num_proc; i++) {
		data->mosaic.FOM += data->mosaic.fitness[i];
	}
	data->mosaic.FOM /= data->mpi.num_proc;
	data->mosaic.FOM_history[0] = data->mosaic.FOM;

	
	if (data->mpi.my_rank == 0) {
		sprintf(filename, "Output/fitness.data");
		if ((fp = fopen(filename, "a+")) == NULL) {
			fprintf(stderr, "can't open file ( %s )\n", filename);
			exit(1);
		}
		
		fprintf(fp, "0: ");
		fprintf(stderr, "0: ");
		fprintf(fp, "%lf, %d pixels", data->mosaic.FOM, totalNN);
		fprintf(stderr, "%lf, %d pixels", data->mosaic.FOM, totalNN);
		if (data->mpi.num_proc != 1) {
			fprintf(fp, " / ");
			fprintf(stderr, " / ");
			for(i = 0; i < data->mpi.num_proc; i++) {
				fprintf(fp, "%lf,", data->mosaic.fitness[i]);
				fprintf(stderr, "%lf,", data->mosaic.fitness[i]);
			}
		}
		fprintf(fp, "\n");
		fprintf(stderr, "\n");
		fclose(fp);
	}

	change = 0;
	for (k=0; k<data->mosaic.DBSitr; k++) {
		if (data->mpi.my_rank == 0) {
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
		}

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(sequence, data->mosaic.DBSseq, totalNN, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		for (l=0; l<totalNN; l++) {
			// initialize fitness
			for (i = 0; i < data->mpi.num_proc; i++) {
				fitness[i] = 0.0;
			}

			// toggle one pixel state
			jj = sequence[l]%Nz;
			ii = sequence[l]/Nz;
			tmpi = (ysym != 0 ? data->mosaic.Ny-1-ii : ii);
			tmpj = (zsym != 0 ? data->mosaic.Nz-1-jj : jj);
			state = data->mosaic.BorW[ii][jj];
			data->mosaic.BorW[ii][jj] = 1-state;
			data->mosaic.BorW[ii][tmpj] = 1-state;
			data->mosaic.BorW[tmpi][jj] = 1-state;
			data->mosaic.BorW[tmpi][tmpj] = 1-state;
			if (yeqzsym != 0) {
				data->mosaic.BorW[jj][ii] = 1-state;
				data->mosaic.BorW[tmpj][ii] = 1-state;
				data->mosaic.BorW[jj][tmpi] = 1-state;
				data->mosaic.BorW[tmpj][tmpi] = 1-state;
			}

			// calc fitness
			calcFEM(data);
			//outputOverlap(data, k);
			fitness[data->mpi.my_rank] = calcFitness(data);

			// Share the simulation result
			MPI_Barrier(MPI_COMM_WORLD);
			/*	one itr. delayed confirmation in order to avoid the abuse of Barrier  //
			//	Confirmation right now is on [l-1] itr, not [l]						  */
			MPI_Allreduce(&change, &sum_change, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			if (sum_change != 0 && sum_change != data->mpi.num_proc) {
				if (data->mpi.my_rank == 0) {
					fprintf(stderr, "Pixel toggling errer among the processes\n");
				}
				exit(1);
			}
			MPI_Allreduce(fitness, data->mosaic.fitness, data->mpi.num_proc, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			temp = 0.0;
			for(i = 0; i < data->mpi.num_proc; i++) {
				temp += data->mosaic.fitness[i];
			}
			temp /= data->mpi.num_proc;

			change = (temp>data->mosaic.FOM ? 1 : 0);
			if (data->mpi.my_rank == 0) {
				sprintf(filename, "Output/fitness.data");
				if ((fp = fopen(filename, "a+")) == NULL) {
					fprintf(stderr, "can't open file ( %s )\n", filename);
					exit(1);
				}
				
				fprintf(fp, "%d: ", k*totalNN+l+1);
				fprintf(stderr, "%d: ", k*totalNN+l+1);
				fprintf(fp, "%lf,%lf,%d", temp, data->mosaic.FOM, change);
				fprintf(stderr, "%lf,%lf,%d", temp, data->mosaic.FOM, change);
				if (data->mpi.num_proc != 1) {
					fprintf(fp, " / ");
					fprintf(stderr, " / ");
					for(i = 0; i < data->mpi.num_proc; i++) {
						fprintf(fp, "%lf,", data->mosaic.fitness[i]);
						fprintf(stderr, "%lf,", data->mosaic.fitness[i]);
					}
				}
				fprintf(fp, "\n");
				fprintf(stderr, "\n");
				fclose(fp);
			}

			if (change==1) {
				data->mosaic.FOM = temp;
			}
			else {
				// go back to the original state
				jj = sequence[l]%Nz;
				ii = sequence[l]/Nz;
				tmpi = (ysym != 0 ? data->mosaic.Ny-1-ii : ii);
				tmpj = (zsym != 0 ? data->mosaic.Nz-1-jj : jj);
				state = data->mosaic.BorW[ii][jj];
				data->mosaic.BorW[ii][jj] = 1-state;
				data->mosaic.BorW[ii][tmpj] = 1-state;
				data->mosaic.BorW[tmpi][jj] = 1-state;
				data->mosaic.BorW[tmpi][tmpj] = 1-state;
				if (yeqzsym != 0) {
					data->mosaic.BorW[jj][ii] = 1-state;
					data->mosaic.BorW[tmpj][ii] = 1-state;
					data->mosaic.BorW[jj][tmpi] = 1-state;
					data->mosaic.BorW[tmpj][tmpi] = 1-state;
				}
			}

			if (data->mpi.my_rank == 0 && change==1) {
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
		}

		if (data->mpi.my_rank == 0) {			
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
		}

		data->mosaic.FOM_history[k+1] = data->mosaic.FOM;
		if (data->mosaic.FOM_history[k+1]-data->mosaic.FOM_history[k] < data->mosaic.termination + 1e-6) {
			fprintf(stderr, "FOM improvement is %lf, which is no greater than %lf", data->mosaic.FOM_history[k+1]-data->mosaic.FOM_history[k], data->mosaic.termination);
			break;
		}
	}

	if (data->mpi.my_rank == 0) {
		sprintf(filename, "Output/fitness.data");
		if ((fp = fopen(filename, "w")) == NULL) {
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

