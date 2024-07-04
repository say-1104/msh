#include "fem3d.h"
#include <omp.h>
#include <sys/time.h>
#include <linux/kernel.h>
#include <sys/sysinfo.h>

// Direct binary search for mozaic pattern
// xz-symmetry is considered

extern int DirectBinarySearch(DataTable *data, int kk){
	int i, j, k;
	int ii, jj;
	int ll, mm;
	long long int nn;
	
	bool xsym = (data->mozaic.symflag == 1 || data->mozaic.symflag == 4), zsym = (data->mozaic.symflag == 4);
	int state;
	int totalNN;
	int starti, startj;

	int total = 0, Nx, Nz, tmpi, tmpj;
	char filename[256];
	FILE *fp;
	int *randtotal;
	int **randtotalMPI;

	fprintf(stderr, "Direct Binary Search\n");

	ALLOCATION(data->mozaic.efficiency_history, double, data->mozaic.DBSitr+1);

	for(i = 0; i < data->mozaic.DBSitr; i++) {
		data->mozaic.efficiency_history[i] = 0;
	}

	// File initialization
	if (data->mozaic.starttotal == 0) {
		sprintf(filename, "IntermediateMatrix.csv");
		if ((fp = fopen(filename, "w")) == NULL) {
			fprintf(stderr, "can't open file ( %s )\n", filename);
			exit(1);
		}
		fclose(fp);

		sprintf(filename, "log-TT.data");
		if ((fp = fopen(filename, "w")) == NULL) {
			fprintf(stderr, "can't open file ( %s )\n", filename);
			exit(1);
		}
		fclose(fp);

		if (data->mozaic.DBSitr != 0) {
			sprintf(filename, "sequence.txt");
			if ((fp = fopen(filename, "w")) == NULL) {
				fprintf(stderr, "can't open file ( %s )\n", filename);
				exit(1);
			}
			fclose(fp);
		}
	}

	// calculate the characteristics of first structure
	analysis2(data, 0);

	calcFitness(data, 0, 0);

	// -------------------------------------------

    if (xsym){
		if(data->mozaic.Nx%2 == 0){
			Nx = data->mozaic.Nx/2;
		}else{
			Nx = data->mozaic.Nx/2+1;
		}
	} 
	
	if (zsym) {
		// z-symmetry is considered
		if(data->mozaic.mozDev == 5){
			//waveguide lense
			Nz = data->mozaic.Nz/3;
		}else{
			if(data->mozaic.Nz%2 == 0){
				Nz = data->mozaic.Nz/2;
			}else{
				Nz = data->mozaic.Nz/2+1;
			}
		}
	}

	totalNN = Nx*Nz;

	if (data->mozaic.PixelReduction == 1) {
		totalNN = 0;
		for (i=0; i<Nx; i++) {
			for (j=0; j<Nz; j++) {
				if (data->mozaic.BorW[kk][i][j] != -1) {
					totalNN++;
				}
			}
		}
	}

	if (data->mozaic.mozDev == 6) {
		if (Nx != Nz) {
			fprintf(stderr, "this structure can't implement Y=Z symmetry, because Nx isn't equal to Nz");
			exit(1);
		}

		totalNN = 0;
		for (i=0; i<Nx; i++) {
			for (j=0; j<Nz; j++) {
				if (j <= i) {
					if (data->mozaic.BorW[kk][i][j] != -1) {
						totalNN++;
					}
				}
			}
		}
	}

	fprintf(stderr, "Nx = %d, Nz = %d\n", Nx, Nz);
	
	ALLOCATION(randtotal, int, totalNN);
	for(i = 0; i < totalNN; i++) {
		randtotal[i] = 0;
	}

	ALLOCATION(randtotalMPI, int *, data->par.num_proc);
	for(mm = 0; mm < data->par.num_proc; mm++) {
		ALLOCATION(randtotalMPI[mm], int, totalNN);
	}
	
	k=0;
	for(mm = 0; mm < data->par.num_proc; mm++) {
		if(data->par.myid == 0) {
			if (data->mozaic.mozDev == 6) {
				for(i = 0; i < Nx*Nz; i++) {
					jj = i%Nz;
					ii = (i-jj)/Nz;
					if (jj <= ii) {
						if (data->mozaic.BorW[kk][ii][jj] != -1) {
							randtotalMPI[mm][k] = i;
							k++;
						}
					}
				}
			} else if (data->mozaic.PixelReduction == 1) {
				for(i = 0; i < Nx*Nz; i++) {
					jj = i%Nz;
					ii = (i-jj)/Nz;
					if (data->mozaic.BorW[kk][ii][jj] != -1) {
						randtotalMPI[mm][k] = i;
						k++;
					}
				}
			} else {
				for(i = 0; i < totalNN; i++) {
					randtotalMPI[mm][i] = i;
				}
			}
		} else {
			for(i = 0; i < totalNN; i++) {
				randtotalMPI[mm][i] = 0;
			}
		}
	}

	total = data->mozaic.starttotal;

	// find the no. to resume
	if (total != 0) {
		sprintf(filename, "IntermediateMatrix.csv");
		if ((fp = fopen(filename, "r")) == NULL) {
			fprintf(stderr, "can't open file ( %s )\n", filename);
			exit(1);
		}
		
		// take note that restarting point is already added (original, previous pixels) + (restarting original). remove 2 originals
		total = -2;
		while ((k = fgetc(fp)) != EOF) {
			if (k == '\n') total++;
		}
		data->mozaic.starttotal = total;

		fprintf(stderr, "DBS restart from %d th itr(%d th pixel)\n",
			data->mozaic.starttotal / totalNN, data->mozaic.starttotal);

		k = data->mozaic.starttotal / totalNN;
	}
	else {
		k = 0;
		data->mozaic.efficiency_history[0] = data->mozaic.efficiency;
	}

	for (; k < data->mozaic.DBSitr; k++) {
	
		fprintf(stderr, "%d th DBS iteration\n", k);

		starti = 0; startj = 0;

		fprintf(stderr, "Nx = %d, Nz = %d, total = %d\n", Nx, Nz, totalNN);

		shuffle(randtotalMPI[data->par.myid], totalNN);

		if(data->par.myid == 0){
			if (k == data->mozaic.starttotal / totalNN && data->mozaic.starttotal != 0) {
				fprintf(stderr, "Restore sequence\n");	
				sprintf(filename, "sequence.txt");
				if ((fp = fopen(filename, "r")) == NULL) {
					fprintf(stderr, "can't open file ( %s )\n", filename);
					exit(1);
				}

				for (i=0; i<(k+1)*totalNN; i++) {
					fscanf(fp, "%d %d", &j, &randtotalMPI[0][i%totalNN]);
				}
				fclose(fp);
				for(i = 0; i < totalNN; i++) {
					jj = randtotalMPI[data->par.myid][i]%Nz;
					ii = (randtotalMPI[data->par.myid][i]-jj)/Nz;
					if (i == data->mozaic.starttotal%totalNN) {
						fprintf(stderr, "-- start point --\n");		
					}
					fprintf(stderr, "randtotal %d Nx %d Nz %d\n", randtotalMPI[0][i], ii, jj);					
				}
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);

		// all pixel ordering is assembled to randtotal
		// randtotalMPI is zero, for data->par.myid != 0
		MPI_Allreduce(randtotalMPI[data->par.myid],
				randtotal, totalNN, 
				MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		if(data->par.myid == 0){
			if (!(k == data->mozaic.starttotal / totalNN && data->mozaic.starttotal != 0)) {
				sprintf(filename, "sequence.txt");
				if ((fp = fopen(filename, "a+")) == NULL) {
					fprintf(stderr, "can't open file ( %s )\n", filename);
					exit(1);
				}
				for(i = 0; i < totalNN; i++) {
					jj = randtotal[i]%Nz;
					ii = (randtotal[i]-jj)/Nz;
					fprintf(stderr, "myid %d randtotal %d Nx %d Nz %d\n", data->par.myid, 
							randtotal[i], ii, jj);
					
					fprintf(fp, "%d %d\n", k*totalNN+i, randtotal[i]);
				}
				fclose(fp);
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);


		// Originally, initialization for loop.
		i = 0;
		if (k == data->mozaic.starttotal / totalNN && data->mozaic.starttotal != 0) {
			i = data->mozaic.starttotal % totalNN;
			fprintf(stderr, "Restart from %d th pixel\n", i);
		}
		for(; i < totalNN; i++) {
			// DBS
			fprintf(stderr, "%d %d\n", i, randtotal[i]); 

			jj = randtotal[i]%Nz;
			ii = (randtotal[i]-jj)/Nz;
			tmpi = (xsym ? data->mozaic.Nx-1-ii : ii);
			tmpj = (zsym ? data->mozaic.Nz-1-jj : jj);
			state = data->mozaic.BorW[kk][ii][jj];

			// fprintf(stderr, "%d %d %d %d\n", ii, jj, tmpi, tmpj); 

			// Inverting the pixel 
			data->mozaic.BorW[kk][ii][jj] = 1-state;
			data->mozaic.BorW[kk][tmpi][jj] = 1-state;
			data->mozaic.BorW[kk][ii][tmpj] = 1-state;
			data->mozaic.BorW[kk][tmpi][tmpj] = 1-state;

			if (data->mozaic.mozDev == 6) {
				data->mozaic.BorW[kk][jj][ii] = 1-state;
				data->mozaic.BorW[kk][jj][tmpi] = 1-state;
				data->mozaic.BorW[kk][tmpj][ii] = 1-state;
				data->mozaic.BorW[kk][tmpj][tmpi] = 1-state;
			}

			fprintf(stderr, "i,j %d %d, ii,jj,tmpi,tmpj %d %d %d %d\n", 
					i, j, ii, jj, tmpi, tmpj);


			MPI_Barrier(MPI_COMM_WORLD);
			// transmission of the structure
			analysis2(data, 1);
			
			calcFitness(data, 1, total+1);

			// -------------------------------
			
			if (data->mozaic.efficiency < data->mozaic.fittemp){
				data->mozaic.efficiency = data->mozaic.fittemp;

				sprintf(filename, "GuessedMatrix-middle.csv");
				if ((fp = fopen(filename, "w")) == NULL) {
					fprintf(stderr, "can't open file ( %s )\n", filename);
					exit(1);
				}
				
				for(mm = 0; mm < data->mozaic.Nx; mm++) {
					for(nn = 0; nn < data->mozaic.Nz; nn++) {
						if(mm == data->mozaic.Nx-1 && nn == data->mozaic.Nz-1){
							fprintf(fp, "%d", data->mozaic.BorW[kk][mm][nn]);
						} else {
							fprintf(fp, "%d,", data->mozaic.BorW[kk][mm][nn]);
						}
					}
				}
				fprintf(fp, "\n");
				fprintf(fp, "ii %d jj %d total %d\n", ii, jj, total);
				
				fclose(fp);

			} else {
				//fprintf(stderr, "(%d %d %d %d)\n", kk, ii, jj, tmpi); 
				
				// return the pixel to original	
				jj = randtotal[i]%Nz;
				ii = (randtotal[i]-jj)/Nz;
				tmpi = (xsym ? data->mozaic.Nx-1-ii : ii);
				tmpj = (zsym ? data->mozaic.Nz-1-jj : jj);
				state = data->mozaic.BorW[kk][ii][jj];

				data->mozaic.BorW[kk][ii][jj] = 1-state;
				data->mozaic.BorW[kk][tmpi][jj] = 1-state;
				data->mozaic.BorW[kk][ii][tmpj] = 1-state;
				data->mozaic.BorW[kk][tmpi][tmpj] = 1-state;
				if (data->mozaic.mozDev == 6) {
					data->mozaic.BorW[kk][jj][ii] = 1-state;
					data->mozaic.BorW[kk][jj][tmpi] = 1-state;
					data->mozaic.BorW[kk][tmpj][ii] = 1-state;
					data->mozaic.BorW[kk][tmpj][tmpi] = 1-state;
				}
			}
					
			total++;
		} // end of loop i

		data->mozaic.efficiency_history[k+1] = data->mozaic.efficiency;

		if (fabs(data->mozaic.efficiency_history[k+1] - data->mozaic.efficiency_history[k]) < data->mozaic.termination) {
			k += 1;		// for file output
			break;
		}

	} // end of loop k


	sprintf(filename, "log-TT.data");
	if ((fp = fopen(filename, "a+")) == NULL) {
		fprintf(stderr, "can't open file ( %s )\n", filename);
		exit(1);
	}
	for (i=0; i <= k; i++) {
		fprintf(fp, "%d-th itr: %lf, ", i, data->mozaic.efficiency_history[i]);
	}
	fprintf(fp, "\n");
	fclose(fp);

	fprintf(stderr, "total = %d\n", total); 

	sprintf(filename, "GuessedMatrix.csv");
	if ((fp = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "can't open file ( %s )\n", filename);
		exit(1);
	}
	 
	for(i = 0; i < data->mozaic.Nx; i++) {
		for(j = 0; j < data->mozaic.Nz; j++) {
			if(i == data->mozaic.Nx-1 && j == data->mozaic.Nz-1){
				fprintf(fp, "%d", data->mozaic.BorW[kk][i][j]);
			}else{
				fprintf(fp, "%d,", data->mozaic.BorW[kk][i][j]);
			}
		}
	}

	fclose(fp);

}