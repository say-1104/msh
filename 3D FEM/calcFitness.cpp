#include "fem3d.h"
#include <omp.h>
#include <sys/time.h>
#include <linux/kernel.h>
#include <sys/sysinfo.h>

void calcFitness(DataTable *data, int flag, int total) {
	int i, j, k;
	int ii, jj, kk;
	int ll, mm, nn;
	int oo, pp, qq;
	//  long long int nn, ne, ns;
	//  int nrhs, count;

	double wl, temp, temp1, temp2;
	double fit1, fit2;
	double alpha;
	
	std::complex<double> cj = std::complex<double>(0, 1);
	std::complex<double> c_a, c_b;

	double T3_ideal, T1_ideal;

	double totalT, dT1, dT2, dT3, ratio1, ratio2, ratio3;

	char filename[256];
	FILE *fp;
	double **tempdata;
	double** im_tempdata;		// 複素振幅を伝達するために用意

	FEM *fem = &(data->fem);
	Param *par = &(data->par);
	PortInfo *port;

	ALLOCATION(tempdata, double*, 2*data->mozaic.numModes);
	ALLOCATION(im_tempdata, double*, 2*data->mozaic.numModes);

	for (i = 0; i < 2*data->mozaic.numModes; i++) {
		ALLOCATION(tempdata[i], double, data->par.num_proc);
		ALLOCATION(im_tempdata[i], double, data->par.num_proc);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// data->mozaic.TT[i][j][k] is calculated in calcPower2 in analysis2.cpp
	// i: wavelength label (when using MPI, it is always 0)
	// j: mode No
	// k: currently, always 0

	if (data->mozaic.mozDev == 1){
		// fitness function of mode divider
		T1_ideal = data->mozaic.gamma;
		T3_ideal = 1.0-data->mozaic.gamma;
		alpha = 0.3;
		temp = temp1 = temp2 = 0.0;

		fit1 = (fabs(data->mozaic.TT[0][3]-T3_ideal)
			+fabs(data->mozaic.TT[0][1]-T1_ideal))
		/(data->mozaic.TT[0][3]+data->mozaic.TT[0][1]);

		fit2 = fabs(data->mozaic.TT[0][3]/data->mozaic.TT[0][1]
			-T3_ideal/T1_ideal)
		/(data->mozaic.TT[0][3]+data->mozaic.TT[0][1]);

		temp1 = (1.0-alpha)*fit1+alpha*fit2;
		// lower fundamental mode
		tempdata[0][data->par.myid] = data->mozaic.TT[0][1];
		// upper higher order mode
		tempdata[1][data->par.myid] = data->mozaic.TT[0][3];

		// all lower TE0 is assembled to tempdata[2]
		MPI_Allreduce(tempdata[0],
			tempdata[2],
			data->par.num_proc,
			MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);//プロセス間の通信をする

		// all higher TE1 is assembled to tempdata[3]
		MPI_Allreduce(tempdata[1],
			tempdata[3],
			data->par.num_proc,
			MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);//プロセス間の通信をする

	} else if (data->mozaic.mozDev == 2){
		// fitness function of mode MUX
		alpha = 0.1;
		data->mozaic.gamma = 1.0;
		temp = temp1 = temp2 = 0.0;

		fit1 = (fabs(data->mozaic.TT[0][2]-data->mozaic.gamma)
			+fabs(data->mozaic.TT[1][3]-data->mozaic.gamma))
			/(data->mozaic.TT[0][2]+data->mozaic.TT[1][3]);

		fit2 = fabs(data->mozaic.TT[1][3]/data->mozaic.TT[0][2]
			-data->mozaic.gamma)
			/(data->mozaic.TT[0][2]+data->mozaic.TT[1][3]);

		temp1 = (1.0-alpha)*fit1+alpha*fit2;
		// TE0 to TE0
		tempdata[0][data->par.myid] = data->mozaic.TT[0][2];
		// TE0 to TE1
		tempdata[1][data->par.myid] = data->mozaic.TT[1][3];

		// all TE0 is assembled to tempdata[2]
		MPI_Allreduce(tempdata[0],
				tempdata[2],
				data->par.num_proc,
				MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);//プロセス間の通信をする

		// all TE1 is assembled to tempdata[3]
		MPI_Allreduce(tempdata[1],
				tempdata[3],
				data->par.num_proc,
				MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);//プロセス間の通信をする

	} else if (data->mozaic.mozDev == 3){
		// fitness function of mode EX
		alpha = 0.3;
		data->mozaic.gamma = 1.0; // ideal transmission
		temp = temp1 = temp2 = 0.0;

		// TT[0]: the transmission of TE0 input
		// TT[1]: the transmission of TE1 input
		// TT[2]: the transmission of TE2 input

		if(data->mozaic.numModes == 2){
			// TT[0]: the transmission of TE0->TE1
			ii = 3;
			// TT[1]: the transmission of TE1->TE0
			jj = 2;

			totalT = data->mozaic.TT[0][ii]+data->mozaic.TT[1][jj];
			dT1 = fabs(data->mozaic.TT[0][ii]-data->mozaic.gamma);
			dT2 = fabs(data->mozaic.TT[1][ii]-data->mozaic.gamma);

			ratio1 = data->mozaic.TT[1][jj]/data->mozaic.TT[0][ii];

			fit1 = (dT1+dT2)/totalT;

			fit2 = fabs(ratio1-data->mozaic.gamma)/totalT;

		} else if(data->mozaic.numModes == 3) {
			// TT[0]: the transmission of TE0->TE1
			// TT[1]: the transmission of TE1->TE2
			// TT[2]: the transmission of TE2->TE0
			ii = 4;
			jj = 5;
			kk = 3;

			totalT = data->mozaic.TT[0][ii]+data->mozaic.TT[1][jj]
				+data->mozaic.TT[2][kk];
			dT1 = fabs(data->mozaic.TT[0][ii]-data->mozaic.gamma);
			dT2 = fabs(data->mozaic.TT[1][ii]-data->mozaic.gamma);
			dT3 = fabs(data->mozaic.TT[2][ii]-data->mozaic.gamma);

			ratio1 = data->mozaic.TT[1][jj]/data->mozaic.TT[0][ii];
			ratio2 = data->mozaic.TT[2][kk]/data->mozaic.TT[1][jj];
			ratio3 = data->mozaic.TT[0][ii]/data->mozaic.TT[2][kk];

			fit1 = (dT1+dT2+dT3)/totalT;

			fit2 = (fabs(ratio1-data->mozaic.gamma)+fabs(ratio2-data->mozaic.gamma)
				+fabs(ratio3-data->mozaic.gamma))/totalT;

			// TE2 to TE0
			tempdata[4][data->par.myid] = data->mozaic.TT[2][kk];

			// all TE2 is assembled to tempdata[5]
			MPI_Allreduce(tempdata[4],
				tempdata[5],
				data->par.num_proc,
				MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);//プロセス間の通信をする
		}

		fprintf(stderr, "%d temp = %lf\n", mm, (1.0-alpha)*fit1+alpha*fit2);

		temp1 = (1.0-alpha)*fit1+alpha*fit2;

		// TE0 to TE1
		tempdata[0][data->par.myid] = data->mozaic.TT[0][ii];
		// TE1 to TE0
		tempdata[1][data->par.myid] = data->mozaic.TT[1][jj];

		// all TE0 is assembled to tempdata[2]
		MPI_Allreduce(tempdata[0],
				tempdata[2],
				data->par.num_proc,
				MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);//プロセス間の通信をする

		// all TE1 is assembled to tempdata[3]
		MPI_Allreduce(tempdata[1],
				tempdata[3],
				data->par.num_proc,
				MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);//プロセス間の通信をする

	} else if (data->mozaic.mozDev == 4){
		// fitness function of mode separator
		alpha = 0.3;
		data->mozaic.gamma = 1.0; // ideal transmission
		temp = temp1 = temp2 = 0.0;

		// TT[0]: the transmission of TE0 input
		// TT[1]: the transmission of TE1 input
		// TT[2]: the transmission of TE2 input

		if (data->mozaic.numModes == 2) {
			// TT[0]: the transmission of TE0-> lower TE0
			ii = 2;
			// TT[1]: the transmission of TE1-> upper TE1
			jj = 4;

			totalT = data->mozaic.TT[0][ii]+data->mozaic.TT[1][jj];
			dT1 = fabs(data->mozaic.TT[0][ii]-data->mozaic.gamma);
			dT2 = fabs(data->mozaic.TT[1][ii]-data->mozaic.gamma);

			ratio1 = data->mozaic.TT[1][jj]/data->mozaic.TT[0][ii];

			fit1 = (dT1+dT2)/totalT;

			fit2 = fabs(ratio1-data->mozaic.gamma)/totalT;

		} else if(data->mozaic.numModes == 3){
			// TT[0]: the transmission of TE0-> upper TE0
			// TT[1]: the transmission of TE1-> lower TE1
			// TT[2]: the transmission of TE2-> middle TE2
			ii = 8;
			jj = 4;
			kk = 7;

			totalT = data->mozaic.TT[0][ii]+data->mozaic.TT[1][jj]
				+data->mozaic.TT[2][kk];
			dT1 = fabs(data->mozaic.TT[0][ii]-data->mozaic.gamma);
			dT2 = fabs(data->mozaic.TT[1][ii]-data->mozaic.gamma);
			dT3 = fabs(data->mozaic.TT[2][ii]-data->mozaic.gamma);

			ratio1 = data->mozaic.TT[1][jj]/data->mozaic.TT[0][ii];
			ratio2 = data->mozaic.TT[2][kk]/data->mozaic.TT[1][jj];
			ratio3 = data->mozaic.TT[0][ii]/data->mozaic.TT[2][kk];

			fit1 = (dT1+dT2+dT3)/totalT;

			fit2 = (fabs(ratio1-data->mozaic.gamma)+fabs(ratio2-data->mozaic.gamma)
				+fabs(ratio3-data->mozaic.gamma))/totalT;

			// TE2 to TE2
			tempdata[4][data->par.myid] = data->mozaic.TT[2][kk];

			// all TE2 is assembled to tempdata[5]
			MPI_Allreduce(tempdata[4],
				tempdata[5],
				data->par.num_proc,
				MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);//プロセス間の通信をする
		}

		fprintf(stderr, "%d temp = %lf\n", mm, (1.0-alpha)*fit1+alpha*fit2);

		temp1 = (1.0-alpha)*fit1+alpha*fit2;

		// TE0 to TE1
		tempdata[0][data->par.myid] = data->mozaic.TT[0][ii];
		// TE1 to TE0
		tempdata[1][data->par.myid] = data->mozaic.TT[1][jj];

		// all TE0 is assembled to tempdata[2]
		MPI_Allreduce(tempdata[0],
				tempdata[2],
				data->par.num_proc,
				MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);//プロセス間の通信をする

		// all TE1 is assembled to tempdata[3]
		MPI_Allreduce(tempdata[1],
				tempdata[3],
				data->par.num_proc,
				MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);//プロセス間の通信をする

	} else if (data->mozaic.mozDev == 5) {
		// fitness function of wavegudie lense
		alpha = 0.3;
		data->mozaic.gamma = 1.0; // ideal transmission
		temp = temp1 = temp2 = 0.0;

		// TT[0]: the transmission of TE0 input
		// TT[1]: the transmission of TE1 input
		// TT[2]: the transmission of TE2 input

		if (data->mozaic.numModes == 1) {
			// TT[0]: the transmission of TE0->TE0
			ii = 1;

			totalT = data->mozaic.TT[0][ii];
			dT1 = fabs(data->mozaic.TT[0][ii]-data->mozaic.gamma);

			fit1 = (dT1)/totalT;

			fit2 = 0.0;

		} else if (data->mozaic.numModes == 2) {
			// TT[0]: the transmission of TE0-> TE0
			ii = 2;
			// TT[1]: the transmission of TE1-> TE1
			jj = 3;

			totalT = data->mozaic.TT[0][ii]+data->mozaic.TT[1][jj];
			dT1 = fabs(data->mozaic.TT[0][ii]-data->mozaic.gamma);
			//dT2 = fabs(data->mozaic.TT[1][ii]-data->mozaic.gamma);
			dT2 = fabs(data->mozaic.TT[1][jj]-data->mozaic.gamma);//11/03偏光
			fprintf(stderr, "dT1 = %lf\ndT2 = %lf\n",dT1,dT2);
			ratio1 = data->mozaic.TT[1][jj]/data->mozaic.TT[0][ii];

			fit1 = (dT1+dT2)/totalT;

			fit2 = fabs(ratio1-data->mozaic.gamma)/totalT;

		} else if(data->mozaic.numModes == 3){
			// TT[0]: the transmission of TE0-> TE0
			// TT[1]: the transmission of TE1-> TE1
			// TT[2]: the transmission of TE2-> TE2
			ii = 3;
			jj = 4;
			kk = 5;

			totalT = data->mozaic.TT[0][ii]+data->mozaic.TT[1][jj]
				+data->mozaic.TT[2][kk];
			dT1 = fabs(data->mozaic.TT[0][ii]-data->mozaic.gamma);
			//dT2 = fabs(data->mozaic.TT[1][ii]-data->mozaic.gamma);
			dT2 = fabs(data->mozaic.TT[1][jj]-data->mozaic.gamma);
			//dT3 = fabs(data->mozaic.TT[2][ii]-data->mozaic.gamma);
			dT3 = fabs(data->mozaic.TT[2][kk]-data->mozaic.gamma);
			fprintf(stderr, "dT1 = %lf\ndT2 = %lf\ndT3 = %lf\n",dT1,dT2,dT3);
			ratio1 = data->mozaic.TT[1][jj]/data->mozaic.TT[0][ii];
			ratio2 = data->mozaic.TT[2][kk]/data->mozaic.TT[1][jj];
			ratio3 = data->mozaic.TT[0][ii]/data->mozaic.TT[2][kk];

			fit1 = (dT1+dT2+dT3)/totalT;

			fit2 = (fabs(ratio1-data->mozaic.gamma)+fabs(ratio2-data->mozaic.gamma)
				+fabs(ratio3-data->mozaic.gamma))/totalT;

			tempdata[1][data->par.myid] = data->mozaic.TT[1][jj];

			// all TE1 is assembled to tempdata[3]
			MPI_Allreduce(tempdata[1], tempdata[3],data->par.num_proc, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);//プロセス間の通信をする

			// TE2 to TE2
			tempdata[4][data->par.myid] = data->mozaic.TT[2][kk];
			fprintf(stderr, "tempdata[4][%d] = %lf\n", data->par.myid, data->mozaic.TT[2][kk]);

			// all TE2 is assembled to tempdata[5]
			MPI_Allreduce(tempdata[4],tempdata[5],data->par.num_proc,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);//プロセス間の通信をする
		}

		fprintf(stderr, "%d temp = %lf\n", mm, (1.0-alpha)*fit1+alpha*fit2);

		temp1 = (1.0-alpha)*fit1+alpha*fit2;

		// TE0 to TE0
		tempdata[0][data->par.myid] = data->mozaic.TT[0][ii];

		// all TE0 is assembled to tempdata[2]
		MPI_Allreduce(tempdata[0],tempdata[2],data->par.num_proc,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);//プロセス間の通信をする

		if(data->mozaic.numModes == 2){
			// TE1 to TE1
			tempdata[1][data->par.myid] = data->mozaic.TT[1][jj];

			// all TE1 is assembled to tempdata[3]
			MPI_Allreduce(tempdata[1], tempdata[3],data->par.num_proc, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);//プロセス間の通信をする
		}
	} else if (data->mozaic.mozDev == 6) {
		// fitness function of Waveguide Crossing
		temp = temp1 = temp2 = 0.0;
		if (data->mozaic.numModes == 4) {
			// arithmetic mean
			temp1 = 0;
			for (i = 0; i < data->mozaic.numModes; i++) {
				temp1 += data->mozaic.TT[i][i+data->mozaic.numModes];
			}
			temp1 /= (double)data->mozaic.numModes;
			// // geometric mean
			// temp1 = 1;
			// for (i = 0; i < data->mozaic.numModes; i++) {
			// 	temp1 *= data->mozaic.TT[i][i+data->mozaic.numModes];
			// }
			// pow(temp1, data->mozaic.numModes);
			// // harmonic mean
			// temp1 = 0;
			// for (i = 0; i < data->mozaic.numModes; i++) {
			// 	temp1 += 1.0 / data->mozaic.TT[i][i+data->mozaic.numModes];
			// }
			// temp1 = (double)data->mozaic.numModes / temp1;
			// // minimum
			// temp1 = 100;
			// for (i = 0; i < data->mozaic.numModes; i++) {
			// 	if (temp1 > data->mozaic.TT[i][i+data->mozaic.numModes]) {
			// 		temp1 = data->mozaic.TT[i][i+data->mozaic.numModes];
			// 	}
			// }

			for (i = 0; i < data->mozaic.numModes; i++) {
				tempdata[i][data->par.myid] = data->mozaic.TT[i][i+data->mozaic.numModes];
			}

			// each data is assembled
			for (i = 0; i < data->mozaic.numModes; i++) {
				MPI_Allreduce(tempdata[i],
						tempdata[i+data->mozaic.numModes],
						data->par.num_proc,
						MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			}
		}
	} else if (data->mozaic.mozDev == 7) {
		// fitness function of nxn Power Splitter
		// temp = temp1 = temp2 = 0.0;
		// if (data->mozaic.numModes == 2) {
		// 	temp1 +=   1 - (fabs(data->mozaic.TT[0][1]-data->mozaic.gamma) + fabs(data->mozaic.TT[0][2]-(1.0-data->mozaic.gamma)))
		// 				- fabs(data->mozaic.TT[0][1] - data->mozaic.TT[0][2]);

		// 	tempdata[0][data->par.myid] = data->mozaic.TT[0][1];
		// 	tempdata[1][data->par.myid] = data->mozaic.TT[0][2];

		// 	// each data is assembled
		// 	MPI_Allreduce(tempdata[0],
		// 			tempdata[2],
		// 			data->par.num_proc,
		// 			MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		// 	MPI_Allreduce(tempdata[1],
		// 			tempdata[3],
		// 			data->par.num_proc,
		// 			MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		// }

		temp = temp1 = temp2 = 0.0;
		if (data->mozaic.numModes == 2) {
			// 回転させた複素数との距離を取る。
			c_a = data->mozaic.c_TT[0][1] * cj;
			c_b = data->mozaic.c_TT[0][2];

			// (1-power_total) = insertion loss + distance between amplitudes (complex) of 2 port
			// temp1 = 1 - (pow(std::abs(c_a), 2) + pow(std::abs(c_b), 2)) 
			// 			+ std::abs(c_a - c_b);
			temp1 = (pow(std::abs(c_a), 2) + pow(std::abs(c_b), 2)) 
						- std::abs(c_a - c_b);

			for (i = 0; i < data->mozaic.numModes; i++) {
				tempdata[i][data->par.myid] = std::real(data->mozaic.c_TT[0][i+1]);
				im_tempdata[i][data->par.myid] = std::imag(data->mozaic.c_TT[0][i+1]);
			}
			
			// each data is assembled
			for (i = 0; i < data->mozaic.numModes; i++) {
				MPI_Allreduce(tempdata[i],
						tempdata[i+data->mozaic.numModes],
						data->par.num_proc,
						MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
				MPI_Allreduce(im_tempdata[i],
						im_tempdata[i+data->mozaic.numModes],
						data->par.num_proc,
						MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			}
		}
	}

	// for MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&temp1,
		&temp,
		1,
		MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	data->mozaic.fittemp = temp/(double)data->par.num_proc;

	if(flag == 0) {
		data->mozaic.efficiency = temp/(double)data->par.num_proc;

		fprintf(stderr, "%d efficiency = %lf\n", total,
			data->mozaic.efficiency);
	}

	fprintf(stderr, "myid = %d, %d efficiency = %lf, Tlamda = %lf\n",
		data->par.myid, total,
		data->mozaic.efficiency, tempdata[0][data->par.myid]);

	if(data->par.myid == 0) {
		sprintf(filename, "log-TT.data");
		if ((fp = fopen(filename, "a+")) == NULL) {
			fprintf(stderr, "can't open file ( %s )\n", filename);
			exit(1);
		}
		fprintf(fp, "%d\t%lf\t%lf\t/\t", total, data->mozaic.efficiency,
			data->mozaic.fittemp);

		if(data->mozaic.mozDev == 1){
			for(mm = 0; mm < data->par.num_proc; mm++) {
				// for mode divider
				fprintf(fp, "%lf %lf ",  tempdata[2][mm], tempdata[3][mm]);
			}
		} else if (data->mozaic.mozDev == 2){
			for(mm = 0; mm < data->par.num_proc; mm++) {
				// for mode MUX
				fprintf(fp, "%lf %lf ",  tempdata[2][mm], tempdata[3][mm]);
			}
		} else if (data->mozaic.mozDev == 3 || data->mozaic.mozDev == 4
			|| data->mozaic.mozDev == 5){
			for(mm = 0; mm < data->par.num_proc; mm++) {
				if(data->mozaic.numModes == 1){
					// for 1 mode EX
					fprintf(fp, "%lf ",  tempdata[2][mm]);
				} else if (data->mozaic.numModes == 2){
					// for 2 mode EX
					fprintf(fp, "%lf %lf ",  tempdata[2][mm], tempdata[3][mm]);
				} else if (data->mozaic.numModes == 3){
					// for 3 mode EX
					fprintf(fp, "%lf\t%lf\t%lf\t",
						tempdata[2][mm], tempdata[3][mm], tempdata[5][mm]);
				}
			}
		} else if (data->mozaic.mozDev == 6) {
			for(mm = 0; mm < data->par.num_proc; mm++) {
				if (data->mozaic.numModes == 4){
					fprintf(fp, "%lf\t%lf\t%lf\t%lf\t",
						tempdata[4][mm], tempdata[5][mm], tempdata[6][mm], tempdata[7][mm]);
				}
			}
		} else if (data->mozaic.mozDev == 7) {
			if (data->par.numStructures != 1) {
				for(mm = 0; mm < data->par.numofwlMPI; mm++) {
					fprintf(fp, "wl%d:\t", mm); 
					for(nn = 0; nn < data->par.numStructures; nn++) {
						if (mm*data->par.numStructures + nn >= data->par.num_proc) break;
						if (data->mozaic.numModes == 2) {
							fprintf(fp, "%lf\t%lf\t%lf\t%lf\t",
								tempdata[2][mm*data->par.numStructures + nn], im_tempdata[2][mm*data->par.numStructures + nn],
								tempdata[3][mm*data->par.numStructures + nn], im_tempdata[3][mm*data->par.numStructures + nn]);
						}
					}
				}
			}
			else {
				for(mm = 0; mm < data->par.num_proc; mm++) {
					if (data->mozaic.numModes == 2) {
						fprintf(fp, "%lf\t%lf\t%lf\t%lf\t",
							tempdata[2][mm], im_tempdata[2][mm],
							tempdata[3][mm], im_tempdata[3][mm]);
					}
				}
			}
		}

		fprintf(fp, "\n");
		fclose(fp);
	}

	// output on the way characteristics of SSC
	if(data->par.myid == 0){
		sprintf(filename, "IntermediateMatrix.csv");
		if ((fp = fopen(filename, "a+")) == NULL) {
			fprintf(stderr, "can't open file ( %s )\n", filename);
			exit(1);
		}

		for(mm = 0; mm < data->mozaic.Nx; mm++) {
			for(nn = 0; nn < data->mozaic.Nz; nn++) {
				fprintf(fp, "%d,", data->mozaic.BorW[0][mm][nn]);
			}
		}

		if (data->mozaic.mozDev == 1) {
			for(mm = 0; mm < data->par.num_proc; mm++) {
				// for mode divider
				fprintf(fp, "%lf %lf ",  tempdata[2][mm], tempdata[3][mm]);
			}
		} else if (data->mozaic.mozDev == 2) {
			for(mm = 0; mm < data->par.num_proc; mm++) {
				// for mode MUX
				fprintf(fp, "%lf %lf ",  tempdata[2][mm], tempdata[3][mm]);
			}
		} else if(data->mozaic.mozDev == 3 || data->mozaic.mozDev == 4
			|| data->mozaic.mozDev == 5) {
			for(mm = 0; mm < data->par.num_proc; mm++) {
				if(data->mozaic.numModes == 1){
					// for 1 mode EX
					fprintf(fp, "%lf ",  tempdata[2][mm]);
				} else if(data->mozaic.numModes == 2) {
					// for 2 mode EX
					fprintf(fp, "%lf %lf ",  tempdata[2][mm], tempdata[3][mm]);
				} else if(data->mozaic.numModes == 3) {
					// for 3 mode EX
					fprintf(fp, "%lf %lf  %lf",
						tempdata[2][mm], tempdata[3][mm], tempdata[5][mm]);
				}
			}
		} else if (data->mozaic.mozDev == 6) {
			for(mm = 0; mm < data->par.num_proc; mm++) {
				if (data->mozaic.numModes == 4){
					fprintf(fp, "%lf %lf %lf %lf ",
						tempdata[4][mm], tempdata[5][mm], tempdata[6][mm], tempdata[7][mm]);
				}
			}
		} else if (data->mozaic.mozDev == 7) {
			for(mm = 0; mm < data->par.num_proc; mm++) {
				if (data->mozaic.numModes == 2){
					fprintf(fp, "%lf %lf ", tempdata[2][mm], tempdata[3][mm]);
				}
			}
		}

		fprintf(fp, "\n");
		fclose(fp);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	for (i = 0; i < 2*data->mozaic.numModes; i++) {
		free(tempdata[i]);
		free(im_tempdata[i]);
	}
	free(tempdata);
	free(im_tempdata);

}