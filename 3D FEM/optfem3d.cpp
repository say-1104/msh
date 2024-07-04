#include "fem3d.h"
#include "function.h"
#include "define.h"
#include <sys/time.h>
#include <linux/kernel.h>
#include <sys/sysinfo.h>
#include <string.h>

int main(int argc, char **argv) {
	fprintf(stderr, "Compile date: %s %s\n", __DATE__, __TIME__ );
	long long int i;
	long long int j;
	static DataTable data;
	static long long int ipc;
	static long long int a;
	struct timeval s, e;
	
	char str[256], command[256];
	char namespase[128];
	long long int nn;
	long long int nc;
	long long int *ia;
	long long int *ja;
	long long int count = 0;
	FILE *fp;
	double time1, time2;
	

	long long int nm = -1;
	long long int np = -1;
	long long int mode = 0;
	long long int *test;
	
	double re, im;
	std::complex<double> cj(0.0, 1.0);

	MPI_Init(&argc, &argv);
	fprintf(stderr, "1\n");
	MPI_Comm_size(MPI_COMM_WORLD, &data.par.num_proc);
	fprintf(stderr, "2\n");
	MPI_Comm_rank(MPI_COMM_WORLD, &data.par.myid); 
	fprintf(stderr, "Number of processes : %d\n", data.par.num_proc);

	srand( (unsigned)time( NULL ) ); 
	randomize();

	checkArgument(&data, argc, argv);

	checkArgument1(&data, argc, argv);

	if(data.calcflag == 0){
		fprintf(stderr, "Normal use of 3D FEM\n");
	}else if(data.calcflag == 1){
		fprintf(stderr, "Inputmode calculation only\n");
	}else if(data.calcflag == 2){
		fprintf(stderr, "Inputmode is read by prepared files\n");
	}

	inputData(&data, data.filename);
	detectElementForIntpol(&data);
		

	createPortInfo(&data); 

	if(data.calcflag == 1){
		exit(0);
	}

	// Prepare for Mozaic simulation
	if(data.mozaic.mozflag == 1){
		inputMozaic(&data);
		//checkMeshforeachVoxel(&data);
		//exit(0);
		checkMeshforMozaic(&data);
	}

  	fprintf(stderr, "---------------- before analysis\n");
  
	if(data.mozaic.DBSflag != 0){
		if(data.mozaic.symflag == 1){
			// symmetric in x direction
			fprintf(stderr, "not inplemented\n");
			exit(1);

		} else if(data.mozaic.symflag == 4) {
			// symmetric in x,z direction
			DirectBinarySearch(&data, 0);
		} else {
			// No symmetry
			fprintf(stderr, "not inplemented\n");
			exit(1);
		}
	} else {
		analysis2(&data, 0);
	}
	
	if (data.par.nofield == 1) {
		MPI_Finalize();
		return 0;
	}


  	fprintf(stderr, "intpol start \n");
	  
	sprintf(command, "mkdir -p ./Field/abs-%.3lf", data.par.wavelength);
	system(command);
	sprintf(command, "mkdir -p ./Field/real-%.3lf", data.par.wavelength);
	system(command);
  
    // Structure and Field distribution
	addCoreframe(&data); setRGBTable(&data);
    for (j = 0; j < data.nrhs ; j++){
#pragma omp parallel for
		for(i = 0; i < 3; i++){
			fprintf(stderr, "%lld-th mode %lld / 3\n", j, i);
			intpol(&data, i, j, 0);
		}
    }
  
	MPI_Finalize();
	return 0;

}
