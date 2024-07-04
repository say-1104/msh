#include "optfem.h"

extern int main(int argc, char **argv)
{  
	int i, j, k;
	int ll, mm;
	DataTable data;
	Param	*par = &(data.par);
	Port	*port = &(data.port);
	double wl;
	char string[256];
	FILE *fp;  

	// Check Output directory
	sprintf(string, "rm -r Output");
	system(string);
	sprintf(string, "mkdir Output");
	system(string);
	printf("Hello1");


	checkArgument1(&data, argc, argv);
	inputData(&data, argv[argc-1]);
	checkArgument2(&data, argc, argv);
	fprintf(stderr, "-------------------------\n modeimflag %d\n", data.par.modEIMflag);

	fprintf(stderr, "Number of port = %d\n", data.port.number);

	allocMatrix(&data);

	// preparation for interpolation
	Preforintpol(&data);

	// Prepare for Mosaic simulation
	if (data.mosaic.mosflag == 1) {
		inputMosaic(&data);
		checkMeshforMosaic(&data);
	}

	inputModifiedEIMData(&data);

	srand((unsigned)time(NULL));
	randomize();

	// set input wave
	if (port->solver == POINTSOURCE) {
		Pointsource(&data);
	} 
	else if (port->solver == SHEETSOURCE) {
		Sheetsource(&data);
	} 
	else if (port->solver == INPUTEIGEN) {
		for (ll = 0; ll < data.par.numoflambda; ll++) {
			par->wavelength = par->w_start+(double)ll*par->w_step;
			for (mm = 0; mm < data.par.num_mode; mm++) {
				fprintf(stderr, "In condition wavelength %lf, input mode is TE%d,\n", par->wavelength, mm);
				SetModifiedEI(&data, par->wavelength, ll, mm);
				for (j = 0; j < data.port.number; j++) {
					fprintf(stderr, "Port No.%d  :  ", j);
					portInfo(&(port->dataMemory[ll][mm][j]), &data, j);
				}
				eigenModeOfPort(&data, par->wavelength, ll, mm);
			}
		}
	}
	else {
		fprintf(stderr, "Light source type is not specified.");
		exit(1);
	}

	ALLOCATION(data.mosaic.TT, double *, data.port.number);
	ALLOCATION(data.mosaic.amp, double_complex *, data.port.number);
	ALLOCATION(data.mosaic.arg, double *, data.port.number);

	for (j = 0; j < data.port.number; j++) {
		ALLOCATION(data.mosaic.TT[j], double, data.par.numoflambda*par->num_mode);
		ALLOCATION(data.mosaic.amp[j], double_complex, data.par.numoflambda*par->num_mode);
		ALLOCATION(data.mosaic.arg[j], double, data.par.numoflambda*par->num_mode);
	}

	// Mosaic simulation
	if (data.mosaic.mzcalcNo != -2 && data.mosaic.DBSflag != -1) {
		// DBS Optimization
		DirectBinarySearch(&data);
	}

	// mosaic pattern simulation for GuessedMatrix.csv
	if (data.mosaic.mzcalcNo == -2) fprintf(stderr, "Guessed Matrix\n");
	
	for (ll = 0; ll < data.par.numoflambda; ll++) {
		par->wavelength = par->w_start+(double)ll*par->w_step;
		fprintf(stderr, "wavelength = %lf\n", par->wavelength);

		for (mm = 0; mm < data.par.num_mode; mm++) {  
			SetModifiedEI(&data, par->wavelength, ll, mm);
			if (port->solver == INPUTEIGEN) {
				setPort(&data, par->wavelength, ll, mm);
			}
			
			FEManalysis(&data, par->wavelength, ll);

			if (port->solver != POINTSOURCE) {
				portField(&data, ll);
				probability(&data, par->wavelength, ll, mm);
			}
			
			intpol(&data, par->wavelength, mm);
			
			intpolXX(&data, par->wavelength, mm);
			intpolYY(&data, par->wavelength, mm);
		}
	}
	outputOverlap(&data, 0);

	if (data.mosaic.mzcalcNo != -1 && data.mosaic.mzcalcNo != -2) {
		for(i = 0; i < data.mosaic.Ny; i++) {
			for(j = 0; j < data.mosaic.Nz; j++) {
				printf("%d ", data.mosaic.BorW[i][j]);
			}
			printf("\n");
		}
		printf("\n");
	}

	return 0;
}
