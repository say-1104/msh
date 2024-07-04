#include "fem3d.h"
#include <getopt.h>

void makeTable(DataTable *data);

long long int checkArgument(DataTable *data, long long int argc, char **argv) {
	long long int i;
	for(i = 0; i < 64; i++) data->arguments[i] = 0;
			
	// data->par.div_max = 319;
	data->par.div_max = 959;
	sprintf(data->filename, "prop");
	fprintf(stderr, "size_t = %lld \n", sizeof(size_t) );
	fprintf(stderr, "int = %lld \n", sizeof(int) );
	
	// data->par.inputFlag = NORMAL;
	// data->par.inputFlag = PERIODIC; //unsupported

	while (1) {
		long long int c;
		int this_option_optind = optind ? optind : 1;
		int option_index = 0;
		static struct option long_options[] = {
			{"nm",      required_argument, 0,  0 },
			{"np",      required_argument, 0,  0 },
			{"file",    required_argument, 0,  0 },
			{"detail",  no_argument,       0,  0 },
			{"ooc",     no_argument,       0,  0 },
			{"output",  no_argument,       0,  0 },
			{"mode",    required_argument, 0,  0 },
			{"pict",    required_argument, 0,  0 },
			{"help",    no_argument,       0,  0 },
			{0,         0,                 0,  0 }
		};
		
		c = getopt_long(argc, argv, "abc:d:012", long_options, &option_index);
		if (c == -1) break;
		if (c == 0){
			fprintf(stderr, "option %s", long_options[option_index].name);
			if (optarg) fprintf(stderr, " with arg %s", optarg);
			fprintf(stderr, "\n");
			
			if(option_index == 0) data->nm = atoi(optarg); 
			// 全体行列をつくるときの並列数(並列数分アロケーションするので注意)
			if(option_index == 1) data->np = atoi(optarg); // 上記以降
			if(option_index == 2) strcpy(data->filename, optarg);
			if(option_index == 3) data->arguments[0] = 1;
			if(option_index == 4) data->arguments[1] = 1;
			if(option_index == 5) data->arguments[2] = 1;
			if(option_index == 6) data->mode = atoi(optarg);
			if(option_index == 7) data->par.div_max = atoi(optarg) - 1;
			if(option_index == 8) {fprintf(stderr, "/home/tsato/Solver/3D-Frequency/optfem3d --file [FILENAME] --nm 64 --np 64 --output\n"); exit(1);}
		}
	} // end of while
	if(data->arguments[2] == 1) fprintf(stderr, "[OUTPUT]\n");
	if(data->mode == 3)         fprintf(stderr, ":::READ mode:::");
	fprintf(stderr, "intpol: %lld px\n", data->par.div_max + 1);
	return 0;
}

extern int checkArgument1(DataTable *data, long long int argc, char **argv)
{
	int i, j;
	
	data->calcflag = 0;
	data->mozaic.DBSflag = 0;
	data->mozaic.mozDev = 0;
	data->mozaic.DBSitr = 10;
	data->mozaic.numModes = 1;
	data->par.numStructures = 1;
	data->par.nofield = 0;
	data->pict.modifiedScale = false;
	data->pict.originalSize = false;


	// data->mozaic.mozDev
	// 1: 2 mode divider
	// 2: 2 mode MUX
	// 3: mode EX 
	// 4: mode separator 
	// 5: waveguide lense
	// 6: waveguide crossing
	// 7: nxn Power Splitter
	// for 3 to 5, the number of modes is designated to data->mozaic.numModes

	for (i = 0; i  < argc; i++) {
		if (argv[i][0] == '-') {
			if (argv[i][1] == 'm' && argv[i][2] == 'o' 
				&& argv[i][3] == 'd' && argv[i][4] == 'e') {
				sscanf(argv[i]+5, "%d", &data->calcflag);
				fprintf(stderr, "calcflag = %d\n", data->calcflag);
			} else if (argv[i][1] == 'm' && argv[i][2] == 'z') {	
				data->mozaic.mozflag = 1;
				fprintf(stderr, "Mozaic pattern simulation\n");
			} else if (argv[i][1] == 'f' && argv[i][2] == 'm' && argv[i][3] == 'z') {
				sscanf(argv[i]+4, "%d", &data->mozaic.mzcalcNo);
				fprintf(stderr, "%d th mozaic structure is calculated\n", 
					data->mozaic.mzcalcNo);
			} else if (argv[i][1] == 'D' && argv[i][2] == 'B' && argv[i][3] == 'S') {
				sscanf(argv[i]+4, "%d", &data->mozaic.DBSflag);
				fprintf(stderr, "Direct Binary Serch optimization");
				if (data->mozaic.DBSflag == 1){
					fprintf(stderr, " from beggining\n");
				}else if (data->mozaic.DBSflag == 2){
					fprintf(stderr, " continued\n");
				}
			} else if (argv[i][1] == 'D' && argv[i][2] == 'I' && argv[i][3] == 'V') {
				fprintf(stderr, "modeDivider is optimized\n");
				data->mozaic.mozDev = 1;
			} else if (argv[i][1] == 'M' && argv[i][2] == 'U' && argv[i][3] == 'X') {
				fprintf(stderr, "modeMUX is optimized\n");
				data->mozaic.mozDev = 2;
			} else if (argv[i][1] == 'E' && argv[i][2] == 'X') {
				sscanf(argv[i]+3, "%d", &data->mozaic.numModes);
				fprintf(stderr, "%d modeEX is optimized\n", data->mozaic.numModes);
				data->mozaic.mozDev = 3;
			} else if (argv[i][1] == 'S' && argv[i][2] == 'E' && argv[i][3] == 'P') {
				sscanf(argv[i]+4, "%d", &data->mozaic.numModes);
				fprintf(stderr, "%d mode separator is optimized\n", 
					data->mozaic.numModes);
				data->mozaic.mozDev = 4;
			} else if (argv[i][1] == 'W' && argv[i][2] == 'L') {
				sscanf(argv[i]+3, "%d", &data->mozaic.numModes);
				fprintf(stderr, "%d mode waveguide lense is optimized\n", 
					data->mozaic.numModes);
				data->mozaic.mozDev = 5;
			} else if (argv[i][1] == 'W' && argv[i][2] == 'C') {
				sscanf(argv[i]+3, "%d", &data->mozaic.numModes);
				fprintf(stderr, "%d mode waveguide crossing is optimized\n", 
					data->mozaic.numModes);
				data->mozaic.mozDev = 6;
			} else if (argv[i][1] == 'P' && argv[i][2] == 'S') {
				sscanf(argv[i]+3, "%d", &data->mozaic.numModes);
				fprintf(stderr, "%d x %d waveguide crossing is optimized\n", 
					data->mozaic.numModes, data->mozaic.numModes);
				data->mozaic.mozDev = 7;
			} else if (argv[i][1] == 'f' && argv[i][2] == 't') {
				sscanf(argv[i]+3, "%d", &data->par.numStructures);
				fprintf(stderr, "%d strucutres are considered to improve fabrication tolerance\n", 
					data->par.numStructures);
			} else if (argv[i][1] == 'F' && argv[i][2] == 'T') {
				sscanf(argv[i]+3, "%d", &data->par.numStructures);
				data->par.numStructures *= -1;
				fprintf(stderr, "%d strucutres are considered to improve fabrication tolerance\n", 
					data->par.numStructures);
			} else if (argv[i][1] == 's') {
				data->par.nofield = 1;
				fprintf(stderr, "No field Output\n");
			} else if (argv[i][1] == 'S' && argv[i][2] == 'C' && argv[i][3] == 'A' && argv[i][4] == 'L' && argv[i][5] == 'E') {
				data->pict.modifiedScale = true;
				fprintf(stderr, "Scaling of the Field is modified.\n");
			} else if (argv[i][1] == 'O' && argv[i][2] == 'R' && argv[i][3] == 'I' && argv[i][4] == 'G' && argv[i][5] == 'I' && argv[i][6] == 'N' && argv[i][7] == 'A' && argv[i][8] == 'L') {
				data->pict.originalSize = true;
				fprintf(stderr, "Original size is used for field outputs.\n");
			}
		} 
	}
	fprintf(stdout, "Number of processes is %d\n", data->par.num_proc);
}

void inputData(DataTable *data, char *name) {
	long long int i, j, k;
	char file[256];
	FILE *fp;
	
	FEM *fem = &(data->fem);
	Param *par = &(data->par);
	Picture *pict = &(data->pict);

	fem->n_node = 10;
	par->mode.target = 1;
	fem->unknown = 1;

	if (data->par.numStructures < 0) {
		// 異なる構造については、中心波長しか考慮しない
		data->par.numStructures *= -1;
		data->par.numofwlMPI = (par->num_proc - (data->par.numStructures - 1) );
		int c_str = data->par.numStructures / 2, c_wl = data->par.numofwlMPI / 2;

		if (data->par.myid < data->par.numofwlMPI) {
			data->par.wlID = data->par.myid % data->par.numofwlMPI;
			data->par.strID = c_str;
		} else {
			data->par.wlID = c_wl;
			data->par.strID = (data->par.myid - data->par.numofwlMPI);
			if (data->par.strID >= c_str) data->par.strID++;
		}
	} else {
		if (par->num_proc % data->par.numStructures != 0) {
			fprintf(stderr, "指定プロセス数は、構造数の整数倍で指定してください。\n");
			fprintf(stderr, "指定プロセス数:%d, 構造数:%d\n", par->num_proc, data->par.numStructures);
			exit(1);
		}
		data->par.wlID = data->par.myid/data->par.numStructures;
		data->par.strID = data->par.myid%data->par.numStructures;
		data->par.numofwlMPI = par->num_proc/data->par.numStructures;
	}
	/* ----------- param read ---------------- */
	if (data->par.numofwlMPI == 1) sprintf(file, "%s.param", name);
	else sprintf(file, "%s-%d.param", name, data->par.wlID);
	if ((fp = fopen(file, "r")) == NULL) {
		fprintf(stderr, "can't open file ( %s )\n", file);
		exit(EXIT_FAILURE);
	}
	
	fscanf(fp, "%*s %lld", &fem->unknown);
	fscanf(fp, "%*s %lf", &par->wavelength);

	MPI_Barrier(MPI_COMM_WORLD);
	fprintf(stderr, "myid = %d (str:%d wl:%d), wavelength = %lf\n", 
		data->par.myid, data->par.strID, data->par.wlID, par->wavelength);

	// -------------------------------------

	fscanf(fp, "%*s %lf", &data->pml.tanD);
	fscanf(fp, "%*s %lf", &data->pml.m);
	
	fscanf(fp, "%*s %lf %lf %lf %lf %lf %lf",
	 &data->pml.dx1, &data->pml.dx2,
	 &data->pml.dy1, &data->pml.dy2,
	 &data->pml.dz1, &data->pml.dz2);
	
	fscanf(fp, "%*s %lf", &par->INPUT_Z);
	fscanf(fp, "%*s %lf", &par->neff);
	fscanf(fp, "%*s %lf", &par->neff_out);
	
	fscanf(fp, "%*s %lld", &par->n_material);
	for (i = 0; i < par->n_material; i++) {
		fscanf(fp, "%lf %lf %lld %lf %lf %lf", &(par->er[i]), &(par->mr[i]),
		 &(par->PMLflag[i]), &(par->PMLangle[i]),
		 &(par->PMLz0[i]), &(par->PMLx0[i]));
		par->er[i] = par->er[i] * par->er[i];
		par->mr[i] = par->mr[i] * par->mr[i];
	}
	
	for (k = 0; k < data->par.n_material; k++) {
		for (i = 0; i < 3; i++) {
			for (j = 0; j < 3; j++) {
				if (i == j) {
					data->par.ermatrix[k][i][j] = data->par.er[k];
					data->par.mrmatrix[k][i][j] = data->par.mr[k];
				} else {
					data->par.ermatrix[k][i][j] = 0;
					data->par.mrmatrix[k][i][j] = 0;
				}
			}
		}
	}

	fscanf(fp, "%*s %lf", &(par->power_adjust));
	fscanf(fp, "%*s %s", par->in_file);
	fscanf(fp, "%*s %lf", &par->out_slv_y_cood);
	fscanf(fp, "%*s %s", data->pml_setting_2D_input);
	fscanf(fp, "%*s %s", data->pml_setting_2D_output);
	fscanf(fp, "%*s %lf", &data->pml_setting_2D_thickness);
	
	//Port extension
	fscanf(fp, "%*s %lld", &data->n_port_io);
	for (i = 0; i < data->n_port_io; i++) {
		fscanf(fp, "%lf %lf %lf %lf", 
		 &(data->port_xmin[i]), &(data->port_xmax[i]),
		 &(data->port_zmin[i]), &(data->port_zmax[i])
		 );
	}
	fscanf(fp, "%*s %lf", &par->out_slv_x_cood);
	fscanf(fp, "%*s %lf", &par->out_slv_z_cood);
	
	// for mozaic?, sato program
	fscanf(fp, "%*s %lf %lf %lf %lf %lf %lf",
	 &par->x1, &par->x2,
	 &par->y1, &par->y2,
	 &par->z1, &par->z2);
	fscanf(fp, "%*s %lld %lld", &par->I, &par->J);
	
	ALLOCATION(par->mat, long long int, par->I * par->J);
	for (i = 0; i < par->I; i++) {
		for (j = 0; j < par->J; j++) {
			fscanf(fp, "%lld", &(par->mat[par->J*i + j]));
			par->mat[par->J*i + j]--;
		}
	}
	
	// For input field from files, added by Fujisawa
	if(data->calcflag == 2){
		fscanf(fp, "%*s %lld", &data->fem.n_port);
		ALLOCATION(data->port, PortInfo, fem->n_port);
		for (i = 0; i < data->fem.n_port; i++) {
			fscanf(fp, "%lld",  &data->port[i].modenum);
		}
	}

	if (fclose(fp) == EOF) {
		fprintf(stderr, "can't close file ( %s )\n", file);
		exit(EXIT_FAILURE);
	}
	
	if (1) {
		fprintf(stderr, "-------- calculation conditions --------\n");
		
		fprintf(stderr, "unknown = %lld\n", fem->unknown);
		fprintf(stderr, "wavelength = %lf\n", par->wavelength);
		fprintf(stderr, "PML : tanD = %lf, m = %lf\n", data->pml.tanD, data->pml.m);
		fprintf(stderr,
			"pml thickness :\n  x_min = %lf; x_max = %lf\n  y_min = %lf; y_max = %lf\n  z_min = %lf; z_max = %lf\n",
			data->pml.dx1, data->pml.dx2,
			data->pml.dy1, data->pml.dy2,
			data->pml.dz1, data->pml.dz2);
		for (i = 0; i < par->n_material; i++) {
			fprintf(stderr,
				"mat[%2lld]: er, mr = %6.3lf, %.1lf, flag = %lld, ", 
				i, par->er[i], par->mr[i], par->PMLflag[i]);
			fprintf(stderr,
				"angle = %3.0lf, z0, x0 = %6.3lf, %6.3lf\n", 
				par->PMLangle[i], par->PMLz0[i], par->PMLx0[i]);
		}
		fprintf(stderr, "Power adjust  : %lf\n", par->power_adjust);
		fprintf(stderr, "Input         : %s(disabled)\n", par->in_file);
		fprintf(stderr, "Outslv x cood : %lf\n", par->out_slv_x_cood);
		fprintf(stderr, "Outslv y cood : %lf\n", par->out_slv_y_cood);
		fprintf(stderr, "Outslv z cood : %lf\n", par->out_slv_z_cood);
		fprintf(stderr, "----------------------------------------\n");
	}


	ALLOCATION(data->name, char, (strlen(name)+1));
	strcpy(data->name, name);
	
	/* --------------- msh read --------------- */
	if (data->par.numStructures != 1) sprintf(file, "%s-%d.msh", name, data->par.strID);
	else sprintf(file, "%s.msh", name);
	if ((fp = fopen(file, "r")) == NULL) {
		fprintf(stderr, "can't open file ( %s )\n", file);
		exit(EXIT_FAILURE);
	}
	
	par->k0 = 2.0 * M_PI / par->wavelength;
	par->k02 = par->k0 * par->k0;
	
	fscanf(fp, "%lld", &fem->ne);
	fscanf(fp, "%lld", &fem->np);
	ALLOCATION(fem->element, Element3D, fem->ne);
	for (i = 0; i < fem->ne; i++) {
		fscanf(fp, "%lld", &(fem->element[i].matID));
		fem->element[i].matID--;
		for (j = 0; j < fem->n_node; j++) {
			fscanf(fp, "%lld", &(fem->element[i].kn[j]));
			fem->element[i].kn[j]--;
		}
	}
	
	ALLOCATION(fem->node, XYZdouble, fem->np);
	fem->min.x = fem->min.y = fem->min.z = 1.0e50;
	fem->max.x = fem->max.y = fem->max.z = -1.0e50;
	fprintf(stderr, "np = %lld \n",fem->np);
	for (i = 0; i < fem->np; i++) {
		fscanf(fp, "%lf %lf %lf", 
		 &(fem->node[i].x), &(fem->node[i].y), &(fem->node[i].z));
		if (fem->node[i].x < fem->min.x) fem->min.x = fem->node[i].x;
		if (fem->node[i].y < fem->min.y) fem->min.y = fem->node[i].y;
		if (fem->node[i].z < fem->min.z) fem->min.z = fem->node[i].z;
		if (fem->node[i].x > fem->max.x) fem->max.x = fem->node[i].x;
		if (fem->node[i].y > fem->max.y) fem->max.y = fem->node[i].y;
		if (fem->node[i].z > fem->max.z) fem->max.z = fem->node[i].z;
	}

	fscanf(fp, "%lld", &fem->order);
	switch (fem->order) {
	case 1:
		fem->n_en = 6;
		fem->n_node_use = 4;
		break;
	case 2:
		fem->n_en = 12;
		fem->n_node_use = 10;
		break;
	case 3:
		fem->n_en = 24;
		fem->n_node_use = 10;
		break;
	}
	fscanf(fp, "%lld", &fem->n_edge);
	fscanf(fp, "%lld", &fem->nr_edge);
	
	fprintf(stderr, "Input msh file : %s\n", file);
	fprintf(stderr, "k0 = %lf\n", par->k0);
	fprintf(stderr, "whole region :\n");
	fprintf(stderr, "  (%6.3lf + %6.3lf < x < %6.3lf - %6.3lf)\n", 
		fem->min.x, data->pml.dx1, fem->max.x, data->pml.dx2);
	fprintf(stderr, "  (%6.3lf + %6.3lf < y < %6.3lf - %6.3lf)\n", 
		fem->min.y, data->pml.dy1, fem->max.y, data->pml.dy2);
	fprintf(stderr, "  (%6.3lf + %6.3lf < z < %6.3lf - %6.3lf)\n", 
		fem->min.z, data->pml.dz1, fem->max.z, data->pml.dz2);
	fprintf(stderr, "total unkowns(nr_edge) = %lld\n", fem->nr_edge);
	fprintf(stderr, "n_edge = %lld\n", fem->n_edge);

	pict->nx = (long long int)((fem->max.x - fem->min.x)*PICTRATIO);
	pict->ny = (long long int)((fem->max.y - fem->min.y)*PICTRATIO);
	pict->nz = (long long int)((fem->max.z - fem->min.z)*PICTRATIO);
	
	if (pict->nz > 10000 && !pict->originalSize) {
		double pictratio = 10;
		pict->nx = (long long int)((fem->max.x - fem->min.x)*pictratio);
		pict->ny = (long long int)((fem->max.y - fem->min.y)*pictratio);
		pict->nz = (long long int)((fem->max.z - fem->min.z)*pictratio);
	}
	fprintf(stderr, "Number of pixels for (x,y,z) direction : (%lld, %lld, %lld)\n", pict->nx, pict->ny, pict->nz);
	/* ---------------------------------------- */

	ALLOCATION(fem->inFlag, long long int, fem->n_edge);
	for (i = 0; i < fem->ne; i++) {
		for (j = 0; j < fem->n_en; j++) {
			fscanf(fp, "%lld", &(fem->element[i].kk[j]));
			fem->element[i].kk[j]--;
		}
	}
	
	if (fclose(fp) == EOF) {
		fprintf(stderr, "can't close file ( %s )\n", file);
		exit(EXIT_FAILURE);
	}
	
	makeTable(data);
	createS3d(data);
	
	return;
}

//added by satonori
void makeTable(DataTable *data){
	fprintf(stderr, "makeTable();\n");
	long long int i, j, k, itmp;
	long long int count;
	double x;
	double y;
	double z;
	FEM *fem = &(data->fem);
	Param *par = &(data->par);
	
	fem->table_n = (long long int**)malloc(sizeof(long long int*) * par->I);
	for (i = 0; i < par->I; i++){
		fem->table_n[i] = (long long int*)malloc(sizeof(long long int) * par->J);
		for (j = 0; j < par->J; j++){
			fem->table_n[i][j] = 0;
			fprintf(stderr, "%lld ", par->mat[par->J*i + j]);
		}
		fprintf(stderr, "\n");
	}
	
	count = 0;
	for (k = 0; k < fem->ne; k++) {
		x = y = z = 0.0;
		for (i = 0; i < 4; i++) {
			itmp = fem->element[k].kn[i];
			x += fem->node[itmp].x;
			y += fem->node[itmp].y;
			z += fem->node[itmp].z;
		}
		x *= 0.250; y *= 0.250; z *= 0.250;
		if(x < par->x1) continue; 
		if(y < par->y1) continue; 
		if(z < par->z1) continue;
		if(x > par->x2) continue; 
		if(y > par->y2) continue; 
		if(z > par->z2) continue;
				
		j = (z - par->z1) / (par->z2 - par->z1) * par->J;
		i = (x - par->x1) / (par->x2 - par->x1) * par->I;
		
		fem->element[k].matID = par->mat[par->J*i + j];
		fem->table_n[i][j]++;
		count++;
	}
	
	fem->table = (long long int***)malloc(sizeof(long long int**) * par->I);
	for (i = 0; i < par->I; i++){
		fem->table[i] = (long long int**)malloc(sizeof(long long int*) * par->J);
		for (j = 0; j < par->J; j++){
			fem->table[i][j] = (long long int*)malloc(sizeof(long long int) * fem->table_n[i][j]);
			fem->table_n[i][j] = 0;
		}
	}
	
	for (k = 0; k < fem->ne; k++) {
		x = y = z = 0.0;
		for (i = 0; i < 4; i++) {
			itmp = fem->element[k].kn[i];
			x += fem->node[itmp].x;
			y += fem->node[itmp].y;
			z += fem->node[itmp].z;
		}
		x *= 0.250; y *= 0.250; z *= 0.250;
		if(x < par->x1) continue; if(y < par->y1) continue; if(z < par->z1) continue;
		if(x > par->x2) continue; if(y > par->y2) continue; if(z > par->z2) continue;
				
		j = (z - par->z1) / (par->z2 - par->z1) * par->J;
		i = (x - par->x1) / (par->x2 - par->x1) * par->I;
		
		fem->table[i][j][fem->table_n[i][j]] = k;
		fem->table_n[i][j] += 1;
	}
		
	fprintf(stderr, "%6.3lf < x < %6.3lf\n", par->x1, par->x2);
	fprintf(stderr, "%6.3lf < y < %6.3lf\n", par->y1, par->y2);
	fprintf(stderr, "%6.3lf < z < %6.3lf\n", par->z1, par->z2);
	fprintf(stderr, "table count = %lld\n", count);
		
	double x1, x2;
	double y1, y2;
	double z1, z2;
	double l, lmin, lmax;
	lmin = 1.0e10;
	lmax = 0.0;
	long long int itmp2;
	for (k = 0; k < fem->ne; k++) {
		for (i = 0; i < 4; i++) {
			itmp = fem->element[k].kn[i];
			itmp2 = fem->element[k].kn[(i+1)%4];
			x1 = fem->node[itmp].x;
			y1 = fem->node[itmp].y;
			z1 = fem->node[itmp].z;
			x2 = fem->node[itmp2].x;
			y2 = fem->node[itmp2].y;
			z2 = fem->node[itmp2].z;
			x = x1-x2;
			y = y1-y2;
			z = z1-z2;
			l = sqrt(x*x + y*y + z*z);
			if (lmin > l) lmin = l;
			if (lmax < l) lmax = l;
		}
	}
	fprintf(stderr, "lmin = %.0lf nm\n", lmin*1000.0);
	fprintf(stderr, "lmax = %.0lf nm\n", lmax*1000.0);
	return;
}
