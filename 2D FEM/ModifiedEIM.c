#include "optfem.h"

extern int inputModifiedEIMData(DataTable *data) {
    if (data->par.modEIMflag == 0) return 0;

    double wl, n_core, n_clad;
    char filename[256];
    FILE *fp;
    Param* par = &(data->par);

    ModifiedEffectiveIndex *ModEI;
    ModEI = &(data->ModEI);

	std::vector<double> ers;
    for (int i = 0; i < par->n_material; i++) {
        if (ers.size() == 0) { 
			ers.push_back(par->er[i]);
        }
		else {
			std::vector<double>::iterator itr;
			itr = std::find(ers.begin(), ers.end(), par->er[i]);
			if (itr == ers.end()) {
				ers.push_back(par->er[i]);
			}
		}
    }

	fprintf(stderr, "材料は全部で%d個あります。\n", ers.size());
    sprintf(filename, "index.pre");
    if ((fp = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "can't open file ( %s )\n", filename);
        exit(1);
    }

	if (data->par.modEIMflag == 1) {
		ModEI->size = 2;
		ALLOCATION(ModEI->er, std::vector<double>*, par->numoflambda);
		for (int ll = 0; ll < par->numoflambda; ll++) {
			ALLOCATION(ModEI->er[ll], std::vector<double>, par->num_mode);
		}

		fprintf(stderr, "Reading ModifiedEIM Data...\n");
		fscanf(fp, "%*s %*s %*s %*s %*s %*s %*s");
		for (int ll = 0; ll < par->numoflambda; ll++) {
			fscanf(fp, "%*s %lf", &wl);

			if (fabs(par->w_start+(double)ll*par->w_step - wl) > 1e-6) {
				fprintf(stderr, "%s内の波長設定が間違っております。\n", filename);
				fprintf(stderr, "wavelength %lf は、正しくは%lfで記載しなくてはなりません。\n", wl, par->w_start+(double)ll*par->w_step);
				fprintf(stderr, "diff: %.9f", fabs(wl - par->w_start+(double)ll*par->w_step));
				exit(1);
			}

			fprintf(stderr,"wavelength: %lf\n", wl);
			for (int mm = 0; mm < par->num_mode; mm++) {
				ModEI->er[ll][mm].resize(2);
				fscanf(fp, "%lf %lf", &(ModEI->er[ll][mm][0]), &(ModEI->er[ll][mm][1]));
				ModEI->er[ll][mm][0] *= ModEI->er[ll][mm][0];
				ModEI->er[ll][mm][1] *= ModEI->er[ll][mm][1];
				std::sort(ModEI->er[ll][mm].begin(), ModEI->er[ll][mm].end(), std::greater<double>());
				
				fprintf(stderr,"TE%d Mode -> core EI:%1.3lf, clad EI:%1.3lf\n", mm, sqrt(ModEI->er[ll][mm][0]), sqrt(ModEI->er[ll][mm][1]));
			}
		}
		
	} else if (data->par.modEIMflag == 2) {
		ModEI->size = ers.size();
		ALLOCATION(ModEI->er, std::vector<double>*, par->numoflambda);

		fprintf(stderr, "Reading ModifiedEIM Data...\n");
		fscanf(fp, "%*s %*s %*s %*s %*s %*s %*s");
		for (int ll = 0; ll < par->numoflambda; ll++) {
			ALLOCATION(ModEI->er[ll], std::vector<double>, par->num_mode);
			fscanf(fp, "%*s %lf", &wl);

			if (fabs(par->w_start+(double)ll*par->w_step - wl) > 1e-6) {
				fprintf(stderr, "%s内の波長設定が間違っております。\n", filename);
				fprintf(stderr, "wavelength %lf は、正しくは%lfで記載しなくてはなりません。\n", wl, par->w_start+(double)ll*par->w_step);
				fprintf(stderr, "diff: %.9f", fabs(wl - par->w_start+(double)ll*par->w_step));
				exit(1);
			}

			fprintf(stderr,"wavelength: %lf\nIndex: ", wl);
			ModEI->er[ll][0].resize(ers.size());
			for (int i=0; i<ers.size(); i++) {
				fscanf(fp, "%lf", &(ModEI->er[ll][0][i]));
				ModEI->er[ll][0][i] *= ModEI->er[ll][0][i];
			}

			std::sort(ModEI->er[ll][0].begin(), ModEI->er[ll][0].end(), std::greater<double>());
			for (int i=0; i<ers.size(); i++) {
				fprintf(stderr, "%d: %1.3lf, ", i, sqrt(ModEI->er[ll][0][i]));
			}	
			
			for (int mm = 0; mm < par->num_mode; mm++) {
				ModEI->er[ll][mm].resize(ers.size());
				for (int i=0; i<ers.size(); i++) {
					ModEI->er[ll][mm][i] = ModEI->er[ll][0][i];
				}
			}
		}
	} else {
		fprintf(stderr, "サポートされていないmodEIMflagです。\n");
		exit(EXIT_FAILURE);
	}
}

extern int SetModifiedEI(DataTable *data, double wl, int ll, int mm) {
    if (data->par.modEIMflag == 0) return 0;

    Param* par = &(data->par);
    ModifiedEffectiveIndex *ModEI;
    ModEI = &(data->ModEI);

	std::vector<double> ers;
    for (int i = 0; i < par->n_material; i++) {
        if (ers.size() == 0) { 
			ers.push_back(par->er[i]);
        }
		else {
			std::vector<double>::iterator itr;
			itr = std::find(ers.begin(), ers.end(), par->er[i]);
			if (itr == ers.end()) {
				ers.push_back(par->er[i]);
			}
		}
    }
	std::sort(ers.begin(), ers.end(), std::greater<double>());
	if (data->par.modEIMflag == 1 && ers.size() >= 3) {
		fprintf(stderr, "error at SetModifiedEI. three materials are found\n");
		exit(1);
	}

    for (int i = 0; i < par->n_material; i++) {
		std::vector<double>::iterator itr;
		itr = std::find(ers.begin(), ers.end(), par->er[i]);
		if (itr == ers.end()) {
			fprintf(stderr, "error at SetModfiedEI. index table is wrong.");
			exit(1);
		}
		const int idx = std::distance(ers.begin(), itr);
		par->er[i] = ModEI->er[ll][mm][idx];
    }

    for (int i = 0; i < par->n_material; i++) {
        par->epT[i].xx = par->epT[i].yy = par->epT[i].zz = par->er[i];
        par->epT[i].xy = par->epT[i].xz = 0.0;
        par->epT[i].yx = par->epT[i].yz = 0.0;
        par->epT[i].zx = par->epT[i].zy = 0.0;
    }

    fprintf(stderr, "Refractive Index has been renewed as follows by ModEIM,\n");
    for (int i = 0; i < par->n_material; i++) {
        fprintf(stderr, "Mat%d: %1.3lf\n", i+1, sqrt(par->er[i]));
    }
}
