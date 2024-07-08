#include "main.hpp"

extern void checkArgument(DataTable *data, int argc, char* argv[]){
    Param *par = &(data->par);
    Flag *flag = &(data->flag);

    flag->pcm = 0;
    par->N_dset = 1;
    par->wl = 1.550;
    par->Leff = 0.0;

    for (int i = 0; i < argc; i++) {
		if (argv[i][0] == '-') {
			if (argv[i][1] == 'p' && argv[i][2] == 'c' && argv[i][3] == 'm') {
                char c = argv[i][4];
                flag->pcm = int(c - '0');
				std::cerr << "Phase Change Material is deposited" << std::endl;
                switch(flag->pcm){
                    case 1:
                        par->N_dset = 2;
                        break;
                    case 2:
                        par->N_dset = 4;
                        break;
                    default:
                        std::cerr << "Argument(-pcm) is incorrect" << std::endl;
                        exit(1);
                        break;
                }
			}
			else if (argv[i][1] == 'w' && argv[i][2] == 'l') {
                i++;
                data->par.wl = atoi(argv[i]);
                std::cerr << std::fixed << std::setprecision(3);
                std::cerr << "Wavelength is " << data->par.wl << std::endl;
			}
			else if (argv[i][1] == 'l' && argv[i][2] == 'e' && argv[i][3] == 'f' && argv[i][4] == 'f') {
                i++;
                data->par.Leff = atoi(argv[i]);
                std::cerr << std::fixed << std::setprecision(3);
                std::cerr << "Leff is " << data->par.Leff << std::endl;
			}
            else {
                std::cerr << "Argument is incorrect" << std::endl;
                exit(1);
            }
		} else {

        }
	}
};

extern void checkConfig(DataTable *data) {
    Param *par = &(data->par);
    Flag *flag = &(data->flag);

    std::ifstream ifs("structure.cfg", std::ios::out);
    if(! ifs) {
		std::cerr << "structure.cfg open error !" << std::endl;
		exit(1);
	}

    std::string buff;

    ifs >> buff >> par->dz;
    ifs >> buff >> par->N_taper;
    ifs >> buff >> buff;
    //bool b_tmp = true;
    for(int i=0; i<par->N_taper+1; i++){
        double w, z;
        ifs >> w >> z;
        /*if(z > par->Leff && b_tmp) {
            double wleff = (par->Leff - par->taper[i-1].second) * (w - par->taper[i-1].first) / (z - par->taper[i-1].second) + par->taper[i-1].first;
            par->taper.push_back(std::make_pair(wleff, par->Leff));
        } */
        par->taper.push_back(std::make_pair(w, z));
    }
    par->wst = par->taper[0].first;
    if(par->taper[par->N_taper].second < par->Leff) {
        std::cerr << "Leff is incorect" << std::endl;
        exit(1);
    }
    ifs.close();
};

extern void inputData(DataTable *data){
    Param *par = &(data->par);
    Flag *flag = &(data->flag);
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3) << data->par.wl;
    std::ifstream ifs(("input_" + oss.str() + ".pre").c_str(), std::ios::out);
    if(! ifs) {
		std::cerr << "input.pre open error !" << std::endl;
		exit(1);
	}
    ALLOCATION(data->dset, Dataset, par->N_dset);

    for(int i=0; i<par->N_dset; i++){
        ifs >> data->dset->N;
        for (int _ = 0; _ < data->dset->N; _++) {
            double w;
            double b1, b2, b3, b4;
            ifs >> w >> b1 >> b2 >> b3 >> b4;
            data->dset->beta_even.append(w, b1);
            data->dset->beta_odd.append(w, b2);
            data->dset->beta_1.append(w, b3);
            data->dset->beta_2.append(w, b4);
        }
        data->dset++;
    }
};

extern void makeZtoW(DataTable *data){
    Param *par = &(data->par);
    int count = 0;
    for (int i=0; i<par->N_taper; i++){
        double l = par->taper[i+1].second - par->taper[i].second;
        int n_div = l / par->dz;
        double wst = par->taper[i].first, wfi = par->taper[i].first;
        for (int step = 0; step < n_div; step++) {
            double cur_z = step * par->dz;
            double cur_w = wst + (cur_z / l) * (wfi - wst);
            if(count*par->dz < par->Leff) par->ZtoW[step] = std::make_tuple(cur_z, cur_w, CPCM);
            else par->ZtoW[step] = std::make_tuple(cur_z, cur_w, APCM);
            count++;
        }
    }
};