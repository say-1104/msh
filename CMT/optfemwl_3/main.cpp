#include "LI.hpp"
#include "main.hpp"

int Div(double x, double y){ return((int)(round(x*1000000)/round(y*1000000)));
}
void checkArgument(DataTable *data, int argc, char* argv[]);
void checkConfig(DataTable *data);
void inputData(DataTable *data);
void makeZtoW(DataTable *data);

//calcCMT.cpp
void calcCMT(DataTable *data);

int main(int argc, char* argv[]){
    DataTable data;
    //引数の確認
    checkArgument(&data, argc, argv);
    std::cerr << "check argument finish" << std::endl;

    checkConfig(&data);

    inputData(&data);

    makeZtoW(&data);
    std::cerr << "makeZtoW finish" << std::endl;
    calcCMT(&data);
    std::cerr << "CMT finish" << std::endl;
    return 0;
}

void checkArgument(DataTable *data, int argc, char* argv[]){
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
				std::cerr << "Phase Change Material(" << flag->pcm <<") is deposited" << std::endl;
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
                data->par.wl = atof(argv[i]);
                std::cerr << std::fixed << std::setprecision(3);
                std::cerr << "Wavelength is " << data->par.wl << std::endl;
			}
			else if (argv[i][1] == 'l' && argv[i][2] == 'e' && argv[i][3] == 'f' && argv[i][4] == 'f') {
                i++;
                data->par.Leff = atof(argv[i]);
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

void checkConfig(DataTable *data) {
    Param *par = &(data->par);

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
    for(int i=0; i<(par->N_taper)+1; i++){
        double w, z;
        ifs >> z >> w;
        par->taper.push_back(std::make_pair(w, z));
    }
    par->wst = par->taper[0].first;
    if(par->taper[par->N_taper].second < par->Leff) {
        std::cerr << "Leff is incorrect" << std::endl;
        std::cerr <<  par->taper[par->N_taper].second << par->Leff << std::endl;
        exit(1);
    }
    ifs.close();
};

void inputData(DataTable *data){
    Param *par = &(data->par);
    
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3) << data->par.wl;
    std::ifstream ifs(("input_" + oss.str() + ".pre").c_str(), std::ios::out);
    if(! ifs) {
		std::cerr << "input.pre open error !" << std::endl;
		exit(1);
	} else std::cerr << "input.pre open" << std::endl;
    //ALLOCATION(data->dset, Dataset, par->N_dset);

    for(int i=0; i<par->N_dset; i++){
        ifs >> data->dset[i].N;
        for (int j = 0; j < data->dset[i].N; j++) {
            double w;
            double b1, b2, b3, b4;
            ifs >> w >> b1 >> b2 >> b3 >> b4;
            data->dset[i].beta_even.append(w, b1);
            data->dset[i].beta_odd.append(w, b2);
            data->dset[i].beta_1.append(w, b3);
            data->dset[i].beta_2.append(w, b4);
        }
    }
};

void makeZtoW(DataTable *data){
    Param *par = &(data->par);
    int count = 0;
    for (int i=0; i<par->N_taper; i++){
        double l = par->taper[i+1].second - par->taper[i].second;
        int n_div = Div(l, par->dz);
        double wst = par->taper[i].first; double wfi = par->taper[i+1].first;
        for (int step = 0; step < n_div; step++) {
            double cur_z = double(count) * par->dz;
            double stepz = double(step) * par->dz;
            double cur_w = wst + (stepz / l) * (wfi - wst);
            if(double(count)*(par->dz) < par->Leff) par->ZtoW.push_back(std::make_tuple(cur_z, cur_w, CPCM));
            else par->ZtoW.push_back(std::make_tuple(cur_z, cur_w, APCM));
            count++;
        }
    }
    double cur_z = double(count) * par->dz; double cur_w = par->taper[par->N_taper].first;
    if(double(count-1)*(par->dz) < par->Leff) par->ZtoW.push_back(std::make_tuple(cur_z, cur_w, CPCM));
    else par->ZtoW.push_back(std::make_tuple(cur_z, cur_w, APCM));
    std::ofstream ofs;

    ofs.open("ZtoW", std::ios::out);
	if(! ofs) {
		std::cerr << "File(output) open error !" << std::endl;
		exit(1);
    }
    ofs << "z\tw\tpcm\n";
    ofs << std::fixed << std::setprecision(10);
    for (int i=0; i<static_cast<int>(par->ZtoW.size()); i++){
        ofs << std::get<0>(par->ZtoW[i]) << " " << std::get<1>(par->ZtoW[i]) << " " << std::get<2>(par->ZtoW[i]) << std::endl;
    }
    ofs.close();
};

void calcCMT(DataTable *data){
    Param *par = &(data->par);
    std::ofstream ofs;

    ofs.open("output", std::ios::out);
	if(! ofs) {
		std::cerr << "File(output) open error !" << std::endl;
		exit(1);
    }
    ofs << "step\tposition_z\tabs(a)\tabs(b)\n";
    ofs << std::fixed << std::setprecision(10);

    Eigen::Vector2cd ab;
    ab << std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0);
    int n_div = Div(par->taper[par->N_taper].second, par->dz);

    int pcmflag = 0;
    for (int step = 0; step < n_div; step++) {
        auto z = std::get<0>(par->ZtoW[step]);
        auto w = std::get<1>(par->ZtoW[step]);
        auto flag = std::get<2>(par->ZtoW[step]);
        if(flag == 1 && pcmflag == 0) pcmflag++;

        double be = data->dset[pcmflag].beta_even[w];
        double bo = data->dset[pcmflag].beta_odd[w];
        double b1 = data->dset[pcmflag].beta_1[w];
        double b2 = data->dset[pcmflag].beta_2[w];

        auto calc_CMT = [&]() -> void { // ラムダ式
            double beta_ave = (b1 + b2) / 2.0;
            double delta = (b1 - b2) / 2.0;
            double q = (be - bo) / 2.0;
            double kappa = ((q * q - delta * delta >= 0.0) ? sqrt(q * q - delta * delta) : 0.0);

            Eigen::Matrix2d P;
            P << delta + q, kappa, kappa, -(delta + q);
            P *= 1.0 / (sqrt(kappa * kappa + (delta + q) * (delta + q)));

            Eigen::Matrix2cd mid;
            mid << exp(-cj * (beta_ave + q) * par->dz), std::complex<double>(0.0, 0.0),
                std::complex<double>(0.0, 0.0), exp(-cj * (beta_ave - q) * par->dz);

            ab = P * ab;
            ab = mid * ab;
            ab = P * ab;
        };

        ofs << step << '\t' << z << '\t' << abs(ab(0)) * abs(ab(0)) << '\t' << abs(ab(1)) * abs(ab(1)) << std::endl;
        calc_CMT();
    }

    ofs << n_div << '\t' << par->taper[par->N_taper].second << '\t' << abs(ab(0)) * abs(ab(0)) << '\t' << abs(ab(1)) * abs(ab(1)) << std::endl;
    ofs.close();
};