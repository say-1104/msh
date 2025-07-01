#include "LI.hpp"
#include "main.hpp"

int Div(double x, double y){ return((int)(round(x*100000)/round(y*100000)));}

void checkArgument(DataTable *data, int argc, char* argv[]);
void checkConfig(DataTable *data);
void inputData(DataTable *data);
void makeZtoW(DataTable *data);
void calcCMT(DataTable *data);

int main(int argc, char* argv[]){
    DataTable data;
    if(IsFileExist("ZtoW")) system("rm ZtoW");

    // 引数の確認
    checkArgument(&data, argc, argv);
    std::cerr << "check argument finish" << std::endl;

    // structure.cfgの読み込み
    checkConfig(&data);
    std::cerr << "check config finish" << std::endl;

    // input_[波長].preの読み込み
    inputData(&data);
    std::cerr << "check inputdata finish" << std::endl;

    // ZtoWテーブルの作成
    makeZtoW(&data);
    std::cerr << "makeZtoW finish" << std::endl;

    // 行列計算
    calcCMT(&data);
    std::cerr << "CMT finish" << std::endl;
    return 0;
}

void checkArgument(DataTable *data, int argc, char* argv[]){
    Param *par = &(data->par);
    Flag *flag = &(data->flag);

    flag->zw = false;
    flag->sbend = false;

    par->N_dset = 1;
    par->wl = 1.550;

    for (int i = 0; i < argc; i++) {
		if (argv[i][0] == '-') {
			if (argv[i][1] == 'p' && argv[i][2] == 'c' && argv[i][3] == 'm') {
                char c = argv[i][4];
                par->N_dset = int(c - '0');
				std::cerr << "Phase Change Material(" << par->N_dset <<") is deposited" << std::endl;
			}
			else if (argv[i][1] == 'w' && argv[i][2] == 'l') {
                i++;
                data->par.wl = atof(argv[i]);
                std::cerr << std::fixed << std::setprecision(3);
                std::cerr << "Wavelength is " << data->par.wl << std::endl;
			}
            else if (argv[i][1] == 'z' && argv[i][2] == 'w') {
                flag->zw = true;
                std::cerr << "Output ZtoW file " << std::endl;
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
    Flag *flag = &(data->flag);

    std::ifstream ifs("structure.cfg", std::ios::out);
    if(! ifs) {
		std::cerr << "structure.cfg open error !" << std::endl;
		exit(1);
	}

    std::string buff;
    ifs >> buff >> par->dz;
    ifs >> buff >> par->gst;
    ifs >> buff >> par->w;
    ifs >> buff >> par->N_elem;

    double tmp_g = par->gst;
    for(int i=0; i<(par->N_elem); i++){
        std::string elem_name;
        ifs >> elem_name;
        if(elem_name == "sbend"){
            flag->sbend = true;
            double L, gfi;
            int pcm;
            ifs >> buff >> L;
            ifs >> buff >> gfi;
            ifs >> buff >> pcm;
            Element elem = {elem_name, L, gfi, pcm};
            par->elem.push_back(elem);
            tmp_g = gfi;
        }else if(elem_name == "dc"){
            double L;
            double gfi = tmp_g;
            int pcm;
            ifs >> buff >> L;
            ifs >> buff >> pcm;
            Element elem = {elem_name, L, gfi, pcm};
            par->elem.push_back(elem);
            tmp_g = gfi;
        }
    }

    ifs.close();
    /*for (int i=0; i<par->N_elem; i++){
        std::cout << par->elem[i].elem_name << std::endl;
        std::cout << par->elem[i].gfi << std::endl;
        std::cout << par->elem[i].L << std::endl;
        std::cout << par->elem[i].pcm << std::endl;
    }*/
};

void inputData(DataTable *data){
    Param *par = &(data->par);
    Flag *flag = &(data->flag);

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3) << data->par.wl;
    std::ifstream ifs(("inputdata/input_" + oss.str() + ".pre").c_str(), std::ios::out);
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
    ifs.close();

    if(flag->sbend) {
        std::ifstream ifs(("inputdata_sbend/input_" + oss.str() + ".pre").c_str(), std::ios::out);
        if(! ifs) {
            std::cerr << "input_sbend.pre open error !" << std::endl;
            exit(1);
        } else std::cerr << "input_sbend.pre open" << std::endl;
        for(int i=0; i<3; i++){
            ifs >> data->dset_sbend[i].N;
            for (int j = 0; j < data->dset_sbend[i].N; j++) {
                double w;
                double b1, b2, b3, b4;
                ifs >> w >> b1 >> b2 >> b3 >> b4;
                data->dset_sbend[i].beta_even.append(w, b1);
                data->dset_sbend[i].beta_odd.append(w, b2);
                data->dset_sbend[i].beta_1.append(w, b3);
                data->dset_sbend[i].beta_2.append(w, b4);
            }
        }
        ifs.close();
    }
};

void makeZtoW(DataTable *data){
    Param *par = &(data->par);
    Flag *flag = &(data->flag);

    int count = 0;
    double tmp_g = par->gst;
    for (int i=0; i<par->N_elem; i++){
        std::string elem_name = par->elem[i].elem_name;
        if(elem_name == "dc"){
            double l = par->elem[i].L;
            double cur_g = tmp_g;
            int pcm = par->elem[i].pcm;
            int n_div = Div(l, par->dz);

            for (int step = 0; step < n_div; step++) {
                double cur_z = double(count) * par->dz;
                par->ZtoW.push_back(std::make_tuple(cur_z, cur_g, 1, pcm));
                count++;
            }
        } else if(elem_name == "sbend"){
            double l = par->elem[i].L;
            double gst = tmp_g;
            double gfi = par->elem[i].gfi;
            double w = par->w;
            int pcm = par->elem[i].pcm;
            int n_div = Div(l, par->dz);
            
            if(gst < gfi){
                double sed = (gfi-gst) / 2;
                double r = l*l / (4*sed) * (1 + sed*sed / (l*l));
                double theta = asin(l/(2*r));
                double z_mid = (r+w/2) * sin(theta);
                
                for (int step = 0; step < n_div; step++) {
                    double cur_z = double(count) * par->dz;
                    double stepz = double(step) * par->dz;
                    if(stepz < z_mid){
                        double cur_theta = asin(stepz / (r+w/2));
                        double cur_g = gst + 2*(r+w/2)*(1 - cos(cur_theta));
                        par->ZtoW.push_back(std::make_tuple(cur_z, cur_g, 2, pcm));
                    } else {
                        double stepz_inv = l - stepz;
                        double cur_theta = asin(stepz_inv / (r-w/2));
                        double cur_g = gst + 2*(sed - (r-w/2)*(1 - cos(cur_theta)));
                        par->ZtoW.push_back(std::make_tuple(cur_z, cur_g, 2, pcm));
                    }
                    count++;
                }
            } else {
                double sed = (gst-gfi) / 2;
                double r = l*l / (4*sed) * (1 + sed*sed / (l*l));
                double theta = asin(l/(2*r));
                double z_mid = (r-w/2) * sin(theta);

                for (int step = 0; step < n_div; step++) {
                    double cur_z = double(count) * par->dz;
                    double stepz = double(step) * par->dz;
                    if(stepz < z_mid){
                        double cur_theta = asin(stepz / (r-w/2));
                        double cur_g = gfi + 2*(sed - (r-w/2)*(1 - cos(cur_theta)));
                        par->ZtoW.push_back(std::make_tuple(cur_z, cur_g, 2, pcm));
                    } else {
                        double stepz_inv = l - stepz;
                        double cur_theta = asin(stepz_inv / (r+w/2));
                        double cur_g = gfi + 2*(r+w/2)*(1 - cos(cur_theta));
                        par->ZtoW.push_back(std::make_tuple(cur_z, cur_g, 2, pcm));
                    }
                    count++;
                }
            }
            tmp_g = gfi;
        } else {
            std::cerr << "element is unknown!" << std::endl;
            exit(1);
        }
    }
    double cur_z = double(count) * par->dz; double cur_g = par->elem[par->N_elem-1].gfi;
    if(par->elem[par->N_elem-1].elem_name == "dc") par->ZtoW.push_back(std::make_tuple(cur_z, cur_g, 1, par->elem[par->N_elem-1].pcm));
    else par->ZtoW.push_back(std::make_tuple(cur_z, cur_g, 2, par->elem[par->N_elem-1].pcm));

    if(flag->zw){
        std::ofstream ofs;
        ofs.open("ZtoW", std::ios::out);
        if(! ofs) {
            std::cerr << "File(output) open error !" << std::endl;
            exit(1);
        }
        ofs << "z[um]\tg[um]\tpcm(0:si, 1:cpcm, 2:apcm)\n";
        ofs << std::fixed << std::setprecision(10);
        for (int i=0; i<static_cast<int>(par->ZtoW.size()); i++){
            ofs << std::get<0>(par->ZtoW[i]) << "\t" << std::get<1>(par->ZtoW[i]) << "\t" << std::get<2>(par->ZtoW[i]) << "\t" << std::get<3>(par->ZtoW[i]) << std::endl;
        }
        ofs.close();
    }
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
    int n_div = Div(std::get<0>(par->ZtoW.back()), par->dz);

    for (int step = 0; step < n_div; step++) {
        auto z = std::get<0>(par->ZtoW[step]);
        auto w = std::get<1>(par->ZtoW[step]);
        auto str_flag = std::get<2>(par->ZtoW[step]);
        auto p = std::get<3>(par->ZtoW[step]);
        int flag = p - 1;
        if(flag < 0) {
            std::cerr << "Flag error !" << std::endl;
		    exit(1);
        }

        double be, bo, b1, b2;
        if(str_flag == 1){
            be = data->dset[flag].beta_even[w];
            bo = data->dset[flag].beta_odd[w];
            b1 = data->dset[flag].beta_1[w];
            b2 = data->dset[flag].beta_2[w];
        }else {
            be = data->dset_sbend[flag].beta_even[w];
            bo = data->dset_sbend[flag].beta_odd[w];
            b1 = data->dset_sbend[flag].beta_1[w];
            b2 = data->dset_sbend[flag].beta_2[w];
        }
        

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

    ofs << n_div << '\t' << std::get<0>(par->ZtoW.back()) << '\t' << abs(ab(0)) * abs(ab(0)) << '\t' << abs(ab(1)) * abs(ab(1)) << std::endl;
    ofs.close();
};