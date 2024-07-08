#include "main.hpp"

extern int calcCMT(DataTable *data){
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
    int n_div = par->taper[par->N_taper].second / par->dz;

    for (int step = 0; step < n_div; step++) {
        auto z = std::get<0>(par->ZtoW[step]);
        auto w = std::get<1>(par->ZtoW[step]);
        auto flag = std::get<2>(par->ZtoW[step]);

        double be = data->dset->beta_even[w];
        double bo = data->dset->beta_odd[w];
        double b1 = data->dset->beta_1[w];
        double b2 = data->dset->beta_2[w];

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

    return 0;
};