#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>

using namespace std;

int main(int argc, char* argv[]){
    string filename = argv[1];
    std::ifstream ifs(("result_" + filename + ".dat").c_str(), std::ios::out);
	if(! ifs) {
		std::cerr << "File(result.dat) open error !" << std::endl;
		exit(1);
    }

    int N_wl = 5;
    int N_t = 6;
    int N_l = 251;
    vector<double> t_ideal = {1.0, 0.8, 0.6, 0.4, 0.2, 0};
    vector<string> wl = {"1.53", "1.54", "1.55", "1.56", "1.57"};
    vector<double> l_min(N_t);
    vector<vector<double>> t_min(N_t, vector<double>(N_wl, 0.0));
    vector<double> md_min(N_t, 1.0);

    double tmp;
    for (int i=0;i<N_wl;i++){
        ifs >> tmp;
    }

    for(int i=0;i<N_l;i++){
        double l; ifs >> l;

        vector<double> t_now;
        for (int j=0;j<N_wl;j++){
            ifs >> tmp;
            t_now.push_back(tmp);
        }
        
        for (int j=0;j<N_t;j++){
            double md = 0.0;
            for (int k=0;k<N_wl;k++) if(abs(t_ideal[j]-t_now[k]) > md) md = abs(t_ideal[j]-t_now[k]);
            if(md < md_min[j]) {
                md_min[j] = md;
                l_min[j] = l;
                for (int k=0;k<N_wl;k++) t_min[j][k] = t_now[k];
            }
        }
    }
    ifs.close();

    fstream ofs;
    ofs.open(("glaph_" + filename + ".dat").c_str(), std::ios::out);
    if(! ofs) {
        std::cerr << "File(glaph_result.dat) open error !" << std::endl;
        exit(1);
    }
    ofs << "\t1.0\t0.8\t0.6\t0.4\t0.2\t0.0" << endl;
    for (int j=0;j<N_wl;j++){
        ofs << wl[j];
        for (int k=0;k<N_t;k++) ofs << "\t" << t_min[k][j];
        ofs << endl;
    }
    for (int k=0;k<N_t;k++) ofs << "\t" << md_min[k]; ofs << endl;
    for (int k=0;k<N_t;k++) ofs << "\t" << l_min[k]; ofs << endl;
    ofs.close();

    return 0;
}