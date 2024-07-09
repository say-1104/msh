#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>

using namespace std;

int main(int argc, char* argv[]){
    fstream ofs;
    ofs.open("structure.cfg", std::ios::out);
	if(! ofs) {
		std::cerr << "File(output) open error !" << std::endl;
		exit(1);
    }
    double dz = atof(argv[1]);
    int N_taper = (argc-1)/2 - 1;
    vector<double> z, w;

    for (int i=0; i<N_taper+1; i++){
        z.push_back(atof(argv[2+i*2]));
        w.push_back(atof(argv[3+i*2]));
    }

    ofs << fixed << setprecision(1);
    ofs << "dz\t" << dz << endl;
    ofs << "N_taper\t" << N_taper << endl;
    ofs << "z\twidth" << endl;
    
    for (int i=0; i<N_taper+1; i++){
        ofs << fixed << setprecision(1);
        ofs << z[i] ;
        ofs << "\t" ;
        ofs << fixed << setprecision(3);
        ofs << w[i] << endl;
    }
    ofs.close();
    return 0;
}