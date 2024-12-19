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

    int i=1;

    double dz = atof(argv[i++]);
    double wst = atof(argv[i++]);
    int N_taper = atof(argv[i++]);

    ofs << fixed << setprecision(1);
    ofs << "dz\t" << dz << endl;
    ofs << fixed << setprecision(3);
    ofs << "wst\t" << wst << endl << endl;

    ofs << "N_taper\t" << N_taper << endl;
    ofs << "z\twidth" << endl;
    for (int _=0; _<N_taper; _++){
        double z = atof(argv[i++]);
        double w = atof(argv[i++]);

        ofs << fixed << setprecision(1);
        ofs << z;
        ofs << "\t" ;
        ofs << fixed << setprecision(3);
        ofs << w << endl;
    }
    ofs  << endl;

    int N_pcm = atof(argv[i++]);
    ofs << "N_pcm\t" << N_pcm << endl;
    ofs << "z\tpcm" << endl;
    for (int _=0; _<N_pcm; _++){
        double z = atof(argv[i++]);
        int pcm = atof(argv[i++]);

        ofs << fixed << setprecision(1);
        ofs << z;
        ofs << "\t" ;
        ofs << pcm << endl;
    }

    ofs.close();
    return 0;
}