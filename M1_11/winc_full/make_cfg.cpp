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
    double gst = atof(argv[i++]);
    double w = atof(argv[i++]);
    int N_element = atof(argv[i++]);

    ofs << fixed << setprecision(1);
    ofs << "dz\t" << dz << endl;
    ofs << fixed << setprecision(3);
    ofs << "gst\t" << gst << endl;
    ofs << "w\t" << w << endl << endl;
    ofs << "N_element\t" << N_element << endl;
    for (int _=0; _<N_element; _++){
        ofs << endl;
        string elem_name = argv[i++];
        ofs << elem_name << endl;
        
        if(elem_name == "sbend") {
            double L = atof(argv[i++]);
            double gfi = atof(argv[i++]);
            int pcm = atof(argv[i++]);

            ofs << fixed << setprecision(1);
            ofs << "L\t" << L << endl;
            ofs << fixed << setprecision(3);
            ofs << "gfi\t" << gfi << endl;
            ofs << "pcm " << pcm << endl;
        } else if(elem_name == "dc"){
            double L = atof(argv[i++]);
            int pcm = atof(argv[i++]);

            ofs << fixed << setprecision(1);
            ofs << "L\t" << L << endl;
            ofs << "pcm " << pcm << endl;
        }
        else {
            std::cerr << "element name is unknown !" << std::endl;
            exit(1);
        }
    }

    ofs.close();
    return 0;
}