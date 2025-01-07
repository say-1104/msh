#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>

using namespace std;

int main(int argc, char* argv[]){
    int N_wl = 201;
    double MDlambda = 0.5;
    double tmpl;
    for (int i=0; i<N_wl; i++){
            double l, t;
            cin >> l >> t;
            if(MDlambda > t) {
                tmpl = l;
                MDlambda = t;
            }
            
    }
    fstream ofs;
    ofs.open("MDlambda", std::ios::out);
	if(! ofs) {
		std::cerr << "File(MDlambda) open error !" << std::endl;
		exit(1);
    }
    ofs << fixed << setprecision(3);
    ofs << tmpl << "\t";
    ofs << fixed << setprecision(10);
    ofs << MDlambda << endl;
    ofs.close();
    return 0;
}