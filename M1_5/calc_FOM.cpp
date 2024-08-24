#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>

using namespace std;

int main(int argc, char* argv[]){
    int N_wl = 5;
    int flag;
    cin >> flag;
    double FOM = 0.0;
    if(flag == 1){
        for (int i=0; i<N_wl; i++){
            double l, w, pa, pb;
            cin >> l >> w >> pb >> pa;
            double psr = pb / (pa + pb);
            FOM += psr;
        }
        for (int i=0; i<N_wl; i++){
            double l, w, pa, pb;
            cin >> l >> w >> pb >> pa;
            double psr = pb / (pa + pb);
            FOM += (1-psr);
        }
    }
    fstream ofs;
    ofs.open("FOM", std::ios::out);
	if(! ofs) {
		std::cerr << "File(FOM) open error !" << std::endl;
		exit(1);
    }

    ofs << fixed << setprecision(10);
    ofs << FOM << endl;
    ofs.close();
    return 0;
}