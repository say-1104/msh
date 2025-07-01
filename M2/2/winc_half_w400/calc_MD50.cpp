#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>

using namespace std;

int main(int argc, char* argv[]){
    int N_wl = 5;
    double MDlambda = 0.0;
    double wish = 0.5;
    for (int i=0; i<N_wl; i++){
            double l, t;
            cin >> l >> t;
            MDlambda = max(MDlambda, abs(t-wish));
    }
    fstream ofs;
    ofs.open("MDlambda50", std::ios::out);
    if(! ofs) {
        std::cerr << "File(MDlambda50) open error !" << std::endl;
        exit(1);
    }
    
    ofs << fixed << setprecision(10);
    ofs << MDlambda << endl;
    ofs.close();

    return 0;
}