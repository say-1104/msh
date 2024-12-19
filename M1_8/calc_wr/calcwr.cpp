#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <iomanip>
using namespace std;

int main(){
    int N_si = 25;
    int N_pcm = 100;
    string buff;
    vector<vector<double>> si(N_si, vector<double>(2));
    vector<vector<double>> pcm(N_pcm, vector<double>(2));
    cin >> buff >> buff;
    for(int i=0; i<N_si; i++){
        si[i][0] = 451+i;
        cin >> buff >> si[i][1];
    }
    cin >> buff >> buff;
    for(int i=0; i<N_pcm; i++){
        pcm[i][0] = 401+i;
        cin >> buff >> pcm[i][1];
    }

    std::ofstream ofs;
    ofs.open("Wh.dat", std::ios::out);
	if(! ofs) {
		std::cerr << "File(output) open error !" << std::endl;
		exit(1);
        
    }
    ofs << "Wr\tWh\n";
    ofs << std::fixed << std::setprecision(3);

    for(int i=0; i<N_si; i++){
        double min = 10000;
        int min_key = 0;
        int wr = si[i][0];
        for(int j=0; j<N_pcm; j++){
            double tmp = abs(si[i][1] - pcm[j][1]);
            if(tmp < min){
                min = tmp;
                min_key = pcm[j][0];
            }

        }
        ofs <<"0." << wr << "\t0." <<  min_key << endl;

    }
    ofs.close();
}