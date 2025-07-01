#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cmath>

using namespace std;

int Div(double x, double y){ return((int)(round(x*100000)/round(y*100000)));}

double curve1_1(double R, double r, double g, double w, double sep, double l, double x){
    double y = g/2 + w/2 + r - sqrt(R*R - x*x);
    return y;
}
double curve1_2(double R, double r, double g, double w, double sep, double l, double x){
    double y = g/2 + w/2 + sep - r + sqrt(R*R - (x-l)*(x-l));
    return y;
}

int main(int argc, char* argv[]){
    int j=1;
    double w = atof(argv[j++]);
    double wpcm = atof(argv[j++]);
    double gst = atof(argv[j++]);
    double gfi = atof(argv[j++]);
    double l = atof(argv[j++]);

    double dz = 0.010;
    int N_x = Div(l, dz);
    cout << N_x << endl;

    fstream ofs;
    ofs.open("wg_param.dat", std::ios::out);
    if(! ofs) {
        std::cerr << "File open error !" << std::endl;
        exit(1);
    }
    
    ofs << fixed << setprecision(3);
    
    if(gst < gfi){
        double sep = (gfi - gst)/2;
        double r = (l*l + sep*sep)/(4 * sep);
        double theta = asin(l/(2*r));
        cout << sep << "\t" << r << "\t" << theta * 180 / M_PI << endl;
        vector<double> R = {r-w/2, r-wpcm/2, r+wpcm/2, r+w/2};

        for (int i=0; i<N_x; i++){
            double x = i * dz;
            
            double g, ws, wp;
            
            if(x < R[0]*sin(theta)){
                g = 2 * curve1_1(R[3], r, gst, w, sep, l, x);
                ws = curve1_1(R[0], r, gst, w, sep, l, x) - curve1_1(R[3], r, gst, w, sep, l, x);
                wp = curve1_1(R[1], r, gst, w, sep, l, x) - curve1_1(R[2], r, gst, w, sep, l, x);
            } else if(x < R[1]*sin(theta)){
                g = 2 * curve1_1(R[3], r, gst, w, sep, l, x);
                ws = curve1_2(R[3], r, gst, w, sep, l, x) - curve1_1(R[3], r, gst, w, sep, l, x);
                wp = curve1_1(R[1], r, gst, w, sep, l, x) - curve1_1(R[2], r, gst, w, sep, l, x);
            } else if(x < R[2]*sin(theta)){
                g = 2 * curve1_1(R[3], r, gst, w, sep, l, x);
                ws = curve1_2(R[3], r, gst, w, sep, l, x) - curve1_1(R[3], r, gst, w, sep, l, x);
                wp = curve1_2(R[2], r, gst, w, sep, l, x) - curve1_1(R[2], r, gst, w, sep, l, x);
            } else if(x < R[3]*sin(theta)){
                g = 2 * curve1_1(R[3], r, gst, w, sep, l, x);
                ws = curve1_2(R[3], r, gst, w, sep, l, x) - curve1_1(R[3], r, gst, w, sep, l, x);
                wp = curve1_2(R[2], r, gst, w, sep, l, x) - curve1_2(R[1], r, gst, w, sep, l, x);
            } else {
                g = 2 * curve1_2(R[0], r, gst, w, sep, l, x);
                ws = curve1_2(R[3], r, gst, w, sep, l, x) - curve1_2(R[0], r, gst, w, sep, l, x);
                wp = curve1_2(R[2], r, gst, w, sep, l, x) - curve1_2(R[1], r, gst, w, sep, l, x);
                
            }
            ofs << g << "\t" << ws << "\t" << wp << endl;


    
        }

    }
    ofs.close();
    return 0;
}