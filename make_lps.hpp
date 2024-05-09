#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
using namespace std;

ofstream ofs;

#define PI		    3.141592653589793

#define ZERO -1.0e-6

void Trans_Gene(double trans, double gene){
	ofs << "TransitionFactor" << endl;
	ofs << trans << " 0\n" << endl;
	ofs << "Generate" << endl;
	ofs << gene << endl;
}

void Filename(string lpsname){
	ofs << "saveMesh" << endl;
	ofs << lpsname << endl;
	ofs << "exportMesh" << endl;
	ofs << lpsname << endl;
}