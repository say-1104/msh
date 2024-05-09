#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
using namespace std;

ofstream ofs;

#define PI		    3.141592653589793

#define ZERO -1.0e-6

typedef struct _Function {
    int flag;
} Function;

void Start_End(Function *func, int num){
	if(func->flag == num) return;
	else {
		switch(func->flag){
			case 1:
				ofs << "Point END\n" << endl;
				break;
			case 2:
				ofs << "Line END\n" << endl;
				break;
			case 3:
				ofs << "Surface END\n" << endl;
				break;
			case 4:
				ofs << "Volume END\n" << endl;
				break;
			case 5:
				ofs << "Copy END\n" << endl;
				break;
			case 6:
				ofs << "Rotation END\n" << endl;
				break;
			case 7:
				ofs << "Material END\n" << endl;
				break;
			case 8:
				ofs << "Unstr_mesh END\n" << endl;
				break;
			case 9:
				ofs << "Str_mesh END\n" << endl;
				break;
		}

		switch(num){
			case 1:
				ofs << "Point" << endl;
				break;
			case 2:
				ofs << "Line" << endl;
				break;
			case 3:
				ofs << "Surface" << endl;
				break;
			case 4:
				ofs << "Volume" << endl;
				break;
			case 5:
				ofs << "Copy" << endl;
				break;
			case 6:
				ofs << "Rotation" << endl;
				break;
			case 7:
				ofs << "Material" << endl;
				break;
			case 8:
				ofs << "Unstr_mesh" << endl;
				break;
			case 9:
				ofs << "Str_mesh" << endl;
				break;
		}
		func->flag = num;
	}
}


void Point(double x, double y, double z, Function *func){
	Start_End(func, 1);
	ofs << x << " " << y << " " << z << endl;
}

void Trans_Gene(double trans, double gene, Function *func){
	Start_End(func, 0);
	ofs << "TransitionFactor" << endl;
	ofs << trans << " 0\n" << endl;
	ofs << "Generate" << endl;
	ofs << gene << "\n" << endl;
}

void Filename(string lpsname, Function *func){
	Start_End(func, 0);
	ofs << "saveMesh" << endl;
	ofs << lpsname << endl;
	ofs << "exportMesh" << endl;
	ofs << lpsname << endl;
}


