#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iomanip>
using namespace std;

ofstream ofs;

#define PI		    3.141592653589793

typedef struct _Function {
    int flag;
	int tp_num;
	vector<vector<int> > tp = vector<vector<int> >(4, vector<int>(1,0));
} Function;


void Tppush(int v, int s, int l, int p, Function *func) {
	int tpv = v + func->tp[0].back();
	int tps = s + func->tp[1].back();
	int tpl = l + func->tp[2].back();
	int tpp = p + func->tp[3].back();
	func->tp[0].push_back(tpv);
	func->tp[1].push_back(tps);
	func->tp[2].push_back(tpl);
	func->tp[3].push_back(tpp);
}
enum Shape {
	Poi,
	Lin,
	Sur,
	Vol
};

enum Axis {
	X,
	Y,
	Z
};

vector<int> Atovec(const char* text){
	string text_str = text;
	stringstream ss(text_str);
	string buf;
	vector<int> v;
	while (std::getline(ss, buf, ' ')) {
		v.push_back(std::stoi(buf));
	}
	return v;
}

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
	ofs << x << '\t' << y << '\t' << z << endl;
}

void Line(double x1, double y1, double z1, double x2, double y2, double z2, Function *func){
	Start_End(func, 2);
	ofs << "0\t ";
	ofs << x1 << '\t' << y1 << '\t' << z1 << "\t ";
	ofs << x2 << '\t' << y2 << '\t' << z2 << endl;
}

void Curve(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, Function *func){
	Start_End(func, 2);
	ofs << "1\t ";
	ofs << x1 << '\t' << y1 << '\t' << z1 << "\t ";
	ofs << x2 << '\t' << y2 << '\t' << z2 << "\t ";
	ofs << x3 << '\t' << y3 << '\t' << z3 << endl;
}

void Surface(const char* text, Function *func){
	vector<int> v = Atovec(text);
	Start_End(func, 3);
	for (int i=0; i<v.size(); i++) {
		ofs << v[i] << ' ';
	}
	ofs << endl;
}

void Surface(vector<int> v, Function *func){
	Start_End(func, 3);
	for (int i=0; i<v.size(); i++) {
		ofs << v[i] + func->tp[1][func->tp_num] << ' ';
	}
	ofs << endl;
}

void Volume(const char* text, Function *func){
	vector<int> v = Atovec(text);
	Start_End(func, 4);
	for (int i=0; i<v.size(); i++) {
		ofs << v[i] << ' ';
	}
	ofs << endl;
}

void Volume(vector<int> v, Function *func){
	Start_End(func, 4);
	for (int i=0; i<v.size(); i++) {
		ofs << v[i] << ' ';
	}
	ofs << endl;
}

void Copy(Shape shape, double x, double y, double z, const char* text, Function *func){
	vector<int> v = Atovec(text);
	Start_End(func, 5);
	int tmp = shape;

	ofs << tmp << ' ' << tmp+1 << " 0" << endl;
	ofs << "0.000\t0.000\t0.000\t ";
	ofs << x << '\t' << y << '\t' << z << endl;
	
	for (int i=0; i<v.size(); i++) {
		ofs << v[i] + func->tp[3-tmp][func->tp_num] << ' ';
	}
	ofs << endl;
}

void Copy(Shape shape, double x, double y, double z, vector<int> v, Function *func){
	Start_End(func, 5);
	int tmp = shape;

	ofs << tmp << ' ' << tmp+1 << " 0" << endl;
	ofs << "0.000\t0.000\t0.000\t ";
	ofs << x << '\t' << y << '\t' << z << endl;
	
	for (int i=0; i<v.size(); i++) {
		ofs << v[i] + func->tp[3-tmp][func->tp_num] << ' ';
	}
	ofs << endl;
}

void Rotate(Shape shape, Axis axis, double x, double y, double z, double angle, const char* text, Function *func){
	vector<int> v = Atovec(text);
	Start_End(func, 6);
	int tmp = shape;

	ofs << tmp << ' ' << tmp+1 << " 0" << endl;
	ofs << (int)axis << "\t ";
	ofs << x << '\t' << y << '\t' << z << "\t ";
	ofs << angle << endl;
	
	for (int i=0; i<v.size(); i++) {
		ofs << v[i] + func->tp[3-tmp][func->tp_num] << ' ';
	}
	ofs << endl;
}

void Rotate(Shape shape, Axis axis, double x, double y, double z, double angle, vector<int> v, Function *func){
	Start_End(func, 6);
	int tmp = shape;
	
	ofs << tmp << ' ' << tmp+1 << " 0" << endl;
	ofs << (int)axis << "\t ";
	ofs << x << '\t' << y << '\t' << z << "\t ";
	ofs << angle << endl;
	
	for (int i=0; i<v.size(); i++) {
		ofs << v[i] + func->tp[3-tmp][func->tp_num] << ' ';
	}
	ofs << endl;
}

void Mat2D(int mat, const char* text, Function *func){
	vector<int> v = Atovec(text);
	Start_End(func, 7);
	int dim = 2;

	ofs << dim << mat << endl;

	for (int i=0; i<v.size(); i++) {
		ofs << v[i] << ' ';
	}
	ofs << endl;
}

void Mat2D(int mat, vector<int> v, Function *func){
	Start_End(func, 7);
	int dim = 2;

	ofs << dim << mat << endl;

	for (int i=0; i<v.size(); i++) {
		ofs << v[i] << ' ';
	}
	ofs << endl;
}

void Mat3D(int mat, const char* text, Function *func){
	vector<int> v = Atovec(text);
	Start_End(func, 7);
	int dim = 3;

	ofs << dim << mat << endl;

	for (int i=0; i<v.size(); i++) {
		ofs << v[i] << ' ';
	}
	ofs << endl;
}

void Mat3D(int mat, vector<int> v, Function *func){
	Start_End(func, 7);
	int dim = 3;

	ofs << dim << mat << endl;

	for (int i=0; i<v.size(); i++) {
		ofs << v[i] << ' ';
	}
	ofs << endl;
}

void Unstr(Shape shape, double unstr, const char* text, Function *func){
	vector<int> v = Atovec(text);
	Start_End(func, 8);
	int tmp = shape;
	
	ofs << tmp << ' ' << unstr << endl;
	
	for (int i=0; i<v.size(); i++) {
		ofs << v[i] << ' ';
	}
	ofs << endl;
}

void Unstr(Shape shape, double unstr, vector<int> v, Function *func){
	Start_End(func, 8);
	int tmp = shape;
	
	ofs << tmp << ' ' << unstr << endl;
	
	for (int i=0; i<v.size(); i++) {
		ofs << v[i] << ' ';
	}
	ofs << endl;
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


