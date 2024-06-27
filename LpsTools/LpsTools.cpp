#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <LpsTools.hpp>

std::vector<int> Plusv(std::vector<int> a, int n){
	for(int i=0; i<a.size(); i++){
		a[i] += n;
	}
	return a;
}

LpsTools::LpsTools(const char* name, double unstr, double trans, double gene) : lpsname(name), unstr(unstr), trans(trans), gene(gene) {
    currentfunc = Function::none;
    y_offset = 0.0;
    step_offset = 0;
    Appendstep(0, 0, 0, 0);

    ofs.open((lpsname + ".lps").c_str(), std::ios::out);
	if(! ofs) {
		std::cerr << "File open error !" << std::endl;
		exit(1);
	}
	ofs << std::fixed << std::setprecision(10);

    std::cerr << "Writing " << lpsname << "..." << std::endl;

}

std::string LpsTools::getFuncname(Function func) {
    std::string buff = "";
    switch(func) {
        case Function::point:
            buff = "Point";
            break;

        case Function::line:
            buff = "Line";
            break;

        case Function::surface:
            buff = "Surface";
            break;
            
        case Function::volume:
            buff = "Volume";
            break;

        case Function::copy:
            buff = "Copy";
            break;

        case Function::rotation:
            buff = "Rotation";
            break;

        case Function::material:
            buff = "Material";
            break;

        case Function::unstr:
            buff = "Unstr_mesh";
            break;
        
        case Function::str:
            buff = "Str_mesh";
            break;

        default:
            buff = "";
            break;
    }
    return buff;
}

void LpsTools::SwitchingFunc(Function func) {
    std::string buff;
    if(currentfunc != func) {
        if(currentfunc != Function::none) {
            buff = getFuncname(currentfunc);
            ofs << buff << " END" << std::endl;
            ofs << std::endl;
        }
        currentfunc = func;
        if(currentfunc != Function::none) {
            buff = getFuncname(currentfunc);
            ofs << buff << std::endl;
        }
    }

}

std::vector<int> LpsTools::Atovec(const char* text){
	std::string text_str = text;
	std::stringstream ss(text_str);
	std::string buf;
	std::vector<int> v;
	while (std::getline(ss, buf, ' ')) {
		v.push_back(std::stoi(buf));
	}
	return v;
}

void LpsTools::Appendstep(int p, int l, int s, int v){
    int tpp = p + step[0].back();
    int tpl = l + step[1].back();
    int tps = s + step[2].back();
    int tpv = v + step[3].back();

    step[0].push_back(tpp);
    step[1].push_back(tpl);
    step[2].push_back(tps);
    step[3].push_back(tpv);
};

void LpsTools::Printstep() {
    std::cerr << "Point\tLine\tSurface\tVolume" << std::endl;
	for(int i=2; i<step[0].size(); i++){
		for(int j=0; j<4; j++){
			std::cerr << step[j][i] << "\t" ;
		}
		std::cerr << std::endl;
	}
}

void LpsTools::Point(double x, double y, double z){
    SwitchingFunc(Function::point);
    ofs << x << "\t" << y << '\t' << z << std::endl;
}

void LpsTools::Point(double x, double z){
    SwitchingFunc(Function::point);
    ofs << x << "\t" << y_offset << '\t' << z << std::endl;
}

void LpsTools::Line(double x1, double y1, double z1, double x2, double y2, double z2){
    SwitchingFunc(Function::line);
    ofs << "0\t ";
	ofs << x1 << '\t' << y1 << '\t' << z1 << "\t ";
	ofs << x2 << '\t' << y2 << '\t' << z2 << std::endl;
}

void LpsTools::Line(double x1, double z1, double x2, double z2){
    SwitchingFunc(Function::line);
    ofs << "0\t ";
	ofs << x1 << '\t' << y_offset << '\t' << z1 << "\t ";
	ofs << x2 << '\t' << y_offset << '\t' << z2 << std::endl;
}

void LpsTools::Surface(const char* text){
    std::vector<int> v = Atovec(text);
    v = Plusv(v, step[1][step_offset]);
    Surface(v);
}

void LpsTools::Surface(std::vector<int> v){
    SwitchingFunc(Function::surface);
    for(int i=0; i<v.size(); i++) {
        ofs << v[i] << " ";
    }
    ofs << std::endl;
}

void LpsTools::Surface(int st, int end, int inc){
    SwitchingFunc(Function::surface);
    for(int i=st; i<=end; i+=inc) {
        ofs << st + step[1][step_offset] << " ";
    }
    ofs << std::endl;
}

void LpsTools::Volume(const char* text){
    std::vector<int> v = Atovec(text);
    v = Plusv(v, step[2][step_offset]);
    Volume(v);
}

void LpsTools::Volume(std::vector<int> v){
    SwitchingFunc(Function::volume);
    for(int i=0; i<v.size(); i++) {
        ofs << v[i] << " ";
    }
    ofs << std::endl;
}

void LpsTools::Volume(int st, int end, int inc){
    SwitchingFunc(Function::volume);
    for(int i=st; i<=end; i+=inc) {
        ofs << st + step[2][step_offset] << " ";
    }
    ofs << std::endl;
}

void LpsTools::Copy(Shape shape, double x, double y, double z, const char* text){
	std::vector<int> v = Atovec(text);
    int n_shape = static_cast<int>(shape);
    v= Plusv(v, step[n_shape][step_offset]);
    Copy(shape, x, y, z, v);
}

void LpsTools::Copy(Shape shape, double x, double y, double z, std::vector<int> v){
	SwitchingFunc(Function::copy);
	int n_shape = static_cast<int>(shape);

	ofs << n_shape << " " << n_shape + 1 << " 0" << std::endl;
	ofs << "0.000\t0.000\t0.000\t ";
	ofs << x << "\t" << y << "\t" << z << std::endl;
	
	for (int i=0; i<v.size(); i++) {
		ofs << v[i] << " ";
	}
	ofs << std::endl;
}

void LpsTools::Copy(Shape shape, double x, double y, double z, int st, int end, int inc){
	SwitchingFunc(Function::copy);
	int n_shape = static_cast<int>(shape);

	ofs << n_shape << " " << n_shape + 1 << " 0" << std::endl;
	ofs << "0.000\t0.000\t0.000\t ";
	ofs << x << "\t" << y << "\t" << z << std::endl;
	
	for (int i=st; i<=end; i+=inc) {
		ofs << st + step[n_shape][step_offset] << " ";
	}
	ofs << std::endl;
}

void LpsTools::Mat2D(int mat, const char* text){
	std::vector<int> v = Atovec(text);
    Mat2D(mat, v);
}

void LpsTools::Mat2D(int mat, std::vector<int> v){
	SwitchingFunc(Function::material);
	int dim = 2;

	ofs << dim << ' ' << mat << std::endl;

	for (int i=0; i<v.size(); i++) {
		ofs << v[i] << ' ';
	}
	ofs << std::endl;
}

void LpsTools::Mat3D(int mat, const char* text){
	std::vector<int> v = Atovec(text);
    Mat3D(mat, v);
}

void LpsTools::Mat3D(int mat, std::vector<int> v){
	SwitchingFunc(Function::material);
	int dim = 3;

	ofs << dim << ' ' << mat << std::endl;

	for (int i=0; i<v.size(); i++) {
		ofs << v[i] << ' ';
	}
	ofs << std::endl;
}


void LpsTools::Unstr(Shape shape, double unstr, const char* text){
    std::vector<int> v = Atovec(text);
    Unstr(shape, unstr, v);
}

void LpsTools::Unstr(Shape shape, double unstr, std::vector<int> v){
	SwitchingFunc(Function::unstr);
	int tmp = static_cast<int>(shape);
	
	ofs << tmp << ' ' << unstr << std::endl;
	
	for (int i=0; i<v.size(); i++) {
		ofs << v[i] << ' ';
	}
	ofs << std::endl;
}

void LpsTools::TransGene() {
    SwitchingFunc(Function::none);
    ofs << "TransitionFactor" << std::endl;
    ofs << trans << " 0" << std::endl;
    ofs << std::endl;

    ofs << "Generate" << std::endl;
    ofs << gene << std::endl;
    ofs << std::endl;
}

void LpsTools::Fileclose() {
    SwitchingFunc(Function::none);
    ofs << "saveMesh" << std::endl;
    ofs << lpsname << std::endl;
    ofs << "exportMesh" << std::endl;
    ofs << lpsname << std::endl;

    Printstep();
    ofs.close();
    std::cerr << "Writing Finish!" << std::endl;
}