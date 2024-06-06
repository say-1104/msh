#include "make_lps.hpp"

void Printvv(vector<vector<int> > a){
	for(int i=0; i<a.size(); i++){
		for(int j=0; j<a[i].size(); j++){
			cout << a[i][j] << "\t" ;
		}
		cout << endl;
	}
	cout << endl;
}

vector<int> Plusv(vector<int> a, int n){
	for(int i=0; i<a.size(); i++){
		a[i] += n;
	}
	return a;
}

int Div(double x, double y){
	return((int)(round(x*10000)/round(y*10000)));
}

int main(int argc, char *argv[]){
	Function func;
	func.flag = 0;
	func.tp_num = 0;
	int i, j, k;
	
	ofs.open("pcmDC.lps");
	if(! ofs) {
		cerr << "File open error !" << endl;
		exit(1);
	}
	ofs << fixed << setprecision(10);
	
	//DCのパラメータ(外部入力)
	double Wh = atof(argv[1]);		//導波路幅(PCM有)
	double Wr = atof(argv[2]);		//導波路幅(PCM無)
	double Wpcm = atof(argv[3]);	//PCM幅
	double g = atof(argv[4]);		//導波路間隔
	double Hpcm = atof(argv[5]);	//PCM層の厚さ
	double Lc = atof(argv[6]);		//PCM層の長さ
	double Leff = atof(argv[7]);	//cPCM層の長さ

	//DCのパラメータ(固定値)
	double Hsi = 0.22;				//Si層の厚さ
	double Hsio2 = 2.0;				//SiO2層の厚さ
	double Wm = 2.0;				//導波路とPML間の距離
	double Wpml = 1.0;				//PML幅
	
	int div = 10;
	double dLc = Lc / div;
	//double dLc = 13.3;
	//Div(Lc, dLc);

	double L = 2 * Wpml + 2 * Wm + Lc;
	double W = 2 * Wpml + 2 * Wm + Wr + g + Wh;
	double H = 2 * Wpml + Hsio2 + Hsi + Hpcm + Wm;



	//メッシュのパラメータ
	double unstr = 0.05, trans = 0.5, gene = 0.5;
	string lpsname = "pcmDC";

	cout << "div: " << div << endl;
	cout << "writing " << lpsname << "..." << endl;

	Tppush(0, 0, 0, 0, &func);

	vector<int> v1, v2, v3, v4, v5, v6;

	auto msh_xy = [&](int tp_num, double Z) -> void {
		func.tp_num = tp_num; 
		Line(0.0, 0.0, Z, 	Wpml, 0.0, Z, &func);
		Copy(Poi, Wm, 0.0, 0.0, 			"2", &func);
		Copy(Poi, Wr, 0.0, 0.0, 			"3", &func);
		Copy(Poi, g, 0.0, 0.0, 				"4", &func);
		Copy(Poi, (Wh-Wpcm)/2, 0.0, 0.0, 	"5", &func);
		Copy(Poi, Wpcm, 0.0, 0.0, 			"6", &func);
		Copy(Poi, (Wh-Wpcm)/2, 0.0, 0.0, 	"7", &func);
		Copy(Poi, Wm, 0.0, 0.0, 			"8", &func);
		Copy(Poi, Wpml, 0.0, 0.0, 			"9", &func);

		Copy(Lin, 0.0, Wpml, 0.0, 			"1 2 3 4 5 6 7 8 9", &func);
		Copy(Lin, 0.0, Hsio2, 0.0, 			"10 11 12 13 14 15 16 17 18", &func);
		Copy(Lin, 0.0, Hsi, 0.0, 			"29 30 31 32 33 34 35 36 37", &func);
		Copy(Lin, 0.0, Hpcm, 0.0, 			"48 49 50 51 52 53 54 55 56", &func);
		Copy(Lin, 0.0, Wm, 0.0, 			"67 68 69 70 71 72 73 74 75", &func);
		Copy(Lin, 0.0, Wpml, 0.0, 			"86 87 88 89 90 91 92 93 94", &func);

		Tppush(0, 54, 123, 70, &func);
	};

	auto msh_z = [&](int tp_num, double Z) -> void {
		func.tp_num = tp_num;
		Copy(Sur, 0.0, 0.0, Z,  			"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54", &func);
		Tppush(54, 123, 70, 0, &func);
	};
	//1-3
	msh_xy(1, 0.0);
	msh_xy(2, Wpml);
	msh_z(1, Wpml);
	//4-5
	msh_xy(4, Wpml+Wm);
	msh_z(2, Wm);

	for(i=0; i<div; i++){
		msh_xy(6+2*i, Wpml + Wm + dLc*(i+1));
		msh_z(4+2*i, dLc);
	}

	msh_xy(6+2*div, L - Wpml);
	msh_z(4+2*div, Wm);

	msh_xy(8+2*div, L);
	msh_z(6+2*div, Wpml);
	
	
	v1 = {1, 2, 3, 4 ,5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
	v2 = {19, 20, 22, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54};
	v3 = {21};
	v4 = {23, 24, 25};
	Mat3D(1, v1, &func);
	Mat3D(2, v2, &func);
	Mat3D(3, v3, &func);
	Mat3D(4, v4, &func);
	v1 = Plusv(v1, func.tp[0][9+2*div]);
	v2 = Plusv(v2, func.tp[0][9+2*div]);
	v3 = Plusv(v3, func.tp[0][9+2*div]);
	v4 = Plusv(v4, func.tp[0][9+2*div]);
	Mat3D(5, v1, &func);
	Mat3D(6, v2, &func);
	Mat3D(7, v3, &func);
	Mat3D(8, v4, &func);

	v1 = {};
	v2 = {};
	v3 = {};
	v4 = {33+func.tp[0][5], 33+func.tp[0][7+2*div]};
	for(i=0;i<div+2;i++){
		v5 = {};
		v6 = {};
		v1.push_back(21+func.tp[0][5+2*i]);
		v2.push_back(23+func.tp[0][5+2*i]);
		v2.push_back(24+func.tp[0][5+2*i]);
		v2.push_back(25+func.tp[0][5+2*i]);
		v5 = Plusv({1, 2, 3, 4 ,5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18}, func.tp[0][5+2*i]);
		copy(v5.begin(),v5.end(), back_inserter(v3));
		v6 = Plusv({19, 20, 22, 26, 27, 28, 29, 30, 31, 32, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54}, func.tp[0][5+2*i]);
		copy(v6.begin(),v6.end(), back_inserter(v4));
	}
	Mat3D(9, v1, &func);
	Mat3D(10, v2, &func);
	Mat3D(11, v3, &func);
	Mat3D(12, v4, &func);
	for (i=0; i<div; i++){
		v1={33+func.tp[0][7+2*i]};
		Mat3D(13+i, v1, &func);
	}
	v1={21, 22, 23, 24, 25, 32, 33, 34};
	for (i=0; i<div+3; i++) {
		v5 = {};
		v5 = Plusv({21, 22, 23, 24, 25, 32, 33, 34}, func.tp[1][2*(i+1)]);
		copy(v5.begin(),v5.end(), back_inserter(v1));
		v6 = {};
		v6 = Plusv({31, 32, 33, 34, 35, 50, 51, 52, 53, 54, 59, 60, 61, 62, 63, 64, 71, 72, 73, 80, 81, 82, 83}, func.tp[1][2*(i+1)+1]);
		copy(v6.begin(),v6.end(), back_inserter(v1));
	}
	Unstr(Sur, unstr, v1, &func);
	
	Trans_Gene(trans, gene, &func);
	Printvv(func.tp);

	Filename(lpsname, &func);

	ofs.close();

	return 0;
}