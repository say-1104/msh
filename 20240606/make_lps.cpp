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

	//DCのパラメータ(固定値)
	double Hsi = 0.22;				//Si層の厚さ
	double Hsio2 = 2.0;				//SiO2層の厚さ
	double Wx = 1.5;				//導波路とPML間の距離
	double Wy = 1.5;				//導波路とPML間の距離
	double Wz = 0.5;				//導波路とPML間の距離
	double Wpml = 1.0;				//PML幅
	
	int div = 10;
	double dLc = Lc / div;
	//double dLc = 13.3;
	//int div = Div(Lc, dLc);

	double L = 2 * Wpml + 2 * Wz + Lc;
	double W = 2 * Wpml + 2 * Wx + Wr + g + Wh;
	double H = 2 * Wpml + Hsio2 + Hsi + Hpcm + Wy;



	//メッシュのパラメータ
	double unstr = 0.05, trans = 0.5, gene = 0.5;
	string lpsname = "pcmDC";

	cout << "writing " << lpsname << "..." << endl;

	Tppush(0, 0, 0, 0, &func);

	vector<int> v1, v2, v3, v4, v5, v6;

	auto lpsxy = [&](int tp_num, double Y) -> void {
		func.tp_num = tp_num;
		Line(0.0, Y, 0.0, 				0.0, Y, Wpml, &func);
		Line(0.0, Y, Wpml, 				0.0, Y, L-Wpml, &func);
		Line(0.0, Y, L-Wpml, 			0.0, Y, L, &func);
		Copy(Lin, Wpml, 0.0, 0.0, 		"1 2 3", &func);
		Copy(Lin, W-2*Wpml, 0.0, 0.0, 	"4 5 6", &func);
		Copy(Lin, Wpml, 0.0, 0.0, 		"11 12 13", &func);
		Tppush(0, 9, 24, 16, &func);
	};

	auto lpsz = [&](int tp_num, double Y) -> void {
		func.tp_num = tp_num;
		Copy(Sur, 0.0, Y, 0.0, 		"1 2 3 4 5 6 7 8 9", &func);
		Tppush(9, 24, 16, 0, &func);
	};

	//1, 2
	lpsxy(1, 0.0);
	lpsxy(2, Wpml);

	//3
	lpsz(0.0, Wpml);
	
	//4, 5
	lpsxy(4, H-Wpml);
	lpsxy(5, H);

	//6
	lpsz(0.0, Wpml);
	
	//7
	func.tp_num = 7;
	vector <int> p1, l1, l2, l3, l4;
	Line(0.0, Wpml+Hsio2, 0.0, 0.0, Wpml+Hsio2, Wpml, &func);
	Copy(Lin, Wpml, 0.0, 0.0, "1", &func);
	Copy(Lin, Wx, 0.0, 0.0, "2", &func);
	Copy(Lin, Wr, 0.0, 0.0, "5", &func);
	Copy(Lin, g, 0.0, 0.0, "8", &func);
	Copy(Lin, Wh, 0.0, 0.0, "11", &func);
	Copy(Lin, Wx, 0.0, 0.0, "14", &func);
	
	Line(0.0, Wpml+Hsio2, L-Wpml, 0.0, Wpml+Hsio2, L, &func);
	Copy(Lin, Wpml, 0.0, 0.0, "17", &func);
	Copy(Lin, Wx, 0.0, 0.0, "18", &func);
	Copy(Lin, Wr, 0.0, 0.0, "21", &func);
	Copy(Lin, g, 0.0, 0.0, "24", &func);
	Copy(Lin, Wh, 0.0, 0.0, "27", &func);
	Copy(Lin, Wx, 0.0, 0.0, "30", &func);
	Copy(Lin, Wpml, 0.0, 0.0, "33", &func);
	
	Copy(Lin, 0.0, 0.0, L-2*Wpml, "4 16", &func);
	
	Copy(Lin, 0.0, 0.0, Wz, "10", &func);
	Copy(Lin, 0.0, 0.0, L-2*Wpml-2*Wz, "49", &func);
	Copy(Lin, 0.0, 0.0, Wz, "52", &func);

	Copy(Poi, 0.0, 0.0, Wm, "10 12", &func);
	Copy(Poi, (Wh-Wpcm)/2, 0.0, 0.0, "37", &func);
	Copy(Poi, Wpcm, 0.0, 0.0, "39", &func);
	Copy(Poi, (Wh-Wpcm)/2, 0.0, 0.0, "40", &func);
	Surface("16 57 59 60 61 58", &func);
	Copy(Lin, 0.0, 0.0, dLc, "59 60 61", &func);
	for (i=0; i<div-1; i++) {
		v1={62+i*7, 63+i*7, 64+i*7};
		Copy(Lin, 0.0, 0.0, dLc, v1, &func);
	}
	p1={37+4*div, 38+4*div};
	l1={62+(div-1)*7, 63+(div-1)*7, 64+(div-1)*7};
	l2={l1[2]+5, l1[2]+6};
	v2={l1[0]+func.tp[2][7], l1[1]+func.tp[2][7], l1[2]+func.tp[2][7], l2[0]+func.tp[2][7], l2[1]+func.tp[2][7], 37+func.tp[2][7]};
	Copy(Poi, 0.0, 0.0, Wm, p1, &func);
	Surface(v2, &func);

	Surface("7 46 28 55 53 50", &func);
	v3={57+func.tp[2][7], 13+func.tp[2][7], 51+func.tp[2][7], 54+func.tp[2][7], 56+func.tp[2][7], 34+func.tp[2][7]};
	v4={58+func.tp[2][7], 19+func.tp[2][7], 47+func.tp[2][7], 40+func.tp[2][7]};
	for (i=0; i<div; i++) {
		v3.push_back(65+i*7+func.tp[2][7]);
		l3.push_back(65+i*7);
		v4.push_back(66+i*7+func.tp[2][7]);
		l4.push_back(66+i*7);
	}
	v3.push_back(65+div*7-3+func.tp[2][7]);
	v4.push_back(66+div*7-3+func.tp[2][7]);
	Surface(v3, &func);
	Surface(v4, &func);
	Tppush(0, 24+div*3, 63+div*7, 40+div*4, &func);
	
	//8
	func.tp_num = 2;
	Copy(Sur, 0.0, Hsio2, 0.0, "1 2 3 7 8 9", &func);
	v1={14+func.tp[2][2], 5+func.tp[2][8], 9+func.tp[2][8], 6+func.tp[2][7], 9+func.tp[2][7], 12+func.tp[2][7], 15+func.tp[2][7], 18+func.tp[2][7]};
	v2={15+func.tp[2][2], 6+func.tp[2][8], 10+func.tp[2][8], 7+func.tp[2][7], 10+func.tp[2][7], 13+func.tp[2][7], 16+func.tp[2][7], 19+func.tp[2][7]};
	v3={16+func.tp[2][2], 7+func.tp[2][8], 11+func.tp[2][8], 28+func.tp[2][7], 31+func.tp[2][7], 34+func.tp[2][7], 37+func.tp[2][7], 40+func.tp[2][7]};
	v4={17+func.tp[2][2], 8+func.tp[2][8], 12+func.tp[2][8], 29+func.tp[2][7], 32+func.tp[2][7], 35+func.tp[2][7], 38+func.tp[2][7], 41+func.tp[2][7]};
	Surface(v1, &func);
	Surface(v2, &func);
	Surface(v3, &func);
	Surface(v4, &func);
	v1={4+func.tp[1][2], 4+func.tp[1][8], 11+func.tp[1][8], 21+func.tp[1][8], 22+func.tp[1][8], 2+func.tp[1][7], 3+func.tp[1][7], 4+func.tp[1][7], 5+func.tp[1][7], 6+func.tp[1][7]};
	v2={5+func.tp[1][2], 5+func.tp[1][8], 12+func.tp[1][8], 22+func.tp[1][8], 23+func.tp[1][8], 17+func.tp[1][7], 18+func.tp[1][7], 19+func.tp[1][7], 20+func.tp[1][7]};
	for (i=0; i<3*div+4; i++) {
		v2.push_back(21+i+func.tp[1][7]);
	}
	v3={6+func.tp[1][2], 6+func.tp[1][8], 13+func.tp[1][8], 23+func.tp[1][8], 24+func.tp[1][8], 9+func.tp[1][7], 10+func.tp[1][7], 11+func.tp[1][7], 12+func.tp[1][7], 13+func.tp[1][7]};
	Volume(v1, &func);
	Volume(v2, &func);
	Volume(v3, &func);
	Tppush(9, 24, 16, 0, &func);

	//9
	func.tp_num = 9;
	Line(0.0, Wpml+Hsio2+Hsi, 0.0, 0.0, Wpml+Hsio2+Hsi, Wpml, &func);
	Copy(Lin, Wpml, 0.0, 0.0, "1", &func);
	Copy(Lin, Wm, 0.0, 0.0, "2", &func);
	Copy(Lin, Wr, 0.0, 0.0, "5", &func);
	Copy(Lin, g, 0.0, 0.0, "8", &func);
	Copy(Lin, Wh, 0.0, 0.0, "11", &func);
	Copy(Lin, Wm, 0.0, 0.0, "14", &func);
	Copy(Lin, Wpml, 0.0, 0.0, "17", &func);
	
	Line(0.0, Wpml+Hsio2+Hsi, L-Wpml, 0.0, Wpml+Hsio2+Hsi, L, &func);
	Copy(Lin, Wpml, 0.0, 0.0, "23", &func);
	Copy(Lin, Wm, 0.0, 0.0, "24", &func);
	Copy(Lin, Wr, 0.0, 0.0, "27", &func);
	Copy(Lin, g, 0.0, 0.0, "30", &func);
	Copy(Lin, Wh, 0.0, 0.0, "33", &func);
	Copy(Lin, Wm, 0.0, 0.0, "36", &func);
	Copy(Lin, Wpml, 0.0, 0.0, "39", &func);
	
	Copy(Lin, 0.0, 0.0, L-2*Wpml, "4 22", &func);
	
	Copy(Lin, 0.0, 0.0, Wm, "10", &func);
	Copy(Lin, 0.0, 0.0, L-2*Wpml-2*Wm, "49", &func);
	Copy(Lin, 0.0, 0.0, Wm, "52", &func);

	Copy(Poi, 0.0, 0.0, Wm, "10 12", &func);
	Copy(Poi, (Wh-Wpcm)/2, 0.0, 0.0, "37", &func);
	Copy(Poi, Wpcm, 0.0, 0.0, "39", &func);
	Copy(Poi, (Wh-Wpcm)/2, 0.0, 0.0, "40", &func);
	Surface("16 57 59 60 61 58", &func);
	Copy(Lin, 0.0, 0.0, dLc, "59 60 61", &func);
	for (i=0; i<div-1; i++) {
		v1={62+i*7, 63+i*7, 64+i*7};
		Copy(Lin, 0.0, 0.0, dLc, v1, &func);
	}
	v2={l1[0]+func.tp[2][9], l1[1]+func.tp[2][9], l1[2]+func.tp[2][9], l2[0]+func.tp[2][9], l2[1]+func.tp[2][9], 37+func.tp[2][9]};
	Copy(Poi, 0.0, 0.0, Wm, p1, &func);
	Surface(v2, &func);

	Surface("7 46 28 55 53 50", &func);
	v3={57+func.tp[2][9], 13+func.tp[2][9], 51+func.tp[2][9], 54+func.tp[2][9], 56+func.tp[2][9], 34+func.tp[2][9]};
	v4={58+func.tp[2][9], 19+func.tp[2][9], 47+func.tp[2][9], 40+func.tp[2][9]};
	for (i=0; i<div; i++) {
		v3.push_back(65+i*7+func.tp[2][9]);
		v4.push_back(66+i*7+func.tp[2][9]);
	}
	v3.push_back(65+div*7-3+func.tp[2][9]);
	v4.push_back(66+div*7-3+func.tp[2][9]);
	Surface(v3, &func);
	Surface(v4, &func);
	Tppush(0, 24+div*3, 63+div*7, 40+div*4, &func);

	//10
	func.tp_num = 7;
	v1={};
	for (i=0; i<3*div+24; i++) {
		v1.push_back(1+i);
	}
	Copy(Sur, 0.0, Hsi, 0.0, v1, &func);
	Tppush(24+div*3, 63+div*7, 40+div*4, 0, &func);

	//11
	func.tp_num = 9;
	v1 = {};
	for (i=0; i<div; i++) {
		v1.push_back(22+i*3);
	}
	Copy(Sur, 0.0, Hpcm, 0.0, v1, &func);
	Tppush(div, 3*div+1+div, 3*div+1+div*2+2, div*2+2, &func);

	//12
	func.tp_num = 4;
	Copy(Sur, 0.0, -(Wm+Hpcm), 0.0, "1 2 3 7 8 9", &func);
	v1={14+func.tp[2][4], 5+func.tp[2][12], 9+func.tp[2][12], 6+func.tp[2][9], 9+func.tp[2][9], 12+func.tp[2][9], 15+func.tp[2][9], 18+func.tp[2][9]};
	v2={15+func.tp[2][4], 6+func.tp[2][12], 10+func.tp[2][12], 7+func.tp[2][9], 10+func.tp[2][9], 13+func.tp[2][9], 16+func.tp[2][9], 19+func.tp[2][9]};
	v3={16+func.tp[2][4], 7+func.tp[2][12], 11+func.tp[2][12], 28+func.tp[2][9], 31+func.tp[2][9], 34+func.tp[2][9], 37+func.tp[2][9], 40+func.tp[2][9]};
	v4={17+func.tp[2][4], 8+func.tp[2][12], 12+func.tp[2][12], 29+func.tp[2][9], 32+func.tp[2][9], 35+func.tp[2][9], 38+func.tp[2][9], 41+func.tp[2][9]};
	Surface(v1, &func);
	Surface(v2, &func);
	Surface(v3, &func);
	Surface(v4, &func);
	v1={4+func.tp[1][4], 4+func.tp[1][12], 11+func.tp[1][12], 21+func.tp[1][12], 22+func.tp[1][12], 2+func.tp[1][9], 3+func.tp[1][9], 4+func.tp[1][9], 5+func.tp[1][9], 6+func.tp[1][9]};
	v2={5+func.tp[1][4], 5+func.tp[1][12], 12+func.tp[1][12], 22+func.tp[1][12], 23+func.tp[1][12], 17+func.tp[1][9], 18+func.tp[1][9], 19+func.tp[1][9], 20+func.tp[1][9]};
	for (i=0; i<div; i++) {
		v2.push_back(21+i*3+func.tp[1][9]);
		v2.push_back(23+i*3+func.tp[1][9]);
	}
	for (i=0; i<4; i++) {
		v2.push_back(21+div*3+i+func.tp[1][9]);
	}
	v2.push_back(1+func.tp[1][11]);
	for (i=0; i<div; i++) {
		v2.push_back(3+i*3+func.tp[1][11]);
		v2.push_back(4+i*3+func.tp[1][11]);
		v2.push_back(3*div+2+i+func.tp[1][11]);
	}
	v2.push_back(3*div-1+func.tp[1][11]);
	
	v3={6+func.tp[1][4], 6+func.tp[1][12], 13+func.tp[1][12], 23+func.tp[1][12], 24+func.tp[1][12], 9+func.tp[1][9], 10+func.tp[1][9], 11+func.tp[1][9], 12+func.tp[1][9], 13+func.tp[1][9]};
	
	Volume(v1, &func);
	Volume(v2, &func);
	Volume(v3, &func);
	Tppush(9, 24, 16, 0, &func);	
	
	Mat3D(1, "1 4 7 19 22 25", &func);
	v1={10, 13, 16, 28, 29, 31, 33, 34, 1+func.tp[0][12], 4+func.tp[0][12], 7+func.tp[0][12]};
	Mat3D(2, v1, &func);
	Mat3D(3, "30", &func);
	Mat3D(4, "32", &func);
	Mat3D(5, "3 6 9 21 24 27", &func);
	v1={12, 15, 18, 35, 36, 38, 40, 41, 3+func.tp[0][12], 6+func.tp[0][12], 9+func.tp[0][12]};
	Mat3D(6, v1, &func);
	Mat3D(7, "37", &func);
	Mat3D(8, "39", &func);
	v1={};
	for (i=0; i<3*div+2; i++){
		v1.push_back(20+i+func.tp[0][10]);
	}
	Mat3D(9, "44 45 46", &func);
	Mat3D(10, v1, &func);
	Mat3D(11, "2 5 8 20 23 26", &func);
	v1={11, 14, 17, 42, 43, 22+3*div+func.tp[0][10], 23+3*div+func.tp[0][10], 24+3*div+func.tp[0][10], 2+func.tp[0][12], 5+func.tp[0][12], 8+func.tp[0][12]};
	Mat3D(12, v1, &func);
	for (i=0; i<div; i++){
		v1={i+1+func.tp[0][11]};
		Mat3D(13+i, v1, &func);
	}
	v1={3+func.tp[1][7], 4+func.tp[1][7], 5+func.tp[1][7], 10+func.tp[1][7], 11+func.tp[1][7], 12+func.tp[1][7], 17+func.tp[1][7], 18+func.tp[1][7], 19+func.tp[1][7], 20+func.tp[1][7]};
	for (i=0; i<3*div+1; i++) {
		v1.push_back(21+i+func.tp[1][7]);
	}
	v1.push_back(23+3*div+func.tp[1][7]);
	v2={5+func.tp[1][10], 8+func.tp[1][10], 9+func.tp[1][10], 10+func.tp[1][10], 11+func.tp[1][10], 12+func.tp[1][10], 13+func.tp[1][10], 14+func.tp[1][10], 15+func.tp[1][10], 16+func.tp[1][10], 27+func.tp[1][10], 30+func.tp[1][10], 31+func.tp[1][10], 32+func.tp[1][10], 33+func.tp[1][10], 34+func.tp[1][10], 35+func.tp[1][10], 36+func.tp[1][10], 37+func.tp[1][10], 38+func.tp[1][10], 49+func.tp[1][10], 50+func.tp[1][10], 51+func.tp[1][10], 52+func.tp[1][10], 53+func.tp[1][10], 54+func.tp[1][10], 55+func.tp[1][10], 56+func.tp[1][10], 57+func.tp[1][10], 58+func.tp[1][10]};
	for (i=0; i<7*div+5; i++) {
		v2.push_back(59+i+func.tp[1][10]);
	}
	v3={3+func.tp[1][9], 4+func.tp[1][9], 5+func.tp[1][9], 10+func.tp[1][9], 11+func.tp[1][9], 12+func.tp[1][9], 17+func.tp[1][9], 18+func.tp[1][9], 19+func.tp[1][9], 20+func.tp[1][9]};
	for (i=0; i<3*div+1; i++) {
		v3.push_back(21+i+func.tp[1][9]);
	}
	v3.push_back(23+3*div+func.tp[1][9]);
	v4={};
	for (i=0; i<3*div+1+div; i++) {
		v4.push_back(1+i+func.tp[1][11]);
	}
	copy(v2.begin(),v2.end(), back_inserter(v1));
	copy(v3.begin(),v3.end(), back_inserter(v1));
	copy(v4.begin(),v4.end(), back_inserter(v1));
	Unstr(Sur, unstr, v1, &func);
	
	Trans_Gene(trans, gene, &func);
	Printvv(func.tp);

	Filename(lpsname, &func);

	ofs.close();

	return 0;
}