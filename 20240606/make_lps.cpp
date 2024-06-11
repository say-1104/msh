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
		Copy(Sur, 0.0, Y, 0.0, 			"1 2 3 4 5 6 7 8 9", &func);
		Tppush(9, 24, 16, 0, &func);
	};

	//1, 2
	lpsxy(1, 0.0);
	lpsxy(2, Wpml);

	//3
	lpsz(1, Wpml);
	
	//4, 5
	lpsxy(4, H-Wpml);
	lpsxy(5, H);

	//6
	lpsz(4, Wpml);
	
	//7
	func.tp_num = 7;
	vector <int> l1, l2;
	Line(0.0, Wpml+Hsio2, 0.0, 		0.0, Wpml+Hsio2, Wpml, 		&func);
	Copy(Lin, Wpml, 0.0, 0.0, 		"1", &func);
	Copy(Lin, Wx+Wr+g, 0.0, 0.0, 	"2", &func);
	Copy(Lin, Wh, 0.0, 0.0, 		"5", &func);
	Copy(Lin, Wx, 0.0, 0.0, 		"8", &func);
	Copy(Lin, Wpml, 0.0, 0.0, 		"11", &func);
	
	Line(0.0, Wpml+Hsio2, L-Wpml, 	0.0, Wpml+Hsio2, L, 		&func);
	Copy(Lin, Wpml, 0.0, 0.0, 		"17", &func);
	Copy(Lin, Wx, 0.0, 0.0, 		"18", &func);
	Copy(Lin, Wr, 0.0, 0.0, 		"21", &func);
	Copy(Lin, g, 0.0, 0.0, 			"24", &func);
	Copy(Lin, Wh, 0.0, 0.0, 		"27", &func);
	Copy(Lin, Wx, 0.0, 0.0, 		"30", &func);
	Copy(Lin, Wpml, 0.0, 0.0, 		"33", &func);
	
	Copy(Lin, 0.0, 0.0, L-2*Wpml, "4 16", &func);
	
	Line(Wpml+Wx, Wpml+Hsio2, Wpml+Wz, 	Wpml+Wx+Wr, Wpml+Hsio2, Wpml+Wz, 		&func);
	Copy(Poi, g, 0.0, 0.0, 			"30", &func);
	Copy(Poi, Wh, 0.0, 0.0, 		"31", &func);

	Copy(Poi, 0.0, 0.0, Wz, "6 8", &func);
	Surface("10 46 45 47", &func);
	Copy(Lin, 0.0, 0.0, dLc, "43 44 45", &func);
	for(i=0; i<div-1; i++) {
		v1={48+i*7, 49+i*7, 50+i*7};
		Copy(Lin, 0.0, 0.0, dLc, v1, &func);
	}
	Copy(Lin, 0.0, 0.0, -Wz, "25 28 31", &func);

	l1={22, 40, 7, 46, 44, 43};
	l2={34, 41, 13, 47};
	l1=Plusv(l1, func.tp[2][7]);
	l2=Plusv(l2, func.tp[2][7]);
	for(i=0; i<div; i++){
		l1.push_back(51+i*7+func.tp[2][7]);
		l2.push_back(54+i*7+func.tp[2][7]);
	}
	l1.push_back(7*div+48+func.tp[2][7]);
	l2.push_back(7*div+51+func.tp[2][7]);
	Surface(l1, &func);
	Surface(l2, &func);
	Tppush(0, 3*div+20, 7*div+51, 4*div+32, &func);
	
	//8
	func.tp_num = 2;
	Copy(Sur, 0.0, Hsio2, 0.0, "1 2 3 7 8 9", &func);
	v1={14+func.tp[2][2], 5+func.tp[2][8], 9+func.tp[2][8], 6+func.tp[2][7], 9+func.tp[2][7], 12+func.tp[2][7]};
	v2={15+func.tp[2][2], 6+func.tp[2][8], 10+func.tp[2][8], 7+func.tp[2][7], 10+func.tp[2][7], 13+func.tp[2][7]};
	v3={16+func.tp[2][2], 7+func.tp[2][8], 11+func.tp[2][8], 22+func.tp[2][7], 25+func.tp[2][7], 28+func.tp[2][7], 31+func.tp[2][7], 34+func.tp[2][7]};
	v4={17+func.tp[2][2], 8+func.tp[2][8], 12+func.tp[2][8], 23+func.tp[2][7], 26+func.tp[2][7], 29+func.tp[2][7], 32+func.tp[2][7], 35+func.tp[2][7]};
	Surface(v1, &func);
	Surface(v2, &func);
	Surface(v3, &func);
	Surface(v4, &func);
	v1={4+func.tp[1][2], 4+func.tp[1][8], 11+func.tp[1][8], 21+func.tp[1][8], 22+func.tp[1][8], 2+func.tp[1][7], 3+func.tp[1][7], 4+func.tp[1][7]};
	v2={5+func.tp[1][2], 5+func.tp[1][8], 12+func.tp[1][8], 22+func.tp[1][8], 23+func.tp[1][8], 15+func.tp[1][7]};
	for (i=0; i<3*div+5; i++) {
		v2.push_back(16+i+func.tp[1][7]);
	}
	v3={6+func.tp[1][2], 6+func.tp[1][8], 13+func.tp[1][8], 23+func.tp[1][8], 24+func.tp[1][8], 7+func.tp[1][7], 8+func.tp[1][7], 9+func.tp[1][7], 10+func.tp[1][7], 11+func.tp[1][7]};
	Volume(v1, &func);
	Volume(v2, &func);
	Volume(v3, &func);
	Tppush(9, 24, 16, 0, &func);

	//9
	func.tp_num = 9;
	Line(0.0, Wpml+Hsio2+Hsi, 0.0, 		0.0, Wpml+Hsio2+Hsi, Wpml, 		&func);
	Copy(Lin, Wpml, 0.0, 0.0, 		"1", &func);
	Copy(Lin, Wx+Wr+g, 0.0, 0.0, 	"2", &func);
	Copy(Lin, Wh, 0.0, 0.0, 		"5", &func);
	Copy(Lin, Wx, 0.0, 0.0, 		"8", &func);
	Copy(Lin, Wpml, 0.0, 0.0, 		"11", &func);
	
	Line(0.0, Wpml+Hsio2+Hsi, L-Wpml, 	0.0, Wpml+Hsio2+Hsi, L, 		&func);
	Copy(Lin, Wpml, 0.0, 0.0, 		"17", &func);
	Copy(Lin, Wx, 0.0, 0.0, 		"18", &func);
	Copy(Lin, Wr, 0.0, 0.0, 		"21", &func);
	Copy(Lin, g, 0.0, 0.0, 			"24", &func);
	Copy(Lin, Wh, 0.0, 0.0, 		"27", &func);
	Copy(Lin, Wx, 0.0, 0.0, 		"30", &func);
	Copy(Lin, Wpml, 0.0, 0.0, 		"33", &func);
	
	Copy(Lin, 0.0, 0.0, L-2*Wpml, "4 16", &func);
	
	Line(Wpml+Wx, Wpml+Hsio2+Hsi, Wpml+Wz, 	Wpml+Wx+Wr, Wpml+Hsio2+Hsi, Wpml+Wz, 		&func);
	Copy(Poi, g, 0.0, 0.0, 			"30", &func);
	Copy(Poi, (Wh-Wpcm)/2, 0.0, 0.0, 		"31", &func);
	Copy(Poi, Wpcm, 0.0, 0.0, 		"32", &func);
	Copy(Poi, (Wh-Wpcm)/2, 0.0, 0.0, 		"33", &func);

	Copy(Poi, 0.0, 0.0, Wz, "6 8", &func);
	Surface("10 48 45 46 47 49", &func);
	Copy(Lin, 0.0, 0.0, dLc, "43 44 45 46 47", &func);
	for(i=0; i<div-1; i++) {
		v1={50+i*11, 51+i*11, 52+i*11, 53+i*11, 54+i*11};
		Copy(Lin, 0.0, 0.0, dLc, v1, &func);
	}
	Copy(Lin, 0.0, 0.0, -Wz, "25 28", &func);
	Copy(Poi, 0.0, 0.0, -Wz, "23", &func);
	v1={31, 11*div+52, 11*div+41, 11*div+42, 11*div+43, 11*div+53};
	v1=Plusv(v1, func.tp[2][9]);
	Surface(v1, &func);

	l1={22, 40, 7, 48, 44, 43};
	l2={34, 41, 13, 49};
	l1=Plusv(l1, func.tp[2][9]);
	l2=Plusv(l2, func.tp[2][9]);
	for(i=0; i<div; i++){
		l1.push_back(55+i*11+func.tp[2][9]);
		l2.push_back(60+i*11+func.tp[2][9]);
	}
	l1.push_back(11*div+50+func.tp[2][9]);
	l2.push_back(11*div+53+func.tp[2][9]);
	Surface(l1, &func);
	Surface(l2, &func);
	Tppush(0, 5*div+20, 11*div+53, 6*div+34, &func);

	//10
	func.tp_num = 7;
	v1={1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
	for (i=0; i<div+1; i++) {
		v1.push_back(16+3*i);
		v1.push_back(17+3*i);
	}
	v1.push_back(19+3*div);
	v1.push_back(20+3*div);
	Copy(Sur, 0.0, Hsi, 0.0, v1, &func);

	v1={45+func.tp[2][7], 31+func.tp[2][10], 32+func.tp[2][10], 45+func.tp[2][9], 46+func.tp[2][9], 47+func.tp[2][9]};
	Surface(v1, &func);
	for (i=0; i<div; i++) {
		v1={50+7*i+func.tp[2][7], 35+4*i+func.tp[2][10], 36+4*i+func.tp[2][10], 52+11*i+func.tp[2][9], 53+11*i+func.tp[2][9], 54+11*i+func.tp[2][9]};
		Surface(v1, &func);
	}

	v1={15+func.tp[1][7], 10+func.tp[1][10], 45+func.tp[1][10], 46+func.tp[1][10], 51+6*div+func.tp[1][10], 15+func.tp[1][9]};
	Volume(v1, &func);
	for (i=0; i<div; i++) {
		v1={18+3*i+func.tp[1][7], 51+6*i+func.tp[1][10], 52+6*i+func.tp[1][10], 51+6*div+i+func.tp[1][10], 52+6*div+i+func.tp[1][10], 18+5*i+func.tp[1][9], 19+5*i+func.tp[1][9], 20+5*i+func.tp[1][9]};
		Volume(v1, &func);
	}
	v1={3*div+18+func.tp[1][7], 7*div+51+func.tp[1][10], 49+6*div+func.tp[1][10], 50+6*div+func.tp[1][10], 31+func.tp[1][10], 18+5*div+func.tp[1][9]};
	Volume(v1, &func);
	Tppush(3*div+20, 7*div+51, 4*div+32, 0, &func);

	//11
	/*func.tp_num = 9;
	v1 = {};
	for (i=0; i<div; i++) {
		v1.push_back(19+i*3);
	}
	Copy(Sur, 0.0, Hpcm, 0.0, v1, &func);
	Tppush(div, 3*div+1+div, 3*div+1+div*2+2, div*2+2, &func);

	//12
	/*func.tp_num = 4;
	Copy(Sur, 0.0, -(Wy+Hpcm), 0.0, "1 2 3 7 8 9", &func);
	v1={14+func.tp[2][4], 5+func.tp[2][12], 9+func.tp[2][12], 6+func.tp[2][9], 9+func.tp[2][9], 12+func.tp[2][9]};
	v2={15+func.tp[2][4], 6+func.tp[2][12], 10+func.tp[2][12], 7+func.tp[2][9], 10+func.tp[2][9], 13+func.tp[2][9]};
	v3={16+func.tp[2][4], 7+func.tp[2][12], 11+func.tp[2][12], 22+func.tp[2][9], 25+func.tp[2][9], 28+func.tp[2][9], 31+func.tp[2][9], 34+func.tp[2][9]};
	v4={17+func.tp[2][4], 8+func.tp[2][12], 12+func.tp[2][12], 23+func.tp[2][9], 26+func.tp[2][9], 29+func.tp[2][9], 32+func.tp[2][9], 35+func.tp[2][9]};
	Surface(v1, &func);
	Surface(v2, &func);
	Surface(v3, &func);
	Surface(v4, &func);
	v1={4+func.tp[1][4], 4+func.tp[1][12], 11+func.tp[1][12], 21+func.tp[1][12], 22+func.tp[1][12], 2+func.tp[1][9], 3+func.tp[1][9], 4+func.tp[1][9]};
	v2={5+func.tp[1][4], 5+func.tp[1][12], 12+func.tp[1][12], 22+func.tp[1][12], 23+func.tp[1][12], 15+func.tp[1][9]};
	for (i=0; i<div; i++) {
		v2.push_back(16+i*5+func.tp[1][9]);
		v2.push_back(17+i*5+func.tp[1][9]);
		v2.push_back(18+i*5+func.tp[1][9]);
		v2.push_back(20+i*5+func.tp[1][9]);
	}
	for (i=0; i<5; i++) {
		v2.push_back(16+div*5+i+func.tp[1][9]);
	}
	v2.push_back(1+func.tp[1][11]);
	for (i=0; i<div; i++) {
		v2.push_back(3+i*3+func.tp[1][11]);
		v2.push_back(4+i*3+func.tp[1][11]);
		v2.push_back(3*div+2+i+func.tp[1][11]);
	}
	v2.push_back(3*div-1+func.tp[1][11]);
	
	v3={6+func.tp[1][4], 6+func.tp[1][12], 13+func.tp[1][12], 23+func.tp[1][12], 24+func.tp[1][12], 7+func.tp[1][9], 8+func.tp[1][9], 9+func.tp[1][9], 10+func.tp[1][9], 11+func.tp[1][9]};
	
	Volume(v1, &func);
	Volume(v2, &func);
	Volume(v3, &func);
	Tppush(9, 24, 16, 0, &func);	
	
	//Material
	/*Mat3D(1, "1 4 7 19 22 25", &func);
	v1={10, 13, 16, 28, 29, 31, 32, 1+func.tp[0][12], 4+func.tp[0][12], 7+func.tp[0][12]};
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
	
	Trans_Gene(trans, gene, &func);*/
	Printvv(func.tp);

	Filename(lpsname, &func);

	ofs.close();

	return 0;
}