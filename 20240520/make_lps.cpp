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
	double Leff = atof(argv[7]);	//cPCM層の長さ

	//DCのパラメータ(固定値)
	double Hsi = 0.22;				//Si層の厚さ
	double Hsio2 = 2.0;				//SiO2層の厚さ
	double Wm = 2.0;				//導波路とPML間の距離
	double Wpml = 1.0;				//PML幅
	double dLc = 1.4;
	int div = Div(Lc, dLc);

	double L = 2 * Wpml + 2 * Wm + Lc;
	double W = 2 * Wpml + 2 * Wm + Wr + g + Wh;
	double H = 2 * Wpml + Hsio2 + Hsi + Hpcm + Wm;



	//メッシュのパラメータ
	double unstr = 0.04, trans = 0.4, gene = 0.4;
	string lpsname = "pcmDC";

	cout << "writing " << lpsname << "..." << endl;

	Tppush(0, 0, 0, 0, &func);

	vector<int> v1, v2, v3, v4, v5, v6;

	//1
	func.tp_num = 1;
	Line(0.0, 0.0, 0.0, 	0.0, 0.0, Wpml, &func);
	Line(0.0, 0.0, Wpml, 	0.0, 0.0, L-Wpml, &func);
	Line(0.0, 0.0, L-Wpml, 	0.0, 0.0, L, &func);
	Copy(Lin, Wpml, 0.0, 0.0, 		"1 2 3", &func);
	Copy(Lin, W-2*Wpml, 0.0, 0.0, 	"4 5 6", &func);
	Copy(Lin, Wpml, 0.0, 0.0, 		"11 12 13", &func);
	Tppush(0, 9, 24, 16, &func);

	//2
	func.tp_num = 2;
	Line(0.0, Wpml, 0.0, 	0.0, Wpml, Wpml, &func);
	Line(0.0, Wpml, Wpml, 	0.0, Wpml, L-Wpml, &func);
	Line(0.0, Wpml, L-Wpml, 0.0, Wpml, L, &func);
	Copy(Lin, Wpml, 0.0, 0.0, 		"1 2 3", &func);
	Copy(Lin, W-2*Wpml, 0.0, 0.0, 	"4 5 6", &func);
	Copy(Lin, Wpml, 0.0, 0.0, 		"11 12 13", &func);
	Tppush(0, 9, 24, 16, &func);

	//3
	func.tp_num = 1;
	Copy(Sur, 0.0, Wpml, 0.0, 		"1 2 3 4 5 6 7 8 9", &func);
	Tppush(9, 24, 16, 0, &func);
	
	//4
	func.tp_num = 4;
	Line(0.0, H-Wpml, 0.0, 		0.0, H-Wpml, Wpml, &func);
	Line(0.0, H-Wpml, Wpml, 	0.0, H-Wpml, L-Wpml, &func);
	Line(0.0, H-Wpml, L-Wpml, 	0.0, H-Wpml, L, &func);
	Copy(Lin, Wpml, 0.0, 0.0, 		"1 2 3", &func);
	Copy(Lin, W-2*Wpml, 0.0, 0.0, 	"4 5 6", &func);
	Copy(Lin, Wpml, 0.0, 0.0, 		"11 12 13", &func);
	Tppush(0, 9, 24, 16, &func);

	//5
	func.tp_num = 5;
	Line(0.0, H, 0.0, 		0.0, H, Wpml, &func);
	Line(0.0, H, Wpml, 		0.0, H, L-Wpml, &func);
	Line(0.0, H, L-Wpml, 	0.0, H, L, &func);
	Copy(Lin, Wpml, 0.0, 0.0, 		"1 2 3", &func);
	Copy(Lin, W-2*Wpml, 0.0, 0.0, 	"4 5 6", &func);
	Copy(Lin, Wpml, 0.0, 0.0, 		"11 12 13", &func);
	Tppush(0, 9, 24, 16, &func);

	//6
	func.tp_num = 4;
	Copy(Sur, 0.0, Wpml, 0.0, "1 2 3 4 5 6 7 8 9", &func);
	Tppush(9, 24, 16, 0, &func);
	
	//7
	func.tp_num = 7;
	Line(0.0, Wpml+Hsio2, 0.0, 0.0, Wpml+Hsio2, Wpml, &func);
	Copy(Lin, Wpml, 0.0, 0.0, "1", &func);
	Copy(Lin, Wm, 0.0, 0.0, "2", &func);
	Copy(Lin, Wr, 0.0, 0.0, "5", &func);
	Copy(Lin, g, 0.0, 0.0, "8", &func);
	Copy(Lin, Wh, 0.0, 0.0, "11", &func);
	Copy(Lin, Wm, 0.0, 0.0, "14", &func);
	Copy(Lin, Wpml, 0.0, 0.0, "17", &func);
	
	Line(0.0, Wpml+Hsio2, L-Wpml, 0.0, Wpml+Hsio2, L, &func);
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
	Copy(Poi, Wh, 0.0, 0.0, "37", &func);
	Surface("16 57 59 58", &func);
	Copy(Lin, 0.0, 0.0, L-2*Wpml-2*Wm, "59", &func);
	Copy(Lin, 0.0, 0.0, Wm, "60", &func);

	Surface("7 46 28 55 53 50", &func);
	Surface("13 51 54 56 34 63 61 57", &func);
	Surface("19 58 62 64 40 47", &func);
	Tppush(0, 25, 64, 40, &func);
	
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
	v2={5+func.tp[1][2], 5+func.tp[1][8], 12+func.tp[1][8], 22+func.tp[1][8], 23+func.tp[1][8], 17+func.tp[1][7], 18+func.tp[1][7], 19+func.tp[1][7], 20+func.tp[1][7], 21+func.tp[1][7], 22+func.tp[1][7], 23+func.tp[1][7], 24+func.tp[1][7], 25+func.tp[1][7]};
	v3={6+func.tp[1][2], 6+func.tp[1][8], 13+func.tp[1][8], 23+func.tp[1][8], 24+func.tp[1][8], 9+func.tp[1][7], 10+func.tp[1][7], 11+func.tp[1][7], 12+func.tp[1][7], 13+func.tp[1][7]};
	Volume(v1, &func);
	Volume(v2, &func);
	Volume(v3, &func);
	Tppush(9, 24, 16, 0, &func);

	//9
	vector <int> p1, l1, l2, l3, l4;
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
	p1={37+4*div, 38+4*div};
	l1={62+(div-1)*7, 63+(div-1)*7, 64+(div-1)*7};
	l2={l1[2]+5, l1[2]+6};
	v2={l1[0]+func.tp[2][9], l1[1]+func.tp[2][9], l1[2]+func.tp[2][9], l2[0]+func.tp[2][9], l2[1]+func.tp[2][9], 37+func.tp[2][9]};
	Copy(Poi, 0.0, 0.0, Wm, p1, &func);
	Surface(v2, &func);

	Surface("7 46 28 55 53 50", &func);
	v3={57+func.tp[2][9], 13+func.tp[2][9], 51+func.tp[2][9], 54+func.tp[2][9], 56+func.tp[2][9], 34+func.tp[2][9]};
	v4={58+func.tp[2][9], 19+func.tp[2][9], 47+func.tp[2][9], 40+func.tp[2][9]};
	for (i=0; i<div; i++) {
		v3.push_back(65+i*7+func.tp[2][9]);
		l3.push_back(65+i*7+func.tp[2][9]);
		v4.push_back(66+i*7+func.tp[2][9]);
		l4.push_back(66+i*7+func.tp[2][9]);
	}
	v3.push_back(65+div*7-3+func.tp[2][9]);
	v4.push_back(66+div*7-3+func.tp[2][9]);
	Surface(v3, &func);
	Surface(v4, &func);
	Tppush(0, 24+div*3, 63+div*7, 40+div*4, &func);

	//10
	func.tp_num = 7;
	Copy(Sur, 0.0, Hsi, 0.0, "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19", &func);
	Copy(Lin, 0.0, Hsi, 0.0, "57 58 63 64", &func);
	v1={59+func.tp[2][7], 37+func.tp[2][10], 38+func.tp[2][10], 59+func.tp[2][9], 60+func.tp[2][9], 61+func.tp[2][9]};
	v2={60+func.tp[2][7], 39+func.tp[2][10], 40+func.tp[2][10], l1[0]+func.tp[2][9], l1[1]+func.tp[2][9], l1[2]+func.tp[2][9]};
	v3={61+func.tp[2][7], 37+func.tp[2][10], 39+func.tp[2][10]};
	v4={62+func.tp[2][7], 38+func.tp[2][10], 40+func.tp[2][10]};
	for (i=0; i<div; i++) {
		v3.push_back(l3[i]);
		v4.push_back(l4[i]);
	}
	Surface(v1, &func);
	Surface(v2, &func);
	Surface(v3, &func);
	Surface(v4, &func);
	
	v1={20+func.tp[1][7], 16+func.tp[1][10], 57+func.tp[1][10], 58+func.tp[1][10], 61+func.tp[1][10], 20+func.tp[1][9]};
	v2={21+func.tp[1][7], 61+func.tp[1][10], 62+func.tp[1][10], 63+func.tp[1][10], 64+func.tp[1][10]};
	for (i=0; i<3*div; i++) {
		v2.push_back(21+i+func.tp[1][9]);
	}
	v3={22+func.tp[1][7], 37+func.tp[1][10], 59+func.tp[1][10], 60+func.tp[1][10], 62+func.tp[1][10], 21+div*3+func.tp[1][9]};
	v4={24+func.tp[1][7], 51+func.tp[1][10], 54+func.tp[1][10], 56+func.tp[1][10], 34+func.tp[1][10], 57+func.tp[1][10], 63+func.tp[1][10], 59+func.tp[1][10], 13+func.tp[1][10], 23+div*3+func.tp[1][9]};
	v5={25+func.tp[1][7], 58+func.tp[1][10], 64+func.tp[1][10], 60+func.tp[1][10], 40+func.tp[1][10], 47+func.tp[1][10], 19+func.tp[1][10], 24+div*3+func.tp[1][9]};
	Volume(v1, &func);
	Volume(v2, &func);
	Volume(v3, &func);
	Copy(Sur, 0.0, Hsi, 0.0, "23", &func);
	Volume(v4, &func);		
	Volume(v5, &func);
	Tppush(25, 64, 40, 0, &func);

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
	v2={5+func.tp[1][4], 5+func.tp[1][12], 12+func.tp[1][12], 22+func.tp[1][12], 23+func.tp[1][12], 17+func.tp[1][8], 18+func.tp[1][8], 19+func.tp[1][8], 20+func.tp[1][8]};
	for (i=0; i<div; i++) {
		v2.push_back(21+i*3+func.tp[1][9]);
		v2.push_back(23+i*3+func.tp[1][9]);
	}
	for (i=0; i<4; i++) {
		v2.push_back(21+div*3+i+func.tp[1][9]);
	}
	v3={6+func.tp[1][4], 6+func.tp[1][12], 13+func.tp[1][12], 23+func.tp[1][12], 24+func.tp[1][12], 9+func.tp[1][9], 10+func.tp[1][9], 11+func.tp[1][9], 12+func.tp[1][9], 13+func.tp[1][9]};
	
	Volume(v1, &func);
	Volume(v2, &func);
	Volume(v3, &func);
	Tppush(9, 24, 16, 0, &func);	
	
	/*Mat3D(1, "1 4 7 19 22 25", &func);
	Mat3D(2, "10 13 16 28 29 31 32 51 54 57", &func);
	Mat3D(3, "30", &func);
	Mat3D(4, "3 6 9 21 24 26", &func);
	Mat3D(5, "12 15 18 33 34 36 38 39 53 56 58", &func);
	Mat3D(6, "35", &func);
	Mat3D(7, "37", &func);
	Mat3D(8, "42", &func);
	Mat3D(9, "44", &func);
	Mat3D(10, "45", &func);
	Mat3D(11, "43 46", &func);
	Mat3D(12, "2 5 8 20 23 27", &func);
	Mat3D(13, "11 14 17 40 41 47 48 52 55 59", &func);
	if(Leff == 0 || Leff == Lc) Mat3D(14, "49 50", &func);
	else {
		Mat3D(14, "49", &func);
		Mat3D(15, "50", &func);
	};
	v1={3+func.tp[1][6], 8+func.tp[1][6], 10+func.tp[1][6], 15+func.tp[1][6], 16+func.tp[1][6], 17+func.tp[1][6], 18+func.tp[1][6], 19+func.tp[1][6]};
	v2={5+func.tp[1][9], 8+func.tp[1][9], 9+func.tp[1][9], 10+func.tp[1][9], 21+func.tp[1][9], 24+func.tp[1][9], 25+func.tp[1][9], 26+func.tp[1][9], 43+func.tp[1][9], 44+func.tp[1][9], 45+func.tp[1][9], 46+func.tp[1][9], 47+func.tp[1][9], 48+func.tp[1][9], 49+func.tp[1][9], 50+func.tp[1][9], 51+func.tp[1][9], 52+func.tp[1][9], 53+func.tp[1][9], 54+func.tp[1][9], 55+func.tp[1][9], 56+func.tp[1][9], 27+func.tp[1][9], 30+func.tp[1][9], 31+func.tp[1][9], 32+func.tp[1][9]};
	v3={3+func.tp[1][8], 8+func.tp[1][8], 10+func.tp[1][8], 15+func.tp[1][8], 16+func.tp[1][8], 17+func.tp[1][8], 18+func.tp[1][8], 19+func.tp[1][8], 20+func.tp[1][8], 21+func.tp[1][8], 22+func.tp[1][8], 23+func.tp[1][8], 24+func.tp[1][8]};
	v4={1+func.tp[1][11], 3+func.tp[1][11], 4+func.tp[1][11], 5+func.tp[1][11], 6+func.tp[1][11], 7+func.tp[1][11], 1+func.tp[1][10], 2+func.tp[1][10]};
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