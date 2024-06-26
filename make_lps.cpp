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

int main(int argc, char *argv[]){
	Function func;
	func.flag = 0;
	func.tp_num = 0;
	
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

	double Lcur = 20;				//曲げ導波路部分の全長
	double Wcur = 8;				//曲げ導波路部分の幅
	Curv curv;
	curv = calcCurv(Lcur, Wcur);

	double L = 2 * Wpml + Wm + Lc + Lcur;
	double W = 2 * Wpml + 2 * Wm + Wr + g + Wh + Wcur;
	double H = 2 * Wpml + Hsio2 + Hsi + Hpcm + Wm;

	//メッシュのパラメータ
	double unstr = 0.04, trans = 0.5, gene = 0.5;
	string lpsname = "pcmDC";

	cout << "writing " << lpsname << "..." << endl;

	//1
	Line(0.0, 0.0, 0.0, 	0.0, 0.0, Wpml, &func);
	Line(0.0, 0.0, Wpml, 	0.0, 0.0, L-Wpml, &func);
	Line(0.0, 0.0, L-Wpml, 	0.0, 0.0, L, &func);
	Copy(Lin, Wpml, 0.0, 0.0, 		"1 2 3", &func);
	Copy(Lin, W-2*Wpml, 0.0, 0.0, 	"4 5 6", &func);
	Copy(Lin, Wpml, 0.0, 0.0, 		"11 12 13", &func);
	Tppush(0, 9, 24, 16, &func);

	//2
	func.tp_num = 1;
	Line(0.0, Wpml, 0.0, 	0.0, Wpml, Wpml, &func);
	Line(0.0, Wpml, Wpml, 	0.0, Wpml, L-Wpml, &func);
	Line(0.0, Wpml, L-Wpml, 0.0, Wpml, L, &func);
	Copy(Lin, Wpml, 0.0, 0.0, 		"1 2 3", &func);
	Copy(Lin, W-2*Wpml, 0.0, 0.0, 	"4 5 6", &func);
	Copy(Lin, Wpml, 0.0, 0.0, 		"11 12 13", &func);
	Tppush(0, 9, 24, 16, &func);

	//3
	func.tp_num = 0;
	Copy(Sur, 0.0, Wpml, 0.0, 		"1 2 3 4 5 6 7 8 9", &func);
	Tppush(9, 24, 16, 0, &func);
	
	//4
	func.tp_num = 3;
	Line(0.0, H-Wpml, 0.0, 		0.0, H-Wpml, Wpml, &func);
	Line(0.0, H-Wpml, Wpml, 	0.0, H-Wpml, L-Wpml, &func);
	Line(0.0, H-Wpml, L-Wpml, 	0.0, H-Wpml, L, &func);
	Copy(Lin, Wpml, 0.0, 0.0, 		"1 2 3", &func);
	Copy(Lin, W-2*Wpml, 0.0, 0.0, 	"4 5 6", &func);
	Copy(Lin, Wpml, 0.0, 0.0, 		"11 12 13", &func);
	Tppush(0, 9, 24, 16, &func);

	//5
	func.tp_num = 4;
	Line(0.0, H, 0.0, 		0.0, H, Wpml, &func);
	Line(0.0, H, Wpml, 		0.0, H, L-Wpml, &func);
	Line(0.0, H, L-Wpml, 	0.0, H, L, &func);
	Copy(Lin, Wpml, 0.0, 0.0, 		"1 2 3", &func);
	Copy(Lin, W-2*Wpml, 0.0, 0.0, 	"4 5 6", &func);
	Copy(Lin, Wpml, 0.0, 0.0, 		"11 12 13", &func);
	Tppush(0, 9, 24, 16, &func);

	//6
	func.tp_num = 3;
	Copy(Sur, 0.0, Wpml, 0.0, "1 2 3 4 5 6 7 8 9", &func);
	Tppush(9, 24, 16, 0, &func);
	
	//7
	func.tp_num = 6;
	Line(0.0, Wpml+Hsio2, 0.0, 0.0, Wpml+Hsio2, Wpml, &func);
	Copy(Lin, Wpml, 0.0, 0.0, "1", &func);
	Copy(Lin, Wm+Wr+g, 0.0, 0.0, "2", &func);
	Copy(Lin, Wh, 0.0, 0.0, "5", &func);
	Copy(Lin, Wm+Wcur, 0.0, 0.0, "8", &func);
	Copy(Lin, Wpml, 0.0, 0.0, "11", &func);
	
	Line(0.0, Wpml+Hsio2, L-Wpml, 0.0, Wpml+Hsio2, L, &func);
	Copy(Lin, Wpml, 0.0, 0.0, "17", &func);
	Copy(Lin, Wm, 0.0, 0.0, "18", &func);
	Copy(Lin, Wr, 0.0, 0.0, "21", &func);
	Copy(Lin, g+Wcur, 0.0, 0.0, "24", &func);
	Copy(Lin, Wh, 0.0, 0.0, "27", &func);
	Copy(Lin, Wm, 0.0, 0.0, "30", &func);
	Copy(Lin, Wpml, 0.0, 0.0, "33", &func);
	
	Copy(Lin, 0.0, 0.0, L-2*Wpml, "4 16", &func);
	
	Copy(Lin, 0.0, 0.0, -(L-2*Wpml-Wm), "25", &func);

	Copy(Poi, 0.0, 0.0, Wm, "6 8", &func);
	Copy(Poi, Wh, 0.0, 0.0, "31", &func);
	Surface("10 46 48 47", &func);
	Copy(Lin, 0.0, 0.0, Lc, "48", &func);
	Rotate(Lin, Y, Wpml+Wm+Wr+g+Wh/2+Wcur-curv.R, Wpml+Hsio2, L-Wpml, curv.angle, "31", &func);
	Rotate(Lin, Y, Wpml+Wm+Wr+g+Wh/2+curv.R, Wpml+Hsio2, Wpml+Wm+Lc, curv.angle, "49", &func);
	

	Surface("7 40 22 44 43 45 28 53 55 50 46", &func);
	Surface("13 41 34 54 56 51 47", &func);
	Tppush(0, 21, 56, 36, &func);

	//8
	func.tp_num = 1;
	Copy(Sur, 0.0, Hsio2, 0.0, "1 2 3 7 8 9", &func);
	vector<int> v1, v2, v3, v4, v5, v6;
	v1={14+func.tp[2][1], 5+func.tp[2][7], 9+func.tp[2][7], 6+func.tp[2][6], 9+func.tp[2][6], 12+func.tp[2][6]};
	v2={15+func.tp[2][1], 6+func.tp[2][7], 10+func.tp[2][7], 7+func.tp[2][6], 10+func.tp[2][6], 13+func.tp[2][6]};
	v3={16+func.tp[2][1], 7+func.tp[2][7], 11+func.tp[2][7], 22+func.tp[2][6], 25+func.tp[2][6], 28+func.tp[2][6], 31+func.tp[2][6], 34+func.tp[2][6]};
	v4={17+func.tp[2][1], 8+func.tp[2][7], 12+func.tp[2][7], 23+func.tp[2][6], 26+func.tp[2][6], 29+func.tp[2][6], 32+func.tp[2][6], 35+func.tp[2][6]};
	Surface(v1, &func);
	Surface(v2, &func);
	Surface(v3, &func);
	Surface(v4, &func);
	v1={4+func.tp[1][1], 4+func.tp[1][7], 11+func.tp[1][7], 21+func.tp[1][7], 22+func.tp[1][7], 2+func.tp[1][6], 3+func.tp[1][6], 4+func.tp[1][6]};
	v2={6+func.tp[1][1], 6+func.tp[1][7], 13+func.tp[1][7], 23+func.tp[1][7], 24+func.tp[1][7], 7+func.tp[1][6], 8+func.tp[1][6], 9+func.tp[1][6], 10+func.tp[1][6], 11+func.tp[1][6]};
	v3={5+func.tp[1][1], 5+func.tp[1][7], 12+func.tp[1][7], 22+func.tp[1][7], 23+func.tp[1][7], 15+func.tp[1][6], 16+func.tp[1][6], 17+func.tp[1][6], 18+func.tp[1][6], 19+func.tp[1][6], 20+func.tp[1][6], 21+func.tp[1][6]};
	Volume(v1, &func);
	Volume(v2, &func);
	Volume(v3, &func);
	Tppush(9, 24, 16, 0, &func);

	//9
	func.tp_num = 8;
	Line(0.0, Wpml+Hsio2+Hsi, 0.0, 0.0, Wpml+Hsio2+Hsi, Wpml, &func);
	Copy(Lin, Wpml, 0.0, 0.0, "1", &func);
	Copy(Lin, Wm+Wr+g, 0.0, 0.0, "2", &func);
	Copy(Lin, Wh, 0.0, 0.0, "5", &func);
	Copy(Lin, Wm+Wcur, 0.0, 0.0, "8", &func);
	Copy(Lin, Wpml, 0.0, 0.0, "11", &func);
	
	Line(0.0, Wpml+Hsio2+Hsi, L-Wpml, 0.0, Wpml+Hsio2+Hsi, L, &func);
	Copy(Lin, Wpml, 0.0, 0.0, "17", &func);
	Copy(Lin, Wm, 0.0, 0.0, "18", &func);
	Copy(Lin, Wr, 0.0, 0.0, "21", &func);
	Copy(Lin, g+Wcur, 0.0, 0.0, "24", &func);
	Copy(Lin, Wh, 0.0, 0.0, "27", &func);
	Copy(Lin, Wm, 0.0, 0.0, "30", &func);
	Copy(Lin, Wpml, 0.0, 0.0, "33", &func);
	
	Copy(Lin, 0.0, 0.0, L-2*Wpml, "4 16", &func);
	
	Copy(Lin, 0.0, 0.0, -(L-2*Wpml-Wm), "25", &func);

	Copy(Poi, 0.0, 0.0, Wm, "6 8", &func);
	Copy(Poi, (Wh-Wpcm)/2, 0.0, 0.0, "31", &func);
	Copy(Poi, Wpcm, 0.0, 0.0, "33", &func);
	Copy(Poi, (Wh-Wpcm)/2, 0.0, 0.0, "34", &func);
	Surface("10 46 48 49 50 47", &func);

	int Leff1;
	if(Leff == 0 || Leff == Lc) Leff1 = Lc/2;
	else Leff1 = Leff;

	Copy(Lin, 0.0, 0.0, Leff1, "48 49 50", &func);
	Copy(Lin, 0.0, 0.0, Lc-Leff1, "51 52 53", &func);

	Rotate(Lin, Y, Wpml+Wm+Wr+g+Wh/2+Wcur-curv.R, Wpml+Hsio2, L-Wpml, curv.angle, "31", &func);
	Rotate(Poi, Y, Wpml+Wm+Wr+g+Wh/2+curv.R, Wpml+Hsio2, Wpml+Wm+Lc, curv.angle, "39 40", &func);
	
	Surface("58 59 60 69 65 68", &func);
	Surface("7 40 22 44 43 45 28 66 68 61 54 46", &func);
	Surface("13 41 34 67 69 62 55 47", &func);

	Tppush(0, 26, 69, 44, &func);

	//10
	func.tp_num = 6;
	Copy(Sur, 0.0, Hsi, 0.0, "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 18", &func);
	Copy(Lin, 0.0, Hsi, 0.0, "46 47 55 56", &func);
	v1={48+func.tp[2][6], 33+func.tp[2][9], 34+func.tp[2][9], 48+func.tp[2][8], 49+func.tp[2][8], 50+func.tp[2][8]};
	v2={49+func.tp[2][6], 35+func.tp[2][9], 36+func.tp[2][9], 58+func.tp[2][8], 59+func.tp[2][8], 60+func.tp[2][8]};
	v3={50+func.tp[2][6], 33+func.tp[2][9], 35+func.tp[2][9], 54+func.tp[2][8], 61+func.tp[2][8]};
	v4={51+func.tp[2][6], 34+func.tp[2][9], 36+func.tp[2][9], 55+func.tp[2][8], 62+func.tp[2][8]};
	Surface(v1, &func);
	Surface(v2, &func);
	Surface(v3, &func);
	Surface(v4, &func);
	v1={16+func.tp[1][6], 10+func.tp[1][9], 49+func.tp[1][9], 50+func.tp[1][9], 53+func.tp[1][9], 16+func.tp[1][8]};
	v2={17+func.tp[1][6], 53+func.tp[1][9], 54+func.tp[1][9], 55+func.tp[1][9], 56+func.tp[1][9], 17+func.tp[1][8], 18+func.tp[1][8], 19+func.tp[1][8], 20+func.tp[1][8], 21+func.tp[1][8], 22+func.tp[1][8]};
	v3={19+func.tp[1][6], 46+func.tp[1][9], 51+func.tp[1][9], 52+func.tp[1][9], 54+func.tp[1][9], 24+func.tp[1][8]};
	v4={20+func.tp[1][6], 7+func.tp[1][9], 40+func.tp[1][9], 22+func.tp[1][9], 44+func.tp[1][9], 43+func.tp[1][9], 45+func.tp[1][9], 28+func.tp[1][9], 47+func.tp[1][9], 51+func.tp[1][9], 55+func.tp[1][9], 49+func.tp[1][9], 25+func.tp[1][8]};
	v5={21+func.tp[1][6], 13+func.tp[1][9], 41+func.tp[1][9], 34+func.tp[1][9], 48+func.tp[1][9], 52+func.tp[1][9], 56+func.tp[1][9], 50+func.tp[1][9], 26+func.tp[1][8]};
	Volume(v1, &func);
	Volume(v2, &func);
	Volume(v3, &func);
	Volume(v4, &func);		
	Volume(v5, &func);
	Tppush(21, 56, 36, 0, &func);

	//11
	func.tp_num = 10;
	Line(Wpml+Wm+Wr+g+(Wh-Wpcm)/2, Wpml+Hsio2+Hsi+Hpcm, Wpml+Wm, Wpml+Wm+Wr+g+(Wh+Wpcm)/2, Wpml+Hsio2+Hsi+Hpcm, Wpml+Wm, &func);
	Copy(Lin, 0.0, 0.0, Leff1, "1", &func);
	Copy(Lin, 0.0, 0.0, Lc-Leff1, "2", &func);
	Tppush(0, 2, 7, 6, &func);

	//12
	func.tp_num = 8;
	Copy(Sur, 0.0, Hpcm, 0.0, "18 21", &func);
	Tppush(2, 7, 6, 0, &func);

	//13
	func.tp_num = 3;
	Copy(Sur, 0.0, -(Wm+Hpcm), 0.0, "1 2 3 7 8 9", &func);
	v1={14+func.tp[2][3], 5+func.tp[2][12], 9+func.tp[2][12], 6+func.tp[2][8], 9+func.tp[2][8], 12+func.tp[2][8]};
	v2={15+func.tp[2][3], 6+func.tp[2][12], 10+func.tp[2][12], 7+func.tp[2][8], 10+func.tp[2][8], 13+func.tp[2][8]};
	v3={16+func.tp[2][3], 7+func.tp[2][12], 11+func.tp[2][12], 22+func.tp[2][8], 25+func.tp[2][8], 28+func.tp[2][8], 31+func.tp[2][8], 34+func.tp[2][8]};
	v4={17+func.tp[2][3], 8+func.tp[2][12], 12+func.tp[2][12], 23+func.tp[2][8], 26+func.tp[2][8], 29+func.tp[2][8], 32+func.tp[2][8], 35+func.tp[2][8]};
	Surface(v1, &func);
	Surface(v2, &func);
	Surface(v3, &func);
	Surface(v4, &func);
	v1={4+func.tp[1][3], 4+func.tp[1][12], 11+func.tp[1][12], 21+func.tp[1][12], 22+func.tp[1][12], 2+func.tp[1][8], 3+func.tp[1][8], 4+func.tp[1][8]};
	v2={6+func.tp[1][3], 6+func.tp[1][12], 13+func.tp[1][12], 23+func.tp[1][12], 24+func.tp[1][12], 7+func.tp[1][8], 8+func.tp[1][8], 9+func.tp[1][8], 10+func.tp[1][8], 11+func.tp[1][8]};
	v4={5+func.tp[1][3], 5+func.tp[1][12], 12+func.tp[1][12], 22+func.tp[1][12], 23+func.tp[1][12], 15+func.tp[1][8], 16+func.tp[1][8], 17+func.tp[1][8], 19+func.tp[1][8], 20+func.tp[1][8], 22+func.tp[1][8], 23+func.tp[1][8], 24+func.tp[1][8], 25+func.tp[1][8], 26+func.tp[1][8], 1+func.tp[1][10], 2+func.tp[1][10], 1+func.tp[1][11], 3+func.tp[1][11], 4+func.tp[1][11], 5+func.tp[1][11], 6+func.tp[1][11], 7+func.tp[1][11]};
	Volume(v1, &func);
	Volume(v2, &func);
	Volume(v4, &func);
	Tppush(9, 24, 16, 0, &func);	
	
	Mat3D(1, "1 4 7 19 22 25", &func);
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
	
	Trans_Gene(trans, gene, &func);
	Printvv(func.tp);

	Filename(lpsname, &func);

	ofs.close();

	return 0;
}