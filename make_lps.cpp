#include "make_lps.hpp"

void Printvec(vector<vector<int> > a){
	for(int i=0; i<a.size(); i++){
		for(int j=0; j<a[i].size(); j++){
			cout << a[i][j] << " " ;
		}
		cout << endl;
	}
	cout << endl;
}

int main(int argc, char *argv[]){
	Function func;
	func.flag = 0;
	func.tp_num = 0;
	func.Tppush(0, 9, 24, 16);
	Printvec(func.tp);
	func.Tppush(0, 9, 24, 16);
	ofs.open("pcmDC.lps");
	if(! ofs) {
		cerr << "File open error !" << endl;
		exit(1);
	}
	ofs << fixed << setprecision(4);
	
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
	double R = 14.5;				//曲げ半径

	double L = 2 * Wpml + 2 * Wm + Lc + R + Wh / 2;
	double W = 2 * Wpml + Wm + Wr + g + Wh / 2 + R;
	double H = 2 * Wpml + Hsio2 + Hsi + Hpcm + Wm;

	//メッシュのパラメータ
	double trans = 0.4, gene = 0.4;
	string lpsname = "pcmDC";

	cout << "writing " << lpsname << "..." << endl;

	//1

	Line(0.0, 0.0, 0.0, 0.0, 0.0, Wpml, &func);
	Line(0.0, 0.0, Wpml, 0.0, 0.0, L-Wpml, &func);
	Line(0.0, 0.0, L-Wpml, 0.0, 0.0, L, &func);
	Copy(Lin, Wpml, 0.0, 0.0, "1 2 3", &func);
	Copy(Lin, W-2*Wpml, 0.0, 0.0, "4 5 6", &func);
	Copy(Lin, Wpml, 0.0, 0.0, "11 12 13", &func);
	
	//2, 3
	Copy(Sur, 0.0, Wpml, 0.0, "1 2 3 4 5 6 7 8 9", &func);

	Printvec(func.tp);
	func.Tppush(9, 24, 16, 0);
	Printvec(func.tp);
	//4
	Line(0.0, H-Wpml, 0.0, 0.0, H-Wpml, Wpml, &func);
	Line(0.0, H-Wpml, Wpml, 0.0, H-Wpml, L-Wpml, &func);
	Line(0.0, H-Wpml, L-Wpml, 0.0, H-Wpml, L, &func);
	func.tp_num = 3;
	Copy(Lin, Wpml, 0.0, 0.0, "1 2 3", &func);
	Copy(Lin, W-2*Wpml, 0.0, 0.0, "4 5 6", &func);
	Copy(Lin, Wpml, 0.0, 0.0, "11 12 13", &func);
	func.Tppush(0, 9, 24, 16);
	Printvec(func.tp);
	//5, 16
	Copy(Sur, 0.0, Wpml, 0.0, "1 2 3 4 5 6 7 8 9", &func);
	func.Tppush(0, 9, 24, 16);
	Printvec(func.tp);
	func.Tppush(9, 24, 16, 0);

	Printvec(func.tp);

	//5
	/*Line(0.0, Wpml+Hsio2, 0.0, 0.0, Wpml+Hsio2, Wpml, &func);
	Copy(Lin, Wpml, 0.0, 0.0, "1", &func);
	Copy(Lin, Wm+Wr+g, 0.0, 0.0, "2", &func);
	Copy(Lin, Wh, 0.0, 0.0, "5", &func);
	Copy(Lin, R-Wh/2, 0.0, 0.0, "8", &func);
	Copy(Lin, Wpml, 0.0, 0.0, "11", &func);
	
	Line(0.0, Wpml+Hsio2, L-Wpml, 0.0, Wpml+Hsio2, L, &func);
	Copy(Lin, Wpml, 0.0, 0.0, "17", &func);
	Copy(Lin, Wm, 0.0, 0.0, "18", &func);
	Copy(Lin, Wr, 0.0, 0.0, "21", &func);
	Copy(Lin, g+R+Wh/2, 0.0, 0.0, "24", &func);
	Copy(Lin, Wpml, 0.0, 0.0, "27", &func);
	
	Copy(Lin, 0.0, 0.0, L-2*Wpml, "4", &func);
	Copy(Lin, 0.0, 0.0, L-2*Wpml-Wm-Wh, "16", &func);
	Copy(Lin, 0.0, 0.0, Wh, "35", &func);
	Copy(Lin, 0.0, 0.0, Wm, "38", &func);
	
	Copy(Lin, 0.0, 0.0, -(L-2*Wpml-Wm), "25", &func);

	Copy(Poi, 0.0, 0.0, Wm, "6 8", &func);
	Copy(Poi, Wh, 0.0, 0.0, "31", &func);
	Surface("10 46 48 47", &func);
	Copy(Lin, 0.0, 0.0, Lc, "48", &func);
	Rotate(Lin, Y, W-Wpml, Wpml+Hsio2, Wpml+Wm+Lc, 90, "49", &func);

	Surface("7 34 22 44 43 45 28 41 52 50 46", &func);
	Surface("13 36 53 51 47", &func);*/




	//Trans_Gene(trans, gene, &func);

	Filename(lpsname, &func);

	ofs.close();

	return 0;
}