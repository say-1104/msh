#include "make_lps.hpp"

void makePML(int argc, char *argv[], double ypml, int in_num, int out_num, int tmpl){

}


int main(int argc, char *argv[]){
	Function func;
	func.flag = 0; 

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
	double R = 14.5;

	//メッシュのパラメータ
	double trans = 0.4, gene = 0.4;
	string lpsname = "pcmDC";

	cout << "writing " << lpsname << "..." << endl;


	Point(0.0, 0.0, 0.0, &func);

	Surface("1 2 3 4 5", &func); 

	Trans_Gene(trans, gene, &func);

	Filename(lpsname, &func);

	ofs.close();

	return 0;
}