#include <LpsTools.hpp>

using namespace std;

int main(int argc, char *argv[]){
    //MMIのパラメータ
	double w = 0.5;			//入出力導波路幅
	double W = 4.5;			//MMI導波路幅
	double L = 20.0;		//MMI導波路長

	//その他パラメータ
	double Wm = 2.0;		//導波路とPML間の距離
	double Wpml = 1.0;		//PML幅
	double pix = 0.5;		//ピクセルサイズ

	double Xmax = 2 * Wpml + 2 * Wm + L, Ymax = 2 * Wpml + 2 * Wm + W;
    double unstr = 0.05, trans = 0.4, gene = 0.1;
    vector<int> v1, v2, v3, v4, v5, v6;

    LpsTools* lt = new LpsTools("mmi", unstr, trans, gene);

	lt->Line(0.0, 0.0, 0.0,				Wpml, 0.0, 0.0);
	lt->Copy(Shape::line, 0.0, Wpml, 0.0, 				"1");
	lt->Copy(Shape::line, 0.0, Wm + (W-w)/2, 0.0, 		"2");
	lt->Copy(Shape::line, 0.0, w, 0.0, 					"5");
	lt->Copy(Shape::line, 0.0, Wm + (W-w)/2, 0.0, 		"8");
	lt->Copy(Shape::line, 0.0, Wpml, 0.0, 				"11");

	lt->Line(Xmax-Wpml, 0.0, 0.0,		Xmax, 0.0, 0.0);
	lt->Copy(Shape::line, 0.0, Wpml, 0.0, 				"17");
	lt->Copy(Shape::line, 0.0, Wm, 0.0, 				"18");
	lt->Copy(Shape::line, 0.0, w, 0.0, 					"21");
	lt->Copy(Shape::line, 0.0, W - w*2, 0.0, 			"24");
	lt->Copy(Shape::line, 0.0, w, 0.0, 					"27");	
	lt->Copy(Shape::line, 0.0, Wm, 0.0, 				"30");	
	lt->Copy(Shape::line, 0.0, Wpml, 0.0, 				"33");

	lt->Copy(Shape::line, Xmax - 2*Wpml, 0.0, 0.0, 		"4 16");

	lt->Line(Wpml+Wm, Wpml+Wm, 0.0,				Wpml+Wm, Wpml+Wm+pix, 0.0);
	for(int i=0; i<8; i++){
		lt->Copy(Shape::point, 0.0, pix, 0.0, 				{30+i});
	}
	lt->Copy(Shape::line, pix, 0.0, 0.0, 					"43 44 45 46 47 48 49 50 51");
	for(int i=0; i<39; i++){
		v1 = {52, 53, 54, 55, 56, 57, 58, 59, 60};
		v1 = Plusv(v1, 19*i);
		lt->Copy(Shape::line, pix, 0.0, 0.0, 				v1);
	}
	lt->Copy(Shape::line, Wm, 0.0, 0.0, 		"10");
	lt->Copy(Shape::line, -Wm, 0.0, 0.0, 		"25 31");
	v1={28, 815, 794, 795, 796, 797, 798, 799, 800, 816};
	v2={814, 22, 40, 7, 812, 46, 45, 44, 43};
	v3={817, 34, 41, 13, 813, 48, 49, 50, 51};
	for(int i=0; i<40; i++){
		v2.push_back(61+i*19);
		v3.push_back(70+i*19);
	}
	lt->Surface(v1);
	lt->Surface(v2);
	lt->Surface(v3);

	lt->Mat2D(1, "1 2 4 5 6 7 9 11 12 13 14");
	lt->Mat2D(2, "3 8 10");
	lt->Mat2D(3, "378 379 380");
	lt->Mat2D(4, "375 376 377");
	v1={};
	for(int i=0; i<360; i++){
		v1.push_back(15+i);
	}
	lt->Mat2D(5, v1);
	v1.push_back(3); v1.push_back(8); v1.push_back(10); 
	v1.push_back(375); v1.push_back(376); v1.push_back(377); 
	lt->Unstr(Shape::surface, v1);

    lt->TransGene();
    lt->Fileclose();

    return 0;
}
