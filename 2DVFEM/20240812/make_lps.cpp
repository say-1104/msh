#include <LpsTools.hpp>

using namespace std;

int main(int argc, char *argv[]){
    //導波路のパラメータ
	double w = atof(argv[1]);		//導波路幅
	double h = 0.22;				//導波路高さ
	double hsio2 = 2.5;

	//その他パラメータ
	double wm = 2.5;		//導波路とメッシュ端の距離

	double Xmax = wm*2 + w, Ymax = wm + h + hsio2;
    double unstr = 0.01, trans = 0.4, gene = 0.1;
    vector<int> v1, v2, v3, v4, v5, v6;

    LpsTools* lt = new LpsTools("wire", unstr, trans, gene);

	lt->Line(0.0, 0.0, 0.0,				Xmax, 0.0, 0.0);
	lt->Line(Xmax, 0.0, 0.0,			Xmax, Ymax, 0.0);
	lt->Line(Xmax, Ymax, 0.0,			0.0, Ymax, 0.0);
	lt->Line(0.0, Ymax, 0.0,			0.0, 0.0, 0.0);

	lt->Line(wm, hsio2, 0.0,			wm+w, hsio2, 0.0);
	lt->Line(wm+w, hsio2, 0.0,			wm+w, hsio2+h, 0.0);
	lt->Line(wm+w, hsio2+h, 0.0,		wm, hsio2+h, 0.0);
	lt->Line(wm, hsio2+h, 0.0,			wm, hsio2, 0.0);

	lt->Surface("5 6 7 8");
	lt->Surface("1 2 3 4 5 6 7 8");

	lt->Mat2D(1, "1");
	lt->Mat2D(2, "2");

	lt->Unstr(Shape::surface, "1");

    lt->TransGene();
    lt->Fileclose();

    return 0;
}
