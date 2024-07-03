#include <LpsTools.hpp>

using namespace std;

int main(int argc, char *argv[]){
    //DCのパラメータ(外部入力)
	double W = atof(argv[1]);		//導波路幅(PCM有)
	double g = atof(argv[2]);		//導波路間隔
	double Lc = atof(argv[3]);		//導波路長

	//DCのパラメータ(固定値)
	double Hsi = 0.22;				//Si層の厚さ
	double Hsio2 = 2.0;				//SiO2層の厚さ
	double Wx = 1.5;				//導波路とPML間の距離
	double Wy = 1.5;				//導波路とPML間の距離
	double Wz = 0.5;				//導波路とPML間の距離
	double Wpml = 1.0;				//PML幅

    int div = 10;
	double dLc = Lc / div;
	double Zmax = 2 * Wpml + 2 * Wz + Lc, Xmax = 2 * Wpml + 2 * Wx + g + 2 * W, Ymax = 2 * Wpml + Hsio2 + Hsi + Wy;
    double unstr = 0.04, trans = 0.3, gene = 0.3;
    vector<int> v1, v2, v3, v4, v5, v6;

    LpsTools* lt = new LpsTools("DC", unstr, trans, gene);

    auto lps1 = [&](int n, double Y) -> void {
		lt->step_offset = n;
        lt->y_offset = Y;
		lt->Line(0.0, 0.0, 				0.0, Wpml);
		lt->Line(0.0, Wpml, 			0.0, Zmax-Wpml);
		lt->Line(0.0, Zmax-Wpml, 		0.0, Zmax);
		lt->Copy(Shape::line, Wpml, 0.0, 0.0, 			"1 2 3");
		lt->Copy(Shape::line, Xmax-2*Wpml, 0.0, 0.0, 	"4 5 6");
		lt->Copy(Shape::line, Wpml, 0.0, 0.0, 			"11 12 13");
		lt->Appendstep(16, 24, 9, 0);
	};

	auto lps2 = [&](int n, double Y) -> void {
		lt->step_offset = n;
		lt->Copy(Shape::surface, 0.0, Y, 0.0, 			"1 2 3 4 5 6 7 8 9");
		lt->Appendstep(0, 16, 24, 9);
	};

	auto lps3 = [&](int n, double Y) -> void {
		lt->step_offset = n;
		lt->y_offset = Y;
		lt->Line(0.0, 0.0, 				0.0, Wpml);
		lt->Copy(Shape::line, Wpml, 0.0, 0.0, 		"1");
		lt->Copy(Shape::line, Wx, 0.0, 0.0, 		"2");
		lt->Copy(Shape::line, W, 0.0, 0.0, 			"5");
		lt->Copy(Shape::line, g, 0.0, 0.0, 			"8");
		lt->Copy(Shape::line, W, 0.0, 0.0, 			"11");
		lt->Copy(Shape::line, Wx, 0.0, 0.0, 		"14");
		lt->Copy(Shape::line, Wpml, 0.0, 0.0, 		"17");

		lt->Line(0.0, Zmax-Wpml, 		0.0, Zmax);
		lt->Copy(Shape::line, Wpml, 0.0, 0.0, 		"23");
		lt->Copy(Shape::line, Wx, 0.0, 0.0, 		"24");
		lt->Copy(Shape::line, W, 0.0, 0.0, 			"27");
		lt->Copy(Shape::line, g, 0.0, 0.0, 			"30");
		lt->Copy(Shape::line, W, 0.0, 0.0, 			"33");
		lt->Copy(Shape::line, Wx, 0.0, 0.0, 		"36");
		lt->Copy(Shape::line, Wpml, 0.0, 0.0, 		"39");
		lt->Copy(Shape::line, 0.0, 0.0, Zmax-2*Wpml, 	"4 22");

		lt->Copy(Shape::line, 0.0, 0.0, Wz, 			"10 13 16");
		for(int i=0; i<div; i++){
			v1 = {49+i*7, 50+i*7, 51+i*7};
			lt->Copy(Shape::line, 0.0, 0.0, dLc, 	Plusv(v1, lt->step[1][lt->step_offset]));
		}
		v1 = {49+7*div, 50+7*div, 51+7*div};
		lt->Copy(Shape::line, 0.0, 0.0, Wz, 	Plusv(v1, lt->step[1][lt->step_offset]));

		v1={28, 46, 7};
		v2={40, 47, 19};
		for(int i=0; i<div+1; i++){
			v1.push_back(52+i*7);
			v2.push_back(55+i*7);
		}
		v1.push_back(56+div*7);
		v2.push_back(59+div*7);
		lt->Surface(Plusv(v1, lt->step[1][lt->step_offset]));
		lt->Surface(Plusv(v2, lt->step[1][lt->step_offset]));
		lt->Appendstep(36+4*div, 59+div*7, 24+3*div, 0);

	};
	
	auto lps4 = [&](int n, double Y) -> void {
		lt->step_offset = n;
		lt->Copy(Shape::surface, 0.0, Y, 0.0, 			1, 3*div+24);
		lt->Appendstep(0, 36+4*div, 59+div*7, 24+3*div);
	};

	auto lps5 = [&](int n1, int n2, int n3, double Y) -> void {
		lt->step_offset = n1;
		lt->Copy(Shape::surface, 0.0, Y, 0.0, 			"1 2 3 7 8 9");
		v1={14+lt->step[1][n1], 5+lt->step[1][n2], 9+lt->step[1][n2], 6+lt->step[1][n3], 9+lt->step[1][n3], 12+lt->step[1][n3], 15+lt->step[1][n3], 18+lt->step[1][n3]};
		v2={15+lt->step[1][n1], 6+lt->step[1][n2], 10+lt->step[1][n2], 7+lt->step[1][n3], 10+lt->step[1][n3], 13+lt->step[1][n3], 16+lt->step[1][n3], 19+lt->step[1][n3]};
		v3={16+lt->step[1][n1], 7+lt->step[1][n2], 11+lt->step[1][n2], 28+lt->step[1][n3], 31+lt->step[1][n3], 34+lt->step[1][n3], 37+lt->step[1][n3], 40+lt->step[1][n3]};
		v4={17+lt->step[1][n1], 8+lt->step[1][n2], 12+lt->step[1][n2], 29+lt->step[1][n3], 32+lt->step[1][n3], 35+lt->step[1][n3], 38+lt->step[1][n3], 41+lt->step[1][n3]};
		lt->Surface(v1);
		lt->Surface(v2);
		lt->Surface(v3);
		lt->Surface(v4);

		v1={4+lt->step[1][n1], 4+lt->step[1][n2], 11+lt->step[1][n2], 21+lt->step[1][n2], 22+lt->step[1][n2], 2+lt->step[1][n3], 3+lt->step[1][n3], 4+lt->step[1][n3], 5+lt->step[1][n3], 6+lt->step[1][n3]};
		v2={5+lt->step[1][n1], 5+lt->step[1][n2], 12+lt->step[1][n2], 22+lt->step[1][n2], 23+lt->step[1][n2]};
		v3={6+lt->step[1][n1], 6+lt->step[1][n2], 13+lt->step[1][n2], 23+lt->step[1][n2], 24+lt->step[1][n2], 9+lt->step[1][n3], 10+lt->step[1][n3], 11+lt->step[1][n3], 12+lt->step[1][n3], 13+lt->step[1][n3]};
		for(int i=0; i<3*div+8; i++){
			v2.push_back(17+i);
		}
		lt->Volume(v1);
		lt->Volume(v2);
		lt->Volume(v3);
		lt->Appendstep(0, 16, 24, 9);
	};

	//step 1, 2, 3
    lps1(1, 0.0); lps1(2, Wpml); lps2(1, Wpml);

	//step 4, 5, 6
	lps1(4, Ymax-Wpml); lps1(5, Ymax); lps2(4, Wpml);

	//step 7, 8, 9
	lps3(7, Wpml+Hsio2); lps3(8, Wpml+Hsio2+Hsi); lps4(7, Hsi);

	//step 10, 11
	lps5(2, 10, 7, Hsio2); lps5(4, 11, 8, -Wy);

    //lt->TransGene();
    lt->Fileclose();
	//lab
    return 0;
}
