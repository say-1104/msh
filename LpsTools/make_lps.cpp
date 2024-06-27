#include <LpsTools.hpp>

using namespace std;

int main(int argc, char *argv[]){
    //DCのパラメータ(外部入力)
	double Wh = atof(argv[1]);		//導波路幅(PCM有)
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
	double L = 2 * Wpml + 2 * Wz + Lc, W = 2 * Wpml + 2 * Wx + g + 2 * Wh, H = 2 * Wpml + Hsio2 + Hsi + Wy;
    double unstr = 0.04, trans = 0.3, gene = 0.3;
    vector<int> v1, v2, v3, v4, v5, v6;

    LpsTools* lt = new LpsTools("DC", unstr, trans, gene);

    auto lps1 = [&](int n, double Y) -> void {
		lt->step_offset = n;
        lt->y_offset = Y;
		lt->Line(0.0, 0.0, 				0.0, Wpml);
		lt->Line(0.0, Wpml, 			0.0, L-Wpml);
		lt->Line(0.0, L-Wpml, 		    0.0, L);
		lt->Copy(Shape::line, Wpml, 0.0, 0.0, 		"1 2 3");
		lt->Copy(Shape::line, W-2*Wpml, 0.0, 0.0, 	"4 5 6");
		lt->Copy(Shape::line, Wpml, 0.0, 0.0, 		"11 12 13");
		lt->Appendstep(0, 9, 24, 16);
	};

	auto lps2 = [&](int n, double Y) -> void {
		lt->step_offset = n;
		lt->Copy(Shape::surface, 0.0, Y, 0.0, 			"1 2 3 4 5 6 7 8 9");
		lt->Appendstep(9, 24, 16, 0);
	};

	//step 1, 2, 3
    lps1(1, 0.0); lps1(2, Wpml); lps2(1, Wpml);

	//step 4, 5, 6
	lps1(4, H-Wpml); lps1(5, H); lps2(4, Wpml);

    lt->TransGene();
    lt->Fileclose();
	//lab
    return 0;
}
