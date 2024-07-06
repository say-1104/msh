#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <omp.h>

int main(int argc, char *argv[]){
	FILE *fp;
	int i, j;

	//DCのパラメータ(外部入力)
	double Wh = atof(argv[1]);		//導波路幅(PCM有)
	double Wr = atof(argv[2]);		//導波路幅(PCM無)
	double Wpcm = atof(argv[3]);	//PCM幅
	double g = atof(argv[4]);		//導波路間隔
	double Hpcm = atof(argv[5]);	//PCM層の厚さ

	//DCのパラメータ(固定値)
	double Hsi = 0.22;				//Si層の厚さ
	double Hsio2 = 2.0;				//SiO2層の厚さ
	double Wm = 1.5;				//導波路とPML間の距離
	double Wpml = 1.0;				//PML幅

	double W = 2 * Wm + Wh + Wr + g;
	double H = Wm + Hpcm + Hsi + Hsio2;



	fp = fopen("wire.lps", "w");

	fprintf(fp,"Line\n");
	fprintf(fp,"0\t 0.000\t0.000\t0.000\t %2.4lf\t0.000\t0.000\n", Wpml);
	fprintf(fp,"0\t %2.4lf\t0.000\t0.000\t %2.4lf\t0.000\t0.000\n", Wpml+W, W+2*Wpml);
	fprintf(fp,"Line END\n\n");

	fprintf(fp,"Copy\n");
	fprintf(fp, "1 2 0\n0.000\t0.000\t0.000\t  0.000\t%2.4lf\t0.000\n1 2\n", Wpml);
	fprintf(fp, "1 2 0\n0.000\t0.000\t0.000\t  0.000\t%2.4lf\t0.000\n3 4\n", Hsio2);
	fprintf(fp, "1 2 0\n0.000\t0.000\t0.000\t  0.000\t%2.4lf\t0.000\n9 10\n", H-Hsio2);
	fprintf(fp, "1 2 0\n0.000\t0.000\t0.000\t  0.000\t%2.4lf\t0.000\n15 16\n", Wpml);
	fprintf(fp, "1 2 0\n0.000\t0.000\t0.000\t  %2.4lf\t0.000\t0.000\n6 24\n", W);

	fprintf(fp, "0 1 0\n0.000\t0.000\t0.000\t  %2.4lf\t0.000\t0.000\n10\n", Wm);
	fprintf(fp, "0 1 0\n0.000\t0.000\t0.000\t  %2.4lf\t0.000\t0.000\n21\n", Wh);
	fprintf(fp, "0 1 0\n0.000\t0.000\t0.000\t  %2.4lf\t0.000\t0.000\n22\n", g);
	fprintf(fp, "0 1 0\n0.000\t0.000\t0.000\t  %2.4lf\t0.000\t0.000\n23\n", Wr);
	fprintf(fp, "0 1 0\n0.000\t0.000\t0.000\t  %2.4lf\t0.000\t0.000\n24\n", Wm);

	fprintf(fp, "1 2 0\n0.000\t0.000\t0.000\t  0.000\t%2.4lf\t0.000\n34\n", Hsi);
	fprintf(fp, "0 1 0\n0.000\t0.000\t0.000\t  0.000\t%2.4lf\t0.000\n21 22\n", Hsi);
	fprintf(fp, "0 1 0\n0.000\t0.000\t0.000\t  %2.4lf\t0.000\t0.000\n27\n", (Wh-Wpcm)/2);
	fprintf(fp, "0 1 0\n0.000\t0.000\t0.000\t  %2.4lf\t0.000\t0.000\n29\n", Wpcm);
	fprintf(fp, "0 1 0\n0.000\t0.000\t0.000\t  %2.4lf\t0.000\t0.000\n30\n", (Wh-Wpcm)/2);
	fprintf(fp, "1 2 0\n0.000\t0.000\t0.000\t  0.000\t%2.4lf\t0.000\n42\n", Hpcm);
	fprintf(fp, "Copy END\n\n");

	fprintf(fp,"Surface\n");
	fprintf(fp,"32 39 40 41 42 43\n");
	fprintf(fp,"12 28 13 35 34 33 32 31\n");
	fprintf(fp,"18 31 39 41 45 44 46 43 40 33 37 36 38 35 19 29\n");
	fprintf(fp,"Surface END\n\n");
	
	fprintf(fp,"Material\n");
	fprintf(fp,"2 1\n1 2 3 4 9 14\n");
	fprintf(fp,"2 2\n5 6 7 8 10 15\n");
	fprintf(fp,"2 3\n11\n");
	fprintf(fp,"2 4\n13\n");
	fprintf(fp,"2 5\n12\n");
	fprintf(fp,"Material END\n\n");
	
	fprintf(fp,"Unstr_mesh\n");
	fprintf(fp,"2 0.01\n11 12 13\n");
	fprintf(fp,"Unstr_mesh END\n\n");
	
	fprintf(fp,"TransitionFactor\n");
	fprintf(fp,"0.4 0\n");
	
	fprintf(fp,"Generate\n");
	fprintf(fp,"0.1\n");
	
	fprintf(fp,"saveMesh\n\"wire\"\nexportMesh\n\"wire\"\n");
}