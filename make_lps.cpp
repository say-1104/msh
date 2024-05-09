#include <make_lps.hpp>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <omp.h>
#include <vector>
int main(int argc, char *argv[]){
	FILE *fp;
	int i, j;

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
	double R = 15.0;

	double W = 2 * Wm + Wh + Wr + g;
	double H = Wm + Hpcm + Hsi + Hsio2;



	fp = fopen("pcmDC.lps", "w");

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
	fprintf(fp, "0 1 0\n0.000\t0.000\t0.000\t  %2.4lf\t0.000\t0.000\n21\n", Wr);
	fprintf(fp, "0 1 0\n0.000\t0.000\t0.000\t  %2.4lf\t0.000\t0.000\n22\n", g);
	fprintf(fp, "0 1 0\n0.000\t0.000\t0.000\t  %2.4lf\t0.000\t0.000\n23\n", Wh);
	fprintf(fp, "0 1 0\n0.000\t0.000\t0.000\t  %2.4lf\t0.000\t0.000\n24\n", Wm);

	fprintf(fp, "1 2 0\n0.000\t0.000\t0.000\t  0.000\t%2.4lf\t0.000\n32\n", Hsi);
	fprintf(fp, "0 1 0\n0.000\t0.000\t0.000\t  0.000\t%2.4lf\t0.000\n23 24\n", Hsi);
	fprintf(fp, "0 1 0\n0.000\t0.000\t0.000\t  %2.4lf\t0.000\t0.000\n27\n", (Wh-Wpcm)/2);
	fprintf(fp, "0 1 0\n0.000\t0.000\t0.000\t  %2.4lf\t0.000\t0.000\n29\n", Wpcm);
	fprintf(fp, "0 1 0\n0.000\t0.000\t0.000\t  %2.4lf\t0.000\t0.000\n30\n", (Wh-Wpcm)/2);
	fprintf(fp, "1 2 0\n0.000\t0.000\t0.000\t  0.000\t%2.4lf\t0.000\n42\n", Hpcm);
	fprintf(fp, "Copy END\n\n");

	fprintf(fp,"Surface\n");
	fprintf(fp,"34 39 41 42 43 40\n");
	fprintf(fp,"12 28 13 35 34 33 32 31\n");
	fprintf(fp,"18 31 37 36 38 33 39 41 45 44 46 43 40 35 19 29\n");
	fprintf(fp,"Surface END\n\n");

	fprintf(fp,"Copy\n");
		fprintf(fp, "2 3 0\n0.000\t0.000\t0.000\t  0.000\t0.000\t%2.4lf\n", Wpml);
		for (i=0; i<15; i++) fprintf(fp,"%d ", i+1);
		fprintf(fp,"\n");
	if(Leff == Lc){
		fprintf(fp, "2 3 0\n0.000\t0.000\t0.000\t  0.000\t0.000\t%2.4lf\n", Leff);
		for (i=0; i<15; i++) fprintf(fp,"%d ", i+62);
		fprintf(fp,"\n");
		fprintf(fp, "2 3 0\n0.000\t0.000\t0.000\t  0.000\t0.000\t%2.4lf\n", Wpml);
		for (i=0; i<15; i++) fprintf(fp,"%d ", i+123);
		fprintf(fp,"\n");
	}else if(Leff > 0){
		fprintf(fp, "2 3 0\n0.000\t0.000\t0.000\t  0.000\t0.000\t%2.4lf\n", Leff);
		for (i=0; i<15; i++) fprintf(fp,"%d ", i+62);
		fprintf(fp,"\n");
		fprintf(fp, "2 3 0\n0.000\t0.000\t0.000\t  0.000\t0.000\t%2.4lf\n", Lc-Leff);
		for (i=0; i<15; i++) fprintf(fp,"%d ", i+123);
		fprintf(fp,"\n");
		fprintf(fp, "2 3 0\n0.000\t0.000\t0.000\t  0.000\t0.000\t%2.4lf\n", Wpml);
		for (i=0; i<15; i++) fprintf(fp,"%d ", i+184);
		fprintf(fp,"\n");
	} else {
		fprintf(fp, "2 3 0\n0.000\t0.000\t0.000\t  0.000\t0.000\t%2.4lf\n", Lc);
		for (i=0; i<15; i++) fprintf(fp,"%d ", i+62);
		fprintf(fp,"\n");
		fprintf(fp, "2 3 0\n0.000\t0.000\t0.000\t  0.000\t0.000\t%2.4lf\n", Wpml);
		for (i=0; i<15; i++) fprintf(fp,"%d ", i+123);
		fprintf(fp,"\n");
	}
	fprintf(fp, "Copy END\n\n");


	fprintf(fp,"Material\n");
	fprintf(fp,"3 1\n1 2 3 4 9 14\n");
	fprintf(fp,"3 2\n5 6 7 8 10 12 15\n");
	fprintf(fp,"3 3\n13\n");
	fprintf(fp,"3 4\n11\n");

	if(Lc > Leff && Leff > 0){
		fprintf(fp,"3 5\n46 47 48 49 54 59\n");
		fprintf(fp,"3 6\n50 51 52 53 55 57 60\n");
		fprintf(fp,"3 7\n58\n");
		fprintf(fp,"3 8\n56\n");
		fprintf(fp,"3 9\n28 43\n");
		fprintf(fp,"3 10\n26 41\n");
		fprintf(fp,"3 11\n16 17 18 19 24 29 31 32 33 34 39 44\n");
		fprintf(fp,"3 12\n20 21 22 23 25 30 35 36 37 38 40 45\n");
		fprintf(fp,"3 13\n27\n");
		fprintf(fp,"3 14\n42\n");
	} else {
		fprintf(fp,"3 5\n31 32 33 34 39 44\n");
		fprintf(fp,"3 6\n35 36 37 38 40 42 45\n");
		fprintf(fp,"3 7\n43\n");
		fprintf(fp,"3 8\n41\n");
		fprintf(fp,"3 9\n28\n");
		fprintf(fp,"3 10\n26\n");
		fprintf(fp,"3 11\n16 17 18 19 24 29\n");
		fprintf(fp,"3 12\n20 21 22 23 25 30\n");
		if(Leff == Lc) fprintf(fp,"3 13\n27\n");
		else fprintf(fp,"3 14\n27\n");	
	}
	fprintf(fp,"Material END\n\n");

	fprintf(fp,"Unstr_mesh\n");
	fprintf(fp,"2 0.04\n");
	i=15;
	fprintf(fp,"11 12 13 %d %d %d %d %d %d %d %d %d %d %d %d %d ", 32+i, 34+i, 36+i, 37+i, 38+i, 39+i, 40+i, 41+i, 42+i, 43+i, 44+i, 45+i, 46+i);
	i=76;
	fprintf(fp,"72 73 74 %d %d %d %d %d %d %d %d %d %d %d %d %d ", 32+i, 34+i, 36+i, 37+i, 38+i, 39+i, 40+i, 41+i, 42+i, 43+i, 44+i, 45+i, 46+i);
	i=137;
	fprintf(fp,"133 134 135 %d %d %d %d %d %d %d %d %d %d %d %d %d ", 32+i, 34+i, 36+i, 37+i, 38+i, 39+i, 40+i, 41+i, 42+i, 43+i, 44+i, 45+i, 46+i);
	fprintf(fp,"194 195 196");
	if(Lc > Leff && Leff > 0){
		i=198;
		fprintf(fp,"194 195 196 %d %d %d %d %d %d %d %d %d %d %d %d %d ", 32+i, 34+i, 36+i, 37+i, 38+i, 39+i, 40+i, 41+i, 42+i, 43+i, 44+i, 45+i, 46+i);
		fprintf(fp,"255 256 257");	
	}
	fprintf(fp,"\n");	

	
	fprintf(fp,"Unstr_mesh END\n\n");
	
	fprintf(fp,"TransitionFactor\n");
	fprintf(fp,"0.4 0\n");

	fprintf(fp,"Generate\n");
	fprintf(fp,"0.4\n");

	fprintf(fp,"saveMesh\n\"pcmDC\"\nexportMesh\n\"pcmDC\"\n");
}