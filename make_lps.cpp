#include "make_lps.hpp"

int main(int argc, char *argv[]){
	double trans = 0.4, gene = 0.4;
	string lpsname = "pcmDC";
	Function func; 

	ofs.open("pcmDC.lps");
	if(! ofs) {
		cerr << "File open error !" << endl;
		exit(1);
	}
	cout << "writing " << lpsname << "..." << endl;

	Point(0.0, 0.0, 0.0, &func);

	Trans_Gene(trans, gene, &func);

	Filename(lpsname, &func);

	ofs.close();

	return 0;
}