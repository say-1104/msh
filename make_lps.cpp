#include "make_lps.hpp"

int main(int argc, char *argv[]){
	double trans = 0.4, gene = 0.4;
	string lpsname = "pcmDC";
	Function func;
	func.flag = 0; 

	ofs.open("pcmDC.lps");
	if(! ofs) {
		cerr << "File open error !" << endl;
		exit(1);
	}
	ofs << fixed << setprecision(4);
	cout << "writing " << lpsname << "..." << endl;

	Point(0.0, 0.0, 0.0, &func);

	vector<int> s = {1, 2, 3, 4, 5};
	Surface(s, &func); 

	Trans_Gene(trans, gene, &func);

	Filename(lpsname, &func);

	ofs.close();

	return 0;
}