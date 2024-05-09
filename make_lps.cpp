#include "make_lps.hpp"

int main(int argc, char *argv[]){
	string lpsname = "pcmDC";
	string filename = "pcmDC.lps";

	ofstream ofs;
	ofs.open( "pcmDC.lps");

	if(! ofs) {
		cerr << "File open error !" << endl;
		exit(1);
	}



	cout << "writing " << filename << "..." << endl;

	for (int i = 0; i<10; i++){
		ofs << i << endl;
	}
	ofs.close();
	return 0;
}