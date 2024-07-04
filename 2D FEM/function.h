/*   FEManalysis.c   */
extern int elementMatrix(DataTable *data, double *x, double *y,
			double_complex ek[6][6], int matID, double wl);
extern int Amatrix(DataTable *data, CGtable *cg_table, double wl);
extern int Bmatrix(PortData *port, DataTable *data, double wl);
extern int Fvector(DataTable *data, double wl, int ll, int ii);
extern int FEManalysis(DataTable *data, double wl, int ll);
extern int rotation(Tensor *ep, Tensor *epR, double cosA, double sinA);

/*   FEMbasic.c   */
extern int matNxNx(double g1[3][3], double wl, double *p);
extern int matNN(double g6[3][3], double wl, double *p);
extern int matNxN(double g6[3][3]);
extern int matNNx(double g6[3][3]);
extern int matNNc(double_complex g6[3][3], double wl, double_complex *p);
extern int matNNxc(double_complex g6[3][3]);

/*   alloc.c   */
extern int allocMatrix(DataTable *data);

/*   bicgstab.c   */
static int sorting(int *aa, int nn);
extern int makeCGmatrixTable(Element3D *element, int ne, int np,
			     CGtable *cg_table, int NK);
extern int columnInCGmatrix(int ii, int jj, CGtable *cg_table, int np);

/*   inputData.c   */
extern int checkArgument1(DataTable *data, int argc, char **argv);
extern int checkArgument2(DataTable *data, int argc, char **argv);
extern int checkBandWidth(FEM *fem, DataTable *data);
extern int inputData(DataTable *data, char *file);
extern int rotationCoordinate(double x0, double y0,
				 double *xr0, double *yr0, double angle);

/*   optfem.c   */
extern int main(int argc, char **argv);

/*   output.c   */
extern int Preforintpol(DataTable *data);
extern int intpolFieldValue(DataTable *data, double x0, double y0, 
			    double xmin, double xmax, double ymin, double ymax,
			    double_complex phi[6]);
extern int intpol(DataTable *data, double w_c, int mm);
extern int intpolXX(DataTable *data, double wl, int mm);
extern int intpolYY(DataTable *data, double wl, int mm);

/*   port.c   */
extern int getMaterialIdOfNode(int nn, DataTable *data, int ii);
extern int portInfo(PortData *port, DataTable *data, int ii);
extern int setPort(DataTable *data, double wl, int ll, int mm);
extern int eigenModeOfPort(DataTable *data, double wl, int ll, int mm);

extern int eigenModeOfPortByFileInput(DataTable *data, PortData *port,
				      double wl, int ll, int mm, int ii);
extern int Pointsource(DataTable *data);
extern int Sheetsource(DataTable *data);
extern double normalize(DataTable *data, PortData *port, int ll, int ii);

/*   probability.c   */
extern int portField(DataTable *data, int ll);
extern int probability(DataTable *data, double wl, int ll, int mm);
extern double_complex overlapIntegral(PortData *port1, PortData *port2, int m);
extern int integral_elem(double xx[3], double_complex FF1[3], double_complex FF2[3], 
		   double le, double pp, double_complex *elementpower);
extern int outputOverlap(DataTable *data, int kk);

// mosaic.c
extern int inputMosaic(DataTable *data);
extern int checkMeshforMosaic(DataTable *data);
extern int shuffle(int *array, int size);

// modifiedEIM.c
extern int inputModifiedEIMData(DataTable *data);
extern int SetModifiedEI(DataTable *data, double wl, int ll, int mm);

// random.c
void advance_random();
int flip(double prob);
void randomize();
double randomperc();
int rnd(int low, int high);
double rndreal(double lo, double hi);
void warmup_random(double random_seed);
double rand_normal( double mu, double sigma );

// ga.c 
extern int DirectBinarySearch(DataTable *data);
extern int calcFEM(DataTable *data);
// callPARDISO.c
extern int callPARDISO(DataTable *data, std::complex<double> **A,
                       std::complex<double> *bb, int nn,
                       std::complex<double> *phi, CGtable *cg_table,
                       int c_count);
