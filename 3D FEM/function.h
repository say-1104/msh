/*   analysis.c   */
void analysis(DataTable *data);
void memoryAllocation(DataTable *data);
long long int checkBandWidth(DataTable *data);
long long int maxValue(long long int *a, long long int n, long long int nr);
long long int minValue(long long int *a, long long int n, long long int nr);
void matrix(DataTable *data);
void reArrangement(DataTable *data);
void mprocess(DataTable *data, std::complex<double> ***array,long long int k);
void inputCondition(DataTable *data);

// analysis2.c
void analysis2(DataTable *data, int flag);
void inputCondition2(DataTable *data);
void solveMatrixEquation2(DataTable *data);
long long int wsmp_pre2(std::complex<double> **A, std::complex<double> *bb, 
			long long int nf, DataTable *data, 
			long long int nrhs);
void calcPower2(DataTable *data);
void calcPortPower2(DataTable *data, long long int ip, 
		    long long int ip_in, long long int mode1, 
		    long long int mode2, long long int jj);

/*   bicgstab.c   */
void makeCGmatrixTable(ELEMENT *element, long long int ne, long long int np, CGtable *cg_table,
		long long int NK);
long long int columnInCGmatrix(long long int ii, long long int jj, CGtable *cg_table, long long int np);

/*   client.c   */
/*
 long long int main1(long long int argc, char **argv, long long int my_rank);
 long long int readMatrixData(C_DataTable *data);
 long long int readRightHandVector(C_DataTable *data);
 long long int solveMatrixEquationInClinet(C_DataTable *data);
 long long int writeSolutionVector(C_DataTable *data);
 long long int readMessageFromServer(C_DataTable *data);
 */

/*   domain.c   */
/*
 long long int ddm(DataTable *data, std::complex<double> **A, std::complex<double> *b,
 CGtable *cg_table, long long int nr, std::complex<double> *x, long long int nd, long long int *sn);
 */

/*   elementMatrix.c   */
/*
 void rotation(std::complex<double> eps[3][3], std::complex<double> epsR[3][3],
 double angle);
 */
void InverseMatrix(double mu[3][3]);
long long int elementMatrix(DataTable *data, double *x, double *y, double *z, long long int matID,
		std::complex<double> sk[24][24], std::complex<double> sm[24][24]);

/*   elementMatrix2D.c   */
void create2dElementParam(double *x1, double *y1, double *b, double *c,
		double *ll, double *ss);
void elementMatrix2D(DataTable *data, double *x1, double *y1, long long int matID,
		std::complex<double> BB[N_EN][N_EN],
		std::complex<double> BBeh[N_EN][N_EN]);
void elementMatrix2DP(DataTable *data, double *x1, double *y1, long long int matID,
		std::complex<double> BB[N_EN][N_EN],
		std::complex<double> BBeh[N_EN][N_EN]);

/*   inputData.c   */
long long int checkArgument(DataTable *data, long long int argc, char **argv);
extern int checkArgument1(DataTable *data, long long int argc, char **argv);
void inputData(DataTable *data, char *name);

/*   intpol.c   */
void intpol(DataTable *data, long long int nn, long long int mode, 
	    long long int realflag);
// void setField(DataTable *data);
void addCoreframe(DataTable *data);
void detectElementForIntpol (DataTable *data);
void WritePPMFile(DataTable *data, long long int mode);
void setRGBTable(DataTable *data);

/*   intpol3D.c   */
void outputField(DataTable *data);
CXYZdouble intpol3D(DataTable *data, long long int n_element, 
		    long long int n_point);
CXYZdouble intpol2D(DataTable *data, long long int nn, double x0, double y0);
void portField(DataTable *data, long long int nn);
void outputField_for_Gid(DataTable *data);

/*   intpolYZ.c   */
void intpolYZ(DataTable *data, long long int nn);

/*   optfem3d.c   */

/*   output.c   */
void outputPortField(DataTable *data, int flag, int ii);
void intpolOutputPortField(DataTable *data);
void intpolAtPoint(DataTable *data, double x0, double y0, double z0,
		   std::complex<double> *Phix, std::complex<double> *Phiy, 
		   std::complex<double> *Phiz, 
		   std::complex<double> *Phix2, std::complex<double> *Phiy2, 
		   std::complex<double> *Phiz2, 
		   long long int mode, long long int nn);
void inputPortField(DataTable *data);


// void inputPeriodicField(DataTable *data, FILE *fp, PortInfo *port);
// void intpolPortField(DataTable *data, PortInfo *port);
// void getFieldValuePeriodic(DataTable *data, PortInfo *port, long long int l_matID,
// 		double x0, double y0, double z0, std::complex<double> *Psix,
// 		std::complex<double> *Psiy, std::complex<double> *Psiz, std::complex<
// 				double> *Phix, std::complex<double> *Phiy,
// 		std::complex<double> *Phiz, long long int flag);
// void createInputFieldInfo(DataTable *data);
// void getFieldValueFrom2DVFEM(FEM2D *fem2d, double xg, double yg, std::complex<
// 		double> *phix, std::complex<double> *phiy, std::complex<double> *phiz);

/*   power-new.c   */
/*
 void output2Dfield(DataTable *data, long long int n_port);
 void output2Dfield2(DataTable *data, long long int n_port);
 */

/*   power.c   */
void calcPower(DataTable *data);
void calcPortPower(DataTable *data, long long int ip, long long int ip_in, 
		   long long int mode1, long long int mode2, long long int jj);

/*   solveMatrixEquation.c   */
void readField(DataTable *data);
void solveMatrixEquation(DataTable *data);

/*   utility.c   */
void createS3d(DataTable *data);
void createS2d(DataTable *data, PortInfo *port, long long int n_port);
void sign2Dand3D(DataTable *data, PortInfo *port);

/*   callPARDISO.c   */
/*
 long long int callPARDISO(std::complex<double> **A, std::complex<double> *bb, long long int nn,
 std::complex<double> *phi, CGtable *cg_table);
 */

/*   wsmp   */
long long int zgsmp_LUS(std::complex<double> **A, std::complex<double> *bb, 
			long long int nn,
			std::complex<double> *phi, CGtable *cg_table);
long long int zssmp_LUS(std::complex<double> **A, std::complex<double> *bb, 
			long long int nn,
			std::complex<double> *phi, CGtable *cg_table);
long long int wsmp_pre(std::complex<double> **A, std::complex<double> *bb, 
		       long long int nf,
		       DataTable *data, long long int nrhs);
long long int pzgsmp_LUS(std::complex<double> **A, std::complex<double> *bb, 
			 long long int nf,
			 std::complex<double> *phi, CGtable *cg_table, 
			 long long int ns, long long int ne);
/*
long long int pzssmp_LUS(std::complex<double> **A, std::complex<double> *bb, long long int nf,
		std::complex<double> *phi, CGtable *cg_table, long long int ns, long long int ne);
*/
long long int pzssmp_LUS(std::complex<double> **A, std::complex<double> *bb, 
			 long long int nf,
			 std::complex<double> *phi, CGtable *cg_table, 
			 long long int ns, long long int ne);

void pwsmp_pre(std::complex<double> **A, std::complex<double> *bb, 
	       long long int nf,
	       DataTable *data, long long int ns, long long int ne, 
	       long long int nrhs);


/*   getlinek.c    */
/*
 long long int getlinek(char s[], const long long int lim);
 */

/* intpol_in_out.c */
long long int intpol_in_out(DataTable *data, long long int in_out_flag);


#ifdef PARDISO_FLOAT
void callPARDISO(long long int nn, long long int nc, long long int *ia, 
		 long long int *ja, std::complex<float> *avals, 
		 std::complex<float> *bb, std::complex<float> *phi, 
		 long long int *arguments, long long int nrhs);
#else
void callPARDISO(long long int nn, long long int nc, long long int *ia, 
		 long long int *ja, std::complex<double> *avals, 
		 std::complex<double> *bb, std::complex<double> *phi, 
		 long long int *arguments, long long int nrhs);
#endif

// mozaic.cpp
extern int inputMozaic(DataTable *data);
extern int checkMeshforMozaic(DataTable *data);
extern int checkMeshforeachVoxel(DataTable *data);
extern int shuffle(int *array, int size);

/*   port.c   */
void createPortInfo(DataTable *data);
void eachPort(DataTable *data, PortInfo *port, long long int n_port);
void savePortFile(DataTable *data, PortInfo *port, long long int n_port);
void normalizeInputPortPower(DataTable *data, PortInfo *port,
			      long long int ip, 
			      long long int mode);

// wfm.cpp
extern int DirectBinarySearch(DataTable *data, int kk);

// calcFitness.cpp
void calcFitness(DataTable *data, int flag, int total);

// random.c
void advance_random();
int flip(double prob);
void randomize();
double randomperc();
int rnd(int low, int high);
double rndreal(double lo, double hi);
void warmup_random(double random_seed);
double rand_normal( double mu, double sigma );
