/*   FEMmatrix.c   */
extern int matrix(DataTable *data, double *xx, double_complex *f,
		  double_complex sk[3][3], double_complex sm[3][3], int matID, 
		  double_complex *Eps, int element);
extern int reArrangment(DataTable *data);
extern int FEMmatrix(DataTable *data);

/* analysis.c */
extern int analysis(DataTable *data, int argc, char **argv);

/*   bhlu1c.c   */
extern int bhlu1c_(double_complex *a, double_complex *ac, int *n, 
	int *ml, int *mu, double *eps, double_complex *wk, int *ip, int *ier);
extern int bhslv1c_(double_complex *a, double_complex *ac, int *n,
	 int *ml, int *mu, double_complex *b, int *ip);

/* calcDispersion.c */
extern int calcDispersion(DataTable *data);

/*   calcNeff.c   */
extern int calcNeff(DataTable *data);

/*   eigenModeCalculation.c   */
extern int eigenModeCalculation(DataTable *data, int ii);

/*   eigv4p.c   */
extern int eigv4p_(double *xa, double *xm, int *n, 
	int *mu, int *ne0, int *ne, int *ialg, double *eval,
	double *eig, double *eigv, double *rn, double *wk, int *iwk, int *ier);

/*   inputData.c   */
extern int checkArgument(DataTable *data, int argc, char **argv);
extern int inputDataMsh(DataTable *data, char *fileName);
extern int MeshGeneration(DataTable *data);
extern int maxValue(int *a, int n);
extern int minValue(int *a, int n);
extern int checkBandWidth(DataTable *data);
extern int allocWorkSpace(DataTable *data);

/* intpol.c */
extern int intpol(DataTable *data, int ii, int number, double *field);

/*   InitialField.c   */
extern int Initialfield(DataTable *data);

/* material.c */
extern double material(DataTable *data, double Delta, double out[2]);

/*   output.c   */
extern int writeFieldDataR(DataTable *data, int ii, char *fileName);
extern int writeFieldDataC(DataTable *data, int ii, char *fileName);

/* power.c */
extern int NormalizeField(DataTable *data, double *field, double neff);
extern double normalize(DataTable *data, double *field, int flag);
extern int integr(DataTable *data, double xx[3], double_complex FF[3], 
                  int kk[3], double le, double NN, double *elementpower);

/*   propagation.c   */
extern int propagation(DataTable *data);

/*   utility.c   */
extern int matNN1Dsecond(double le, double ff[3][3]);
extern int matNyNy1Dsecond(double le, double ff[3][3]);
extern int matrNN(double r1, double r2, double le, 
				  double u, double v, double ff[3][3]);
extern int matrNrNr(double r1, double r2, double le, 
				    double u, double v, double ff[3][3]);
extern int matrNrN(double r1, double r2, double le, 
				    double u, double v, double ff[3][3]);
extern int matrinvNN(double r1, double r2, double le, 
				    double u, double v, double ff[3][3]);
extern int matrinvNN2(double r1, double r2, double le, double r0, 
				    double u, double v, double ff[3][3]);
extern int multifrontalFree();
