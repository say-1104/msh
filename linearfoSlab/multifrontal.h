#include "./iteration.h"

extern int makeCGmatrixTable(ELEMENT *element, int ne, int np,
                            CGtable *cg_table, int NK);
extern int columnInCGmatrix(int ii, int jj, CGtable *cg_table, int np);
extern int multifrontal(double_complex **A, double_complex *b,
                        int nr, double_complex *x, CGtable *cg_table);
extern int multifrontalF(double_complex **A, double_complex *b,
                        int nr, double_complex *x, CGtable *cg_table);
extern int multifrontalS(double_complex **A, double_complex *b,
                        int nr, double_complex *x, CGtable *cg_table);
extern int multifrontalRF(double **A, double *b,
			int nr, double *x, CGtable *cg_table);
extern int multifrontalRS(double **A, double *b,
			int nr, double *x, CGtable *cg_table);
