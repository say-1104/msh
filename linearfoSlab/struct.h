/*
	XZ 型
*/
typedef struct _XZdouble {
  double x, z;
} XZdouble;

typedef struct _XYdouble {
  double x, y;
} XYdouble;

/*
	諸パラメータ
*/
typedef struct _Param {

  XZdouble grid;
  int    *division, numOflayer;
  double *thickness;
  double wavelength;
  int numoflamda;
  double power;
  double lamda1, lamda2, step; 
  double k0, k02;
  double n0, n02;
  double **neff, **beta, *wl;
  double ***field, ***Gamma;
  double Xmin, Xmax, Zstart, Zend;

  int symFlag;
  int div;
  double dx;
  double *Ex, *Ey, *Ez, *Hx, *Hy, *Hz, *x;
  
  double x0;

} Param;

/*
	入力波形データ
*/
typedef struct _InputWave {
  int flag, eigv, number;
  char file[256];
  XYdouble center, width;
  int ne, np;
  int modeID, modeNo;
  double beta;
  Element *element;
  double *x0, *y0;
  double_complex *field;
  double c_neff, p_neff, c_beta;
  double xmin, ymin, xmax, ymax; 
} InputWave;

/*
	有限要素法データ
*/
typedef struct _FEMparam {
  int ne, np, nbw, nr;
  double_complex epsilon[20], epsilonInv[20];
  double rho[20], velocity[20];
  double *xx;
  double_complex *amat, *bmat;
  double_complex *field, *fieldn;
  int *matID;
#ifdef CG
  Element *element;
  CGtable *cg_table;
  double_complex **aa, **bb;
  double_complex **Pplus, **Pminus;
#else
  int **element;
  double *aa_half, *bb_half;
  double *aa, *bb, *vv;
  double_complex *Pplus, *Pminus;
  double_complex *cac, *cwk;
  int *iwp, ml, mu;
#endif
  int symFlag;
  double symAxisX, symAxisY;
  double xmin, xmax, ymin, ymax;
  int **table_n, ***table;
} FEMparam;

/*
	TBC 用のテーブル
*/
typedef struct _TBCtable {
  int el, ed;
} TBCtable;

typedef struct _PMLregion {
  double x0, y0, x1, y1, dx, dy;
} PMLregion;


typedef struct _PMLdata {
  int m;
  double tanD;
  int nx, ny;
  PMLregion x_data[N_PML], y_data[N_PML];
} PMLdata;

/*
	非線形導波路用データ
*/

typedef struct _Sampling {
  double_complex nn[7];
} Sampling;

typedef struct _Nonlinear {
  int cw_flag, flag, inputflag;
  double input, scale;
  double kerr[20], sat[20], tpa[20];
  Sampling *Eps;
/*
    int    bpmNP, bpmNE;
    int    element[2000][3];
    double *YY;
    double_complex *norm_phi;
    double Pin, Pout;
*/
} Nonlinear;

typedef struct _Probability {
/* CWの時に透過率，反射率を計算するのに使うパラメータ */
    int    bpmNP, bpmNE;
    int    element[2000][3];
    double *YY, *nn;
    double_complex *norm_phi, *diffz, *inputphi;
    double Pin, Pout;
    int    numOfref_y, numOfref_x;
    double ref_x[10], ref_y[10];
    double y1, y2, x1, x2;
} Probability;

/*
	データテーブル
*/
typedef struct _DataTable {

  int problem;
  char file[256];
  InputWave input;
  Param	par;
  FEMparam fem;
  PMLdata pml;
  Nonlinear nonlinear;
  Probability prob;

} DataTable;


