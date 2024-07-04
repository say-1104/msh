typedef struct _XYZdouble {
  double x, y, z;
} XYZdouble;

typedef struct _XYdouble {
  double x, y;
} XYdouble;

typedef struct _CXYZdouble {
  std::complex<double> x, y, z;
} CXYZdouble;

typedef struct Element2D {
  long long int kn[20];
  long long int ke[20];
  long long int matID;
} Element2D;

typedef struct Element2DS {
  long long int kk[15];
  long long int matID;
} Element2DS;

typedef struct _Matrix {
  long long int nn, nc;
  long long int *ia, *ja;
  std::complex<double> *avals, *b;
} Matrix;

typedef struct _FEM {
  long long int unknown;
  long long int ne, np, n_edge, nr_edge, nbw;
  long long int order, n_en, n_node, n_node_use;
  long long int n_port;
  XYZdouble min, max;
  Element3D *element;
  XYZdouble *node;
  std::complex<double> *phi;
  std::complex<double> *phi_in;
  long long int *inFlag;
  long long int *s3d;
  /* for CG */
  CGtable *cg_table;
  std::complex<double> **G;
  long long int **Gtable;
  long long int nr;
  /* ------ */
  std::complex<double> **ok;
  std::complex<double> **ok12;
  long long int *ip;
  
  // added by satonori for makeTable
  long long int **table_n;
  long long int ***table;
} FEM;

typedef struct _FEM2D {
  long long int nt, nz, np, nr, ne, order;
  long long int l_edge, l_node, l_en;
  Element2DS *element;
  XYdouble *node;
  long long int *iz, *izz;
  std::complex<double> *vv;
  double beta;
} FEM2D;

typedef struct _Mode {
  long long int target;
} Mode;

typedef struct _RGBcolor {
	int r, g, b;
} RGBcolor;

typedef struct _Picture {
	long long int nx, ny, nz;
	XYZdouble **absxz;
	XYZdouble **realxz;
	XYZdouble **str;
	RGBcolor rgb[PPMCOLOR + 1];
	RGBcolor rgb2[PPMCOLOR + 1];
	long long int **element, **f_flag;
	long long int ecount;

	// 入出力導波路とPMLの境界にてノイズが入りやすいので、更にPMLの厚さ分だけ内側の枠内で、フィールド分布の最高値・最低値の規格化を行う
	bool modifiedScale;
	bool originalSize;
} Picture;

// 波長ごとに屈折率を設定するのが面倒なので
typedef class _Param {
public:
	double wavelength, k0, k02;
	long long int n_material;
	double er[MAX_MATERIAL], mr[MAX_MATERIAL];
	long long int PMLflag[MAX_MATERIAL];
	double PMLangle[MAX_MATERIAL];
	double PMLz0[MAX_MATERIAL];
	double PMLx0[MAX_MATERIAL];
	double PMLthick[MAX_MATERIAL];
	double ermatrix[MAX_MATERIAL][3][3], mrmatrix[MAX_MATERIAL][3][3];
	Mode mode;
	double INPUT_Z;
	double out_slv_x_cood;
	double out_slv_y_cood;
	double out_slv_z_cood;
	long long int inputFlag;
	long long int pml_type;
	double neff, neff_out;
	double power_adjust;
	char in_file[256];
	int nofield;
	
	// added by satonori for makeTable
	double x1, x2;
	double y1, y2;
	double z1, z2;
	long long int I, J;
	long long int *mat;
	long long int div_max;

	int numofwlMPI;
	int myid, num_proc;
	int numStructures;
	int strID, wlID;
} Param;


typedef struct _PortInfo {
  int modenum;
  long long int ne, np, n_edge, nr;
  long long int l_edge, l_node, l_order, l_en;
  double x, y, z, angle;
  Element2D *element;
  long long int *toG, *toLn, *toLe;
  long long int *s2d, *s23d;
  std::complex<double> **m_func;
  std::complex<double> *beta;
  std::complex<double> *vv;
  double power, amp, phase;
} PortInfo;


typedef struct _PMLdata {
  double m;
  double tanD;
  double dx1, dx2, dy1, dy2, dz1, dz2;
} PMLdata;

typedef struct _Mozaic {
  int  Nx, Ny, Nz; // The number of mozaic in y and z directions
  double zs, ys, xs; // left lower x,y coordinate of the mozaic region
  double dx, dy, dz; // the length of one mozaic region in  y and z directions

  int symflag; // symmetry flag 0: no symmetry, 1: symmetry in y direction
  int totalN; // The number of mozaic pattern to be generated
  
  int mozaicmatID; // matID for mozaic region 
  int coreID, cladID; // matID used for core and clad in mozaic region
  int ***BorW; // 2D pattern of mozaic, 0: white, 1: black
  // BorW[ii][i][j]:ii No, i z-direction, j y-direction 

  double **z0, **z1, **y0, **y1, **x0, **x1;
  // minimum and max z and y coordinate of each mozaic

  int *matID; // matID of all the elements, used in 2D-FEM
  int *ii_pixel, *jj_pixel; // (ii, jj) pixel for all elements
  // elements outside the mozaic region, ii = jj = -1 

  int mozflag; // 1: do mozaic simulation
  int mozDev;

  int **ne_pixel, ***kk_pixel; 
  // the number of element and element number information of pixel

  int mzcalcNo;
  int DBSitr, DBSflag;
  double efficiency, gamma;
  double termination;
  double *efficiency_history;
  double **TT;
  int PixelReduction;
  std::complex<double> **c_TT;

  int numofT;
  int starti, startj, starttotal;

  char filename[256];


  double ***eta_mozaicMPI;

  double fittemp;
  int numModes;

} Mozaic;

typedef struct _DataTable {
  char *name;
  Param par;
  Picture pict;
  FEM fem;
  PortInfo *port;
  Matrix matrix;
  PMLdata pml;
  char pml_setting_2D_output[16];
  char pml_setting_2D_input[16];
  double pml_setting_2D_thickness;
  long long int nm;
  long long int np;
  long long int n_port_io;
  double port_xmin[16];
  double port_xmax[16];
  double port_zmin[16];
  double port_zmax[16];
  
  long long int mode;
  long long int arguments[64];
  char filename[128];
  long long int nrhs;
#ifdef PARDISO_FLOAT
  std::complex<float> *avals;
  std::complex<float> *bb;
  std::complex<float> *x;
  std::complex<float> *phi;
#else
  std::complex<double> *avals;
  std::complex<double> *bb;
  std::complex<double> *x;
  std::complex<double> *phi;
#endif

  long long int calcflag;

  Mozaic mozaic;
} DataTable;

typedef struct _C_DataTable {
  long long int nr;
  std::complex<double> **A;
  CGtable *cg_table;
  std::complex<double> *b, *x;
  long long int soc;
} C_DataTable;
