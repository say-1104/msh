#define MAXMODE         200
#define MAXNP           2001

#define MAXLAYER        100
#define MAX_MATERIAL	50

#define N_PML	100

typedef struct Element3D {
	int kn[MAX_COMPONENT];
	int kk[MAX_COMPONENT];
	int pflag[MAX_COMPONENT];
	int matID;
} Element3D;

typedef struct _CGtable {
	int n;
	int *col;
} CGtable;

/* �󼡸���ɸ */
typedef struct _XYdouble {
	double x, y;
} XYdouble;

/* �󼡸���ɸ��ʣ�ǿ��� */
typedef struct _XYCZdouble {
	double x, y;
	double_complex z;
} XYCZdouble;

typedef struct _Tensor {
	double_complex xx, xy, xz;
	double_complex yx, yy, yz;
	double_complex zx, zy, zz;
} Tensor;

/* ���Ϥ��Ѥ���ѥ�᡼�� */
typedef struct _Param {
	double w_start, w_end, w_step; // ��Ĺ(����, ����, ���ƥå�)
	int numoflambda;
	double wavelength;			   // ��Ĺ
	double k0, k02;
	double power;
	int    modeID;                 // TE = 1 or TM = 2
	int modeNo;
	int    f_step;                 // ������ (���֤�̤����)
	int    n_material;             // ������ (number of materials)
	double er[MAX_MATERIAL], mr[MAX_MATERIAL]; // Ͷ��Ψ, Ʃ��Ψ
	Tensor epT[MAX_MATERIAL];      // �ƥ󥽥�Ͷ��Ψ
	Tensor muT[MAX_MATERIAL];      // �ƥ󥽥�Ʃ��Ψ
	int    PMLflag[MAX_MATERIAL];      // PML ���ֹ� (0 ���̾�β����ΰ�)
	double PMLthick[2][MAX_MATERIAL];  // PML �θ���
	double PMLangle[MAX_MATERIAL];  // PML �θ����Ƥ������
	double PMLx0[MAX_MATERIAL];     // PML �������� x ��ɸ
	double PMLy0[MAX_MATERIAL];     // PML �������� y ��ɸ
	double reflect;                 // ����ȿ�ͷ����� log10
	int    PMLtype;
	double pml_thick_x, pml_thick_y;
	int    inputPort;
	int    num_mode;
	int    modEIMflag;
	
	int    mo_effect; 

	// field related parameters
	int div_x, div_y, numx, numy;
	double *center_x, *center_y;
	int division_x, division_y;

	// added for TD-BPM msh
	int inflag[MAX_MATERIAL];
	int symFlag;
	double symAxisX, symAxisY;

} Param;

/* ͭ������ˡ���Ϥ��Ѥ��������� */
typedef struct _FEM {
	int ne;         // ���ǿ�
	int np;         // ������
	int nr;         // �ºݤ˲�������
	int nbw;        // ����ΥХ����
	Element3D *element;   // ���Ǿ���
	XYdouble *node;      // ������ɸ
	double_complex *phi;   // ��٥��ȥ�
	double_complex *b;     // ���ե٥��ȥ�
	CGtable *cg_table;  // ������ʬ�γ�Ǽ���ξ���
	double_complex **A;        // ����ΰ��̷�����
	double_complex **A12;      // ���о�������
	double_complex **B12;      // ���о�������
	
	double xmax, xmin, ymax, ymin;  // ��ɸ�κ����͡��Ǿ���
	int **table_n, ***table;
	int *inFlag;    // ���������Ͷ����κ����ɤ���ˤ��뤫�Υե饰
	double_complex *FinR, *FinL;  // ��¦(��¦)���ͥ٥��ȥ�

	int symFlag;
} FEM;

/* �ġ��������ϥݡ��Ȥ˴ؤ���ǡ��� */
typedef struct _PortData {
	int     np;              // �ݡ��Ⱦ��������
	int     *kp;             // �ݡ��Ⱦ�����Ǿ���
	double  *xp;             // �ݡ��Ⱦ��������ɸ
	double  *yp;             // �ݡ��Ⱦ��������ɸ
	double  x0, y0;          // �ݡ��Ȥ��濴��ɸ
	int     n_layer;         // �ݡ��Ȥ����ؤ��޼��ǹ�������Ƥ��뤫
	int     nn[MAXLAYER];    // ���ؤ�ʬ���
	double  pp[MAXLAYER];    // ��ʬ�������η��� (������)
	double  qq[MAXLAYER];    // ��ʬ�������η��� (������)
	Tensor  epT[MAXLAYER];   // �ƥ󥽥�Ͷ��Ψ
	Tensor  muT[MAXLAYER];   // �ƥ󥽥�Ʃ��Ψ
	double  width[MAXLAYER]; // ���ؤθ���
	XYdouble normal;              // Normal Vector(indicate the direction of wave propagation)


	double	angle;
	
	/* for creating right vector */
	double_complex **NN, **NNy;          
	double_complex Amp;  // �ƥ⡼�ɤ����Ϳ���
	
	/* for calculating power */
	double_complex  *phi;   // �ݡ��Ⱦ�γ�ʬ��
	double_complex  FinR[MAXMODE][MAXNP];  // ���ͥ⡼�� (��)
	double_complex  FinL[MAXMODE][MAXNP];  // ���ͥ⡼�� (��)
	double_complex  PhiR[MAXMODE][MAXNP];  // ��ͭ�⡼�� (��)
	double_complex  PhiL[MAXMODE][MAXNP];  // ��ͭ�⡼�� (��)
	double_complex  PsiR[MAXMODE][MAXNP];  // �б�����⡼�ɴؿ� (��)
	double_complex  PsiL[MAXMODE][MAXNP];  // �б�����⡼�ɴؿ� (��)
	
	/* for analytic */
	double beta2[MAXMODE];     // ��������Σ���(���Хͥå���Ȥ��θ���뤿��
	double_complex beta[MAXMODE];      // ������� (����Ū�ˤϻȤ�ʤ�)

	char   file[256];          // ��ͭ�⡼�ɤ���¸����Ƥ���ե�����̾

	int ne;  
} PortData;

/* �����ϥݡ��Ȥ˴ؤ���ǡ��� */
typedef struct _Port {
	int      number;             // how many ports there are.
	int      inputDirection;     // SINGLE or BOTH
	double forward, backward;  // ���ͥѥ (���ȡ�����)
	SOLVER     solver;          // ���ͥ⡼�ɤλ���ˡ(���ϲ򡤣ƣţ�)
	PortData **data;              // Container of PortData
	PortData ***dataMemory;       // Stores instance of portdata

	// for point and sheet source
	int centerkp;
	double x0, y0;

	int xflag;
} Port;

/*  */
typedef struct _Nonlinear {
	double kerr[MAXLAYER], sat[MAXLAYER];
} Nonlinear;

typedef struct _PMLregion {
	double x0, y0, x1, y1, dx, dy;
} PMLregion;

typedef struct _PMLdata {
	int m;
	double tanD;
	int nx, ny;
	PMLregion x_data[N_PML], y_data[N_PML];

	// for cylindrical PML
	double thick, x0, y0, r0;

} PMLdata;

typedef struct _Mosaic {
	int  Ny, Nz; // The number of mosaic in y and z directions
	double zs, ys; // left lower x,y coordinate of the mosaic region
	double dy, dz; // the length of one mosaic region in  y and z directions

	int symflag; // symmetry flag 0: no symmetry, 1: symmetry in y direction

	int PixelReduction;
	
	
	int mosaicmatID; // matID for mosaic region 
	int coreID, cladID; // matID used for core and clad in mosaic region

	// 2D pattern of mosaic, 0: white, 1: black, -1 : not considered
	int **BorW; 	// Ny*Nz

	double **z0, **z1, **y0, **y1;
	// minimum and max z and y coordinate of each mosaic

	int *matID; // matID of all the elements, used in 2D-FEM
	int *ii_pixel, *jj_pixel; // (ii, jj) pixel for all elements
	// elements outside the mosaic region, ii = jj = -1 

	int mosflag; // 1: do mosaic simulation

	int **ne_pixel, ***kk_pixel; 
	// the number of element and element number information of pixel

	int mzcalcNo;


	int DBSflag;
	int DBSitr;
	int* DBSseq;
	double *FOM_history;
	double FOM;
	int evaluation;

	// TTは電力で、ampは複素振幅係数を表す。
	// クロストーク計算において、amp^2は、電力を表さない。
	double termination, gamma;
	double **TT;
	double_complex **amp;
	double **arg;

	double *zz, *yy;

	DEVICETYPE Device;
} Mosaic;

typedef struct _ModifiedEffectiveIndex {
	int size;
	std::vector<double> **er;         // index(0: core 1: cladding, 2:hole)
} ModifiedEffectiveIndex;

/* ���Ϥ��Ѥ�����ǡ��� */
typedef struct _DataTable {
	char file[256];   // �ե�����̾
	Param par;         // �ѥ�᡼��
	FEM fem;         // ͭ�����ǥ�å������ξ���
	Port port;        // �����ϥݡ��Ȥ˴ؤ������
	Nonlinear nonlinear;
	PMLdata     pml;
	long long int nm;
	long long int np;
	int realflag;
	Mosaic mosaic;
	ModifiedEffectiveIndex ModEI;
} DataTable;
