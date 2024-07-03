#define MAX_COMPONENT  24
#define COMPLEX
#define DIM3
#define ELEMENT Element3D

#ifdef DIM3
typedef struct _ELEMENT {
    int kn[MAX_COMPONENT];
    int kk[MAX_COMPONENT];
	int pflag[MAX_COMPONENT];
    int matID;
} ELEMENT;
#else
typedef struct _ELEMENT {
    int kk[MAX_COMPONENT];
    int matID;
} ELEMENT;
#endif

typedef struct _CGtable {
	int n;
	int *col;
} CGtable;

#ifdef COMPLEX
#define DOUBLE_COMPLEX double_complex
#define CONJ(a)        (a)
#define ABS(a)         abs(a)
#else
#define DOUBLE_COMPLEX double
#define CONJ(a)        a
#define ABS(a)         fabs(a)
#endif

