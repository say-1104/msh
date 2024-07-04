#define MAX_COMPONENT  24
#define COMPLEX
#define DIM3
#define ELEMENT Element3D


#ifdef DIM3
typedef struct _ELEMENT {
	long long int kn[MAX_COMPONENT];
	long long int kk[MAX_COMPONENT];
	long long int matID;
}ELEMENT;
#else
typedef struct _ELEMENT {
	long long int kk[MAX_COMPONENT];
	long long int matID;
}ELEMENT;
#endif

typedef struct _CGtable {
	long long int n;
	long long int *col;
} CGtable;
