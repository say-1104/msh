#define ALLOCATION(a,c,n)  if (n == 0) {				\
  } else {								\
    if ((a = (c *)calloc((n),sizeof(c))) == NULL) {			\
      fprintf(stderr, "can't allocate memory\n");			\
      exit(EXIT_FAILURE);						\
    }									\
  }

#define ZERO 1.0e-6

#define ROT1(x,z,cosA,sinA) (x*cosA-z*sinA)


#define MAX_CLIENT 25
#define PORT       9000
#define MAX_SIZE_1 5000000 //500000

// #define SVEA

#define CG
// #define LU

#define N_NODE 6
#define N_EDGE 9
#define N_EN (N_NODE+N_EDGE)

#define NORMAL 0
#define PERIODIC 1
#define PERIODIC_OTHER 2

#define MAX_MATERIAL 50

#define PPMCOLOR 64
#define PICTRATIO 100.0

#define MASTER_WORK 100

// #define SIMPLE

// added by sato
// #define PARDISO_FLOAT
