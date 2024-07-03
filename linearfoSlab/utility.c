#include "fem.h"

extern int matNN1Dsecond(double le, double ff[3][3])
{
	int		i, j, k;
	double	tmp1, tmp2;

	tmp1 = le/30.0;

	ff[0][0] = 4.0*tmp1;
	ff[0][1] = -1.0*tmp1;
	ff[0][2] = 2.0*tmp1;

	ff[1][0] = -1.0*tmp1;
	ff[1][1] = 4.0*tmp1;
	ff[1][2] = 2.0*tmp1;

	ff[2][0] = 2.0*tmp1;
	ff[2][1] = 2.0*tmp1;
	ff[2][2] = 16.0*tmp1;
	
}

extern int matNyNy1Dsecond(double le, double ff[3][3])
{
	int		i, j, k;
	double	tmp1, tmp2;

	tmp1 = 1.0/(3.0*le);

	ff[0][0] = 7.0*tmp1;
	ff[0][1] = 1.0*tmp1;
	ff[0][2] = -8.0*tmp1;

	ff[1][0] = 1.0*tmp1;
	ff[1][1] = 7.0*tmp1;
	ff[1][2] = -8.0*tmp1;

	ff[2][0] = -8.0*tmp1;
	ff[2][1] = -8.0*tmp1;
	ff[2][2] = 16.0*tmp1;
	
}

extern int matrNN(double r1, double r2, double le, 
				  double u, double v, double ff[3][3])
{
	int		i, j, k;
	double	tmp1, tmp2;

	tmp1 = r1*le/60.0;
	tmp2 = r2*le/60.0;

	ff[0][0] = tmp1*7.0+tmp2*1.0;
	ff[0][1] = tmp1*(-1.0)+tmp2*(-1.0);
	ff[0][2] = tmp1*4.0+tmp2*0.0;

	ff[1][0] = tmp1*(-1.0)+tmp2*(-1.0);
	ff[1][1] = tmp1*1.0+tmp2*7.0;
	ff[1][2] = tmp1*0.0+tmp2*4.0;

	ff[2][0] = tmp1*4.0+tmp2*0.0;
	ff[2][1] = tmp1*0.0+tmp2*4.0;
	ff[2][2] = tmp1*16.0+tmp2*16.0;
	
}

extern int matrNrNr(double r1, double r2, double le, 
				    double u, double v, double ff[3][3])
{
	int		i, j, k;
	double	tmp1, tmp2;

	tmp1 = r1/(le*6.0);
	tmp2 = r2/(le*6.0);

	ff[0][0] = tmp1*11.0+tmp2*3.0;
	ff[0][1] = tmp1*(1.0)+tmp2*(1.0);
	ff[0][2] = tmp1*(-12.0)+tmp2*(-4.0);

	ff[1][0] = tmp1*(1.0)+tmp2*(1.0);
	ff[1][1] = tmp1*3.0+tmp2*11.0;
	ff[1][2] = tmp1*(-4.0)+tmp2*(-12.0);

	ff[2][0] = tmp1*(-12.0)+tmp2*(-4.0);
	ff[2][1] = tmp1*(-4.0)+tmp2*(-12.0);
	ff[2][2] = tmp1*16.0+tmp2*16.0;
	
}

extern int matrNrN(double r1, double r2, double le, 
				    double u, double v, double ff[3][3])
{
	int		i, j, k;
	double	tmp1, tmp2;

	tmp1 = r1/30.0;
	tmp2 = r2/30.0;

	ff[0][0] = tmp1*(-13.0)+tmp2*(-2.0);
	ff[0][1] = tmp1*(2.0)+tmp2*(3.0);
	ff[0][2] = tmp1*(-14.0)+tmp2*(-6.0);

	ff[1][0] = tmp1*(-3.0)+tmp2*(-2.0);
	ff[1][1] = tmp1*2.0+tmp2*13.0;
	ff[1][2] = tmp1*(6.0)+tmp2*(14.0);

	ff[2][0] = tmp1*(16.0)+tmp2*(4.0);
	ff[2][1] = tmp1*(-4.0)+tmp2*(-16.0);
	ff[2][2] = tmp1*8.0+tmp2*(-8.0);
	
}

extern int matrinvNN1(double r1, double r2, double le, 
				    double u, double v, double ff[3][3])
{
	int		i, j, k;
	double A11, A12, A13;
	double A21, A22, A23;
	double A31, A32, A33;
	double	tmp1, tmp2;

	A11 = -4.0*u*u*u-10.0*u*u-(25.0/3.0)*u-5.0/2.0
	      +(u+1.0)*(u+1.0)*(2.0*u+1.0)*(2.0*u+1.0)*v;
	A12 = -4.0*u*u*u-6.0*u*u-(7.0/3.0)*u-1.0/6.0
	      +(u+1.0)*(2.0*u+1.0)*(2.0*u+1.0)*u*v;
	A13 = 2.0*u*u*u+4.0*u*u+(13.0/6.0)*u+1.0/6.0
	      -(u+1.0)*(u+1.0)*(2.0*u+1.0)*u*v;
	A22 = -4.0*u*u*u-2.0*u*u-(1.0/3.0)*u+1.0/6.0
	      +(2.0*u+1.0)*(2.0*u+1.0)*u*u*v;
	A23 = 2.0*u*u*u+2.0*u*u+(1.0/6.0)*u
	      -(u+1.0)*(2.0*u+1.0)*u*u*v;
	A33 = (-1.0)*u*u*u-(3.0/2.0)*u*u-(1.0/3.0)*u+1.0/12.0
	      +(u+1.0)*(u+1.0)*u*u*v;

	ff[0][0] = A11;
	ff[0][1] = A12;
	ff[0][2] = 4.0*A13;

	ff[1][0] = A12;
	ff[1][1] = A22;
	ff[1][2] = 4.0*A23;

	ff[2][0] = 4.0*A13;
	ff[2][1] = 4.0*A23;
	ff[2][2] = 16.0*A33;
	
}

extern int matrinvNN2(double r1, double r2, double le, double r0, 
				    double u, double v, double ff[3][3])
{
	int		i, j, k;
	double	tmp1, tmp2;

	tmp1 = le/(30.0*r0);

	ff[0][0] = 4.0*tmp1;
	ff[0][1] = -1.0*tmp1;
	ff[0][2] = 2.0*tmp1;

	ff[1][0] = -1.0*tmp1;
	ff[1][1] = 4.0*tmp1;
	ff[1][2] = 2.0*tmp1;

	ff[2][0] = 2.0*tmp1;
	ff[2][1] = 2.0*tmp1;
	ff[2][2] = 16.0*tmp1;
	
}
