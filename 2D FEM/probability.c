#include "optfem.h"

/* ----------------- extract field on port from FEM results --------------- */
extern int portField(DataTable *data, int ll)
{
	int ip, i;
	PortData *port;
	
	for (ip = 0; ip < data->port.number; ip++) {
		port = data->port.data[ip];
		for (i = 0; i < port->np; i++) {
			port->phi[i] = data->fem.phi[port->kp[i]];
		}
	}
	
	if (data->port.inputDirection == BOTH_DIRECTION) {
		port = data->port.data[data->par.inputPort];
		for (i = 0; i < port->np; i++) {
			port->phi[i] -= port->PhiR[ll][i];
		}
	}
}

extern int probability(DataTable *data, double wl, int ll, int mm)
{
	int i, k, ip, im, n_mode;
	double_complex amp, amp_n;
	char string[256];
	FILE *fp;
	double_complex beta_in, beta_out;

	n_mode = data->par.num_mode;
	
	for (ip = 0; ip < data->port.number; ip++) {
		beta_in = sqrt(data->port.data[data->par.inputPort]->beta2[ll]);
		beta_out = sqrt(data->port.data[ip]->beta2[ll]);

		if (data->port.data[ip]->beta2[ll] > 0.0) {
			amp = overlapIntegral(data->port.data[ip], data->port.data[ip], ll);
			amp_n = normalize(data, data->port.data[ip], ll, ip);
			amp /= amp_n;
			data->mosaic.TT[ip][ll*n_mode+mm] = pow(std::abs(amp), 2.0);
			data->mosaic.amp[ip][ll*n_mode+mm] = amp;
			data->mosaic.arg[ip][ll*n_mode+mm] = atan2(std::imag(amp), std::real(amp));
		}
	}
}

// function for calulating scaling factor 
// *added by muratsubaki
// Changed name from overlap_fujisawa into overlapIntegral
// this function calculates the overlap integral of
// InputField (of port1) and currentField (of port2)
extern double_complex overlapIntegral(PortData *port1, PortData *port2, int m)
{
	// accord with the direction of the port
	if (port1->normal.x*port2->normal.x 
			+ port1->normal.y*port2->normal.y != 1.0) {
		fprintf(stderr, "Unable to calculate overlap integral since directions of two ports aren't in agreement\n");
		exit(1);
	}

	int i, j, k, n = 0;
	int matID, ne;
	double pos[3], le;
	double_complex FF1[3], FF2[3];
	double pp;
	double_complex amp, elem_amp;

	amp = 0.0;
	ne = port1->ne;
	
	k = 1;
	n = 0;
	for(i = 0; i < ne; i++) {
		elem_amp = 0.0;

		if (k-1 == port1->nn[n]) n++;

		if (port1->normal.x == 1.0) {
			pos[0] = port1->yp[k-1];
			pos[1] = port1->yp[k+1];
			pos[2] = port1->yp[k];
		}
		else {
			pos[0] = port1->xp[k-1];
			pos[1] = port1->xp[k+1];
			pos[2] = port1->xp[k];
		}

		FF1[0] = port1->PhiR[m][k-1];
		FF1[1] = port1->PhiR[m][k+1];
		FF1[2] = port1->PhiR[m][k];

		FF2[0] = port2->phi[k-1];
		FF2[1] = port2->phi[k+1];
		FF2[2] = port2->phi[k];

		le = pos[1]-pos[0];

		pp = port1->pp[n];
	
		integral_elem(pos, FF1, FF2, le, pp, &elem_amp);
 
		amp += elem_amp;
		k += 2;
	}  
	return(amp);
}

#define GZAI1  -0.90617985
#define GZAI2  -0.53846931
#define GZAI3  0.0
#define GZAI4  0.53846931
#define GZAI5  0.90617985
 
extern int integral_elem(double xx[3], 
			 double_complex FF1[3], double_complex FF2[3], 
			 double le, double pp, double_complex *elem_amp)
{
	int i,j,k;
	double  gzai, w, L1, L2;
	double  N[3];
	double  Ex_square;
	double_complex phi1, phi2;

	for(k = 0; k < 5; k++) {
		switch(k) {
		case 0:   
			gzai = GZAI1; w = 0.23692689; break;
		case 1:   
			gzai = GZAI2; w = 0.47862867; break;
		case 2:   
			gzai = GZAI3; w = 0.56888889; break;
		case 3:   
			gzai = GZAI4; w = 0.47862867; break;
		case 4:   
			gzai = GZAI5; w = 0.23692689; break;
		}
		
		L1 = (1-gzai)/2;  L2 = (1+gzai)/2;
		
		N[0] = L1*(2*L1-1);
		N[1] = L2*(2*L2-1);
		N[2] = 4*L1*L2;

		phi1 = FF1[0]*N[0]+FF1[1]*N[1]+FF1[2]*N[2];
		phi2 = FF2[0]*N[0]+FF2[1]*N[1]+FF2[2]*N[2];
		
		*elem_amp += 0.5*le*w*pp*conj(phi1)*phi2;
	}
}

extern int outputOverlap(DataTable *data, int kk)
{
	int i, k, ip, ll, mm, n_mode;
	double lambda;
	char string[256];
	FILE *fp;

	n_mode = data->par.num_mode;
	
	sprintf(string, "Output/Overlap");
	if (kk == 0) {
		if ((fp = fopen(string, "w")) == NULL) {
			fprintf(stderr, "can't open file (%s)\n", string);
			return 0;
		}
	}
	else {
		if ((fp = fopen(string, "a+")) == NULL) {
			fprintf(stderr, "can't open file (%s)\n", string);
			return 0;
		}
	}
	
	fprintf(fp, "wl\tInputMode\tOutputMode\tTransmission\tAmplitude\tPhase[rad]\tRe\tIm\n");
	
	for (ll = 0; ll < data->par.numoflambda; ll++) {
		for (mm = 0; mm < n_mode; mm++) {
			for (ip = 0; ip < data->port.number; ip++) {
				lambda = data->par.w_start+(double)ll*data->par.w_step;
				if (ip == data->par.inputPort) fprintf(fp, "%lf\t%d\tref\t", lambda, mm);
				else fprintf(fp, "%lf\t%d\t%d\t", lambda, mm, ip);
				
				fprintf(fp, "%.15lf\t%.15lf\t%.15lf\t%.15lf\t%.15lf\n", data->mosaic.TT[ip][ll*n_mode + mm], 
					sqrt(data->mosaic.TT[ip][ll*n_mode + mm]), data->mosaic.arg[ip][ll*n_mode + mm],
					imag(data->mosaic.amp[ip][ll*n_mode + mm]), real(data->mosaic.amp[ip][ll*n_mode + mm]));
			}
		}
	}

	fclose(fp);
}
