#include "optfem.h"

#define COMMON
#define MODE_NUMBER 10

#define PRINT_BETA(a,b)	(a > 0.0) ? \
printf("(%18.16lf, %18.16lf)\n", sqrt(a)/(2.0*PI/b), 0.0) : \
printf("(%18.16lf, %18.16lf)\n", 0.0, -sqrt(fabs(a))/(2.0*PI/b))

// *added by muratsubaki what if the return statement is not processed?
extern int getMaterialIdOfNode(int nn, DataTable *data, int ii)
{
	int i, j;
	FEM *fem = &(data->fem);
	Param *par = &(data->par);
	double xg, yg, x1, x2, x3, y1, y2, y3, temp;

	for (i = 0; i < fem->ne; i++) {
		x1 = data->fem.node[fem->element[i].kk[0]].x;
		x2 = data->fem.node[fem->element[i].kk[1]].x;
		x3 = data->fem.node[fem->element[i].kk[2]].x;
		y1 = data->fem.node[fem->element[i].kk[0]].y;
		y2 = data->fem.node[fem->element[i].kk[1]].y;
		y3 = data->fem.node[fem->element[i].kk[2]].y;

		xg = (x1+x2+x3)/3.0;
		yg = (y1+y2+y3)/3.0;

		if (ii >= 0 && ii < data->par.numx){
			temp = data->par.center_y[ii]-yg;
		}
		else {
			temp = data->par.center_x[ii-data->par.numx]-xg;
		}

		for (j = 3; j < 6; j++) {
			if (fem->element[i].kk[j] == nn) {
				if (temp > 0.0) {
					return fem->element[i].matID;
				}
			}
		}
	}
}

// ii : port number
extern int portInfo(PortData *port, DataTable *data, int ii)
{
	int i;
	int layNo = -1;
	int m_id = -1, this_id;
	XYdouble node1, node2;
	double width;
		
	for (i = 0; i < MAXLAYER; i++) port->width[i] = 0.0;

	for (i = 0; i < port->np-2; i += 2) {
		this_id = getMaterialIdOfNode(port->kp[i+1], data, ii);

		node1 = data->fem.node[port->kp[i]];
		node2 = data->fem.node[port->kp[i+2]];

		if (this_id != m_id) {
			m_id = this_id;
			layNo++;
			port->muT[layNo] = data->par.muT[m_id];
			port->epT[layNo] = data->par.epT[m_id];
			port->width[layNo] += width = LENGTH(node1, node2);
			port->nn[layNo] = i;
		}
		else {
			port->width[layNo] += width = LENGTH(node1, node2);
		}

	}
	fprintf(stderr, "Number of layer = %d\n", layNo);
	if (port->normal.x == 1) {
		fprintf(stderr, " y = %lf , %lf\n", port->yp[0], port->yp[port->np-1]);
	}
	else {
		fprintf(stderr, " x = %lf , %lf\n", port->xp[0], port->xp[port->np-1]);
	}

	for (i = 0; i < layNo+1; i++) { 
		port->nn[i] = port->nn[i+1];
	}

	port->nn[layNo] = port->np;
	port->n_layer = layNo+1;

	for (i = 0; i < port->n_layer; i++) {
		if (data->par.modeID == 2) {
			port->pp[i] = 1.0/real(port->epT[i].zz);
			port->qq[i] = real(port->muT[i].zz);
		} else {
			port->pp[i] = 1.0/real(port->muT[i].zz);
			port->qq[i] = real(port->epT[i].zz);
		}
	}

	for (i = 0; i < port->n_layer; i++) { 
		fprintf(stderr, "Port->nn %d, er %lf\n", port->nn[i], port->pp[i]);
	}
}

extern int setPort(DataTable *data, double wl, int ll, int mm) {
	int i;
	for (i=0; i<data->port.number; i++) {
		data->port.data[i] = &(data->port.dataMemory[ll][mm][i]);
	}
	fprintf(stderr, "Port has been set appropriately");
}

extern int eigenModeOfPort(DataTable *data, double wl, int ll, int mm)
{
	int i;
	Port *port = &(data->port);
	int inputPort = data->par.inputPort;
	
	if (port->solver == INPUTEIGEN) {
		for (i = 0; i < port->number; i++) {
			eigenModeOfPortByFileInput(data, &(port->dataMemory[ll][mm][i]), wl, ll, mm, i);
		}
	}	
}

/* -------------------- Read eigenmode solver Solution ------------------ */
extern int eigenModeOfPortByFileInput(DataTable *data, PortData *port,
							double wl, int ll, int mm, int ii)
{
	int i, j, k;
	int div;
	double re, im;
	double beta, *xx, pos, maximum;
	double_complex *ff, slope;
	double xmin, tmp = 0.0;
	double neff, k0, ZZ, rr;
	FILE *fp;
	char string[256];
	double_complex cj(0.0, 1.0);

	// read input field
	if (data->par.num_mode != 1) {
		if (ii == data->par.inputPort || data->par.modEIMflag == 1) {
			sprintf(string, "./Input/InputField-%d_%d-%1.3f", ii, mm, wl);
		}
		else {
			sprintf(string, "./Input/InputField-%d-%1.3f", ii, wl);
		}
	}
	else {
		sprintf(string, "./Input/InputField-%d-%1.3f", ii, wl);
	}

	if ((fp = fopen(string, "r")) == NULL) {
		fprintf(stderr, "can't open file (%s)\n", string);
		exit(0);
	}

	fscanf(fp, "%*s %d %lf", &(div), &(beta));

	ALLOCATION(xx, double, div);
	ALLOCATION(ff, double_complex, div);

	for(i = 0; i < div; i++) {
		fscanf(fp, "%lf %lf", &(xx[i]), &re);
		ff[i] = re+cj*0.0;
	}

	xmin = xx[0];

	// beta2[0] is for mode number, not port number 
	port->beta[ll] = beta;
	port->beta2[ll] = beta*beta;

	fclose(fp);

	// interpolation field value
	for (i = 0; i < port->np; i++) port->PhiR[ll][i] = 0.0;
	maximum = 0.0;

	for(i = 0; i < port->np; i++) {
		if (port->normal.x == 1) pos = port->yp[i];
		else pos = port->xp[i];

		for(j = 0; j < div-1; j++) {
			if (fabs(pos-xx[j]) < 1.0e-9) {
				port->PhiR[ll][i] = ff[j];
				break;
			}
			else if (fabs(pos-xx[j+1]) < 1.0e-9) {
				port->PhiR[ll][i] = ff[j+1];
				break;
			}
			else if (pos > xx[j] && pos < xx[j+1]) {
				slope = (ff[j+1]-ff[j])/(xx[j+1]-xx[j]);
				port->PhiR[ll][i] = slope*(pos-xx[j])+ff[j];
				break;
			}
		}
		if (fabs(real(port->PhiR[ll][i])) > fabs(maximum)) maximum = real(port->PhiR[ll][i]);
	}

	// maximum is used for sign modification
	if (maximum < 0) maximum = -1.0;
	else maximum = 1.0;

	tmp = normalize(data, port, ll, ii);

	ZZ = 376.730313461;
	k0 = 2.0 * M_PI / wl;
	if (data->par.modeID == 2) {
		rr = k0 / ZZ / beta;
	} else {
		rr = k0 * ZZ / beta;
	}
			
	tmp /= rr * 2.0;
	tmp = sqrt(tmp) * maximum;
	for (i = 0; i < port->np; i++) {
		port->PhiR[ll][i] /= tmp;
	}

	sprintf(string, "Output/portfield\0");
	if (ii == 0 && ll == 0) {
		fp = fopen(string, "w");
	}
	else {
		fp = fopen(string, "a+");
	}

	if (port->normal.x == 1) {
		for (i = 0; i < port->np; i++) {
		fprintf(fp, "%lf %lf\n", port->yp[i], real(port->PhiR[ll][i]));
		}
	}
	else {
		for (i = 0; i < port->np; i++) {
		fprintf(fp, "%lf %lf\n", port->xp[i], real(port->PhiR[ll][i]));
		}
	}
	fprintf(fp, "\n");
	fclose(fp);

	for (i = 0; i < port->np; i++) {
		port->PhiL[ll][i] = port->PhiR[ll][i];
		port->PsiR[ll][i] = port->PhiR[ll][i];
		port->PsiL[ll][i] = port->PhiR[ll][i];
		port->FinR[ll][i] = port->PhiR[ll][i];
		port->FinL[ll][i] = port->PhiR[ll][i];
	}

	free(xx); free(ff);
}

/* -------------------- Pointsource ------------------ */
extern int Pointsource(DataTable *data)
{
	int i;
	// check x0, y0 exist in msh 
	data->port.centerkp = -1;
	for (i = 0; i < data->fem.np; i++) {
		if (fabs(data->port.x0-data->fem.node[i].x) < 1.0e-9 &&
			 fabs(data->port.y0-data->fem.node[i].y) < 1.0e-9) {
			data->port.centerkp = i;
			break;
		}
	}

	if (data->port.centerkp == -1) {
		fprintf(stderr, "Point x0, y0 do not exist!\n");
		exit(0);
	}
}


/* -------------------- Sheetsource ------------------ */

extern int Sheetsource(DataTable *data)
{
	int i;

	// check orientation
	// xflag = 1, source normal to x direction
	data->port.centerkp = -1;

	if (data->port.y0 == 1e50) {
		data->port.xflag = 1;
	}
	else {
		data->port.xflag = 0;
	}

	for (i = 0; i < data->fem.np; i++) {
		if (data->port.xflag == 1) {
			if (fabs(data->port.x0-data->fem.node[i].x) < 1.0e-9) {
				data->port.centerkp = i;
				break;
			}
		}
		else {
			if (fabs(data->port.y0-data->fem.node[i].y) < 1.0e-9) {
				data->port.centerkp = i;
				break;
			}
		}
	}

	if (data->port.centerkp == -1) {
		fprintf(stderr, "Point x0, y0 do not exist!\n");
		exit(0);
	}
}

// function for calulating scaling factor
extern double normalize(DataTable *data, PortData *port, int ll, int ii)
{
	int i, j, k, n = 0;
	int matID, ne;
	double pos[3], le;
	double_complex FF[3];
	double pp;
	double power;
	double_complex elementpower;

	power = 0.0;
	ne = port->ne;
	
	k = 1;
	n = 0;
	for(i = 0; i < ne; i++) {
		elementpower = 0.0;

		if (k-1 == port->nn[n]) n++;

		if (port->normal.x == 1) {
			pos[0] = port->yp[k-1];
			pos[1] = port->yp[k+1];
			pos[2] = port->yp[k];
		}
		else {
			pos[0] = port->xp[k-1];
			pos[1] = port->xp[k+1];
			pos[2] = port->xp[k];
		}

		FF[0] = port->PhiR[ll][k-1];
		FF[1] = port->PhiR[ll][k+1];
		FF[2] = port->PhiR[ll][k];

		le = pos[1]-pos[0];

		pp = port->pp[n];
		integral_elem(pos, FF, FF, le, pp, &elementpower);		
 
		power += real(elementpower);	// imagは0なので無視して良い

		k += 2;
	}
	
	fprintf(stderr, "power = %e\n", power);
	
	return(power);
}