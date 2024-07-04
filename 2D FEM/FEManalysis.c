// input field in portX isn't permitted(considered)
#include "optfem.h"

extern int elementMatrix(DataTable *data, double *x, double *y,
			double_complex ek[6][6], int matID, double wl)
{
	int i, j, k, ip, flag;
	double A = 1.0/3.0;
	double B = 0.05971587;
	double C = 0.47014206;
	double D = 0.79742669;
	double E = 0.10128651;
	double_complex sk[6][6], sm[6][6];
	double k0 = 2.0*PI/wl;
	double k02 = k0*k0;

	double L1, L2, L3, W;
	double N[6], NL1[6], NL2[6];
	double J11, J12, J21, J22;
	double Nx[6], Ny[6];
	double NxNx[6][6], NyNy[6][6], NN[6][6], NxNy[6][6], NyNx[6][6];
	double dj, coef1;
	double_complex Pxx, Pyy, Pxy, Pyx, Qzz;
	double xg, yg;
	double cosA, sinA, angle, x0, y0, xr0, yr0, xrg, yrg;
	double_complex co, delta, deltaInv;
	double ss, thick;
	double_complex cj(0.0, 1.0);
	Tensor ep, mu, epR, muR, Rot;
	double_complex sigma;
	double xmax, xmin, ymax, ymin;
	double thick_x, thick_y;

	/* PML */
	double s_max = 1.0/(2.0*data->par.k0)*log(1.0/data->pml.tanD);
	double_complex s0, sx, sy, sz;
	double_complex Sx, Sy, Sz;
	double dd, rho;
	double mm = data->pml.m;
	double zg, zmax, zmin;

	double_complex pp, qq;

	if (data->par.PMLtype == 0) {
	thick_x = data->par.pml_thick_x;
	thick_y = data->par.pml_thick_y;
	xmax = data->fem.xmax-thick_x;
	xmin = data->fem.xmin+thick_x;
	ymax = data->fem.ymax-thick_y;
	ymin = data->fem.ymin+thick_y;
	}

	if ( data->par.modeID == 1 ) {
		pp = 1.0;
		qq = data->par.er[matID];
	} else {
		pp = 1.0/data->par.er[matID];
		qq = 1.0;
	}

	for (i = 0; i < 6; i++) {
		for (j = 0; j < 6; j++) {
			sk[i][j] = 0.0; sm[i][j] = 0.0; ek[i][j] = 0.0;
		}
	}

	for (k = 0; k < 7; k++) {
		switch (k) {
			case 0: L1 = A; L2 = A; L3 = A; W = 0.225;      break;
			case 1: L1 = B; L2 = C; L3 = C; W = 0.13239415; break;
			case 2: L1 = C; L2 = B; L3 = C; W = 0.13239415; break;
			case 3: L1 = C; L2 = C; L3 = B; W = 0.13239415; break;
			case 4: L1 = D; L2 = E; L3 = E; W = 0.12593918; break;
			case 5: L1 = E; L2 = D; L3 = E; W = 0.12593918; break;
			case 6: L1 = E; L2 = E; L3 = D; W = 0.12593918; break;
		}

		N[0] = L1*(2.0*L1-1.0); NL1[0] = 4.0*L1-1.0; NL2[0] = 0.0;
		N[1] = L2*(2.0*L2-1.0); NL1[1] = 0.0;		 NL2[1] = 4.0*L2-1.0;
		N[2] = L3*(2.0*L3-1.0); NL1[2] = 1.0-4.0*L3; NL2[2] = 1.0-4.0*L3;
		N[3] = 4.0*L1*L2;       NL1[3] = 4.0*L2;     NL2[3] = 4.0*L1;
		N[4] = 4.0*L2*L3;       NL1[4] = -4.0*L2;    NL2[4] = 4.0*(L3-L2);
		N[5] = 4.0*L3*L1;       NL1[5] = 4.0*(L3-L1);NL2[5] = -4.0*L1;

		J11 = J12 = J21 = J22 = 0.0;
		for (i = 0; i < 6; i++) {
			J11 += NL1[i]*x[i];
			J12 += NL1[i]*y[i];
			J21 += NL2[i]*x[i];
			J22 += NL2[i]*y[i];
		}

		dj = J11*J22-J12*J21;
		xg = L1*x[0]+L2*x[1]+L3*x[2];
		yg = L1*y[0]+L2*y[1]+L3*y[2];

		// PML settings
		sx = sy = sz = 1.0;
		for (i = 0; i < data->pml.nx; i++) {
			if ((data->pml.x_data[i].x0-xg)*(data->pml.x_data[i].x1-xg) < 0.0 &&
				(data->pml.x_data[i].y0-yg)*(data->pml.x_data[i].y1-yg) < 0.0) {
				dd = fabs(data->pml.x_data[i].dx);
				rho = fabs(xg-data->pml.x_data[i].x0);
				s0 = 1.0-cj*s_max/dd*pow(rho/dd, mm)*(1.0+mm);
				sx = s0;
			}
		}
		for (i = 0; i < data->pml.ny; i++) {
			if ((data->pml.y_data[i].x0-xg)*(data->pml.y_data[i].x1-xg) < 0.0 &&
				(data->pml.y_data[i].y0-yg)*(data->pml.y_data[i].y1-yg) < 0.0) {
				dd = fabs(data->pml.y_data[i].dy);
				rho = fabs(yg-data->pml.y_data[i].y0);
				s0 = 1.0-cj*s_max/dd*pow(rho/dd, mm)*(1.0+mm);
				sy = s0;
			}
		}
		Sx = sy*sz/sx;
		Sy = sz*sx/sy;
		Sz = sx*sy/sz;

		for (i = 0; i < 6; i++) {
			Nx[i] = J22*NL1[i]-J12*NL2[i];
			Ny[i] = -J21*NL1[i]+J11*NL2[i];
		}

		for (i = 0; i < 6; i++) {
			for (j = 0; j < 6; j++) {
				NxNx[i][j] = Nx[i]*Nx[j];
				NyNy[i][j] = Ny[i]*Ny[j];
				NN[i][j]   = N[i]*N[j];
				NxNy[i][j] = Nx[i]*Ny[j];
				NyNx[i][j] = Ny[i]*Nx[j];
			}
		}

		ep = data->par.epT[matID];
		mu = data->par.muT[matID];

		for(i = 0; i < 6; i++) {
			for(j = 0; j < 6; j++) {
				sk[i][j] = W/dj*0.5*(pp*NyNy[i][j]/Sx+pp*NxNx[i][j]/Sy);
				sm[i][j] = W*dj*0.5*Sz*qq*NN[i][j];
			}
		}

		for (i = 0; i < 6; i++) {
			for (j = 0; j < 6; j++) {
				ek[i][j] += (sk[i][j]-k02*sm[i][j]);
			}
		}
	}
}

extern int Amatrix(DataTable *data, CGtable *cg_table, double wl)
{
	int     i, j, k;
	int     ii, jj, it;
	double  x[6], y[6];
	int     matID;
	FEM     *fem = &(data->fem);
	Param	*par = &(data->par);
	Port	*port = &(data->port);
	double_complex ek[6][6];
	double  xg, yg;
	double_complex cj(0.0, 1.0);
	double xin, yin, xr, yr, angle;
	int pmlFlag;

	for (i = 0; i < data->fem.nr; i++) {
		for (j = 0; j < cg_table[i].n; j++) {
			fem->A[i][j] = 0.0;
		}
	}

	for (k = 0; k < fem->ne; k++) {
		for (j = 0; j < 6; j++) {
			i = fem->element[k].kk[j];
			x[j] = fem->node[i].x;
			y[j] = fem->node[i].y;
		}
		matID = fem->element[k].matID;
		pmlFlag = data->par.PMLflag[matID];

		if (data->mosaic.mosflag == 1) {
			if (matID == data->mosaic.mosaicmatID) {
				ii = data->mosaic.ii_pixel[k];
				jj = data->mosaic.jj_pixel[k];

				if (data->mosaic.BorW[ii][jj] == 1) {
					matID = data->mosaic.coreID;
				}
				else if (data->mosaic.BorW[ii][jj] == 0 && data->par.modEIMflag != 2) {
					matID = data->mosaic.cladID;
				}
			}
		}

		elementMatrix(data, x, y, ek, matID, wl);

		for (i = 0; i < 6; i++) {
			for (j = 0; j < 6; j++) {
				ii = fem->element[k].kk[i];
				jj = fem->element[k].kk[j];
				if (ii < fem->nr && jj < fem->nr) {
					it = columnInCGmatrix(ii, jj, cg_table, data->fem.nr);
					fem->A[ii][it] += ek[i][j];
				}
			}
		}
	}
}

extern int Bmatrix(PortData *port, DataTable *data, double wl)
{
	int     i, j, k;
	double  dl;
	int n;
	double_complex sigma;
	double_complex px, py;

	for (i = 0; i < port->np; i++) {
		for (j = 0; j < port->np; j++) {
			port->NN[i][j] = 0.0;
			port->NNy[i][j] = 0.0;
		}
	}  

	n = -1;
	for (k = 0; k < port->np-2; k += 2) {
		if (n == -1 || k == port->nn[n]) { // in each layer
			n++;
			if (data->par.modeID == 2) {
				px = 1.0/port->epT[n].xx;
				py = px;
			} 
			else {
				px = py = 1.0;
			}
		}

		// the length of line element
		dl = fabs(port->yp[k+2]-port->yp[k]);

		port->NN[k][k] += dl*2.0/15.0*px;
		port->NN[k][k+1] += dl/15.0*px;
		port->NN[k][k+2] -= dl/30.0*px;
		port->NN[k+1][k] += dl/15.0*px;
		port->NN[k+1][k+1] += dl*8.0/15.0*px;
		port->NN[k+1][k+2] += dl/15.0*px;
		port->NN[k+2][k] -= dl/30.0*px;
		port->NN[k+2][k+1] += dl/15.0*px;
		port->NN[k+2][k+2] += dl*2.0/15.0*px;

		port->NNy[k][k] += -3.0/6.0*py;
		port->NNy[k][k+1] += -1.0/6.0*py;
		port->NNy[k][k+2] += 4.0/6.0*py;
		port->NNy[k+1][k] += 1.0/6.0*py;
		port->NNy[k+1][k+1] += 3.0/6.0*py;
		port->NNy[k+1][k+2] += -4.0/6.0*py;
		port->NNy[k+2][k] += -4.0/6.0*py;
		port->NNy[k+2][k+1] += 4.0/6.0*py;
		port->NNy[k+2][k+2] += 0.0*py;
	}
}

extern int Fvector(DataTable *data, double wl, int ll, int ii)
{
	int i, j, m, n;
	PortData *port = data->port.data[ii];
	double_complex f1, ft;
	double_complex cj(0.0, 1.0);
	double_complex sigma, px[MAXMODE], py[MAXMODE];

	for (i = 0; i < port->np; i++) {
		f1 = 0.0; n = 0;
		for (j = 0; j < port->np; j++) {
			if (port->beta2[ll] > 0.0) {
				// *added take care for number of port is not 1
				f1 += cj*sqrt(port->beta2[ll])*port->Amp
					*(port->NN[i][j]*port->PsiR[ll][j]);	
			}
		}    
		data->fem.b[port->kp[i]] += f1;
	}
	if (port->kp[0] > data->fem.nr) data->fem.b[port->kp[0]] = 0.0;
	if (port->kp[port->np-1] > data->fem.nr) data->fem.b[port->kp[port->np-1]] = 0.0;
}

extern int FEManalysis(DataTable *data, double wl, int ll)
{
	int     i, j, k, l, ip, ii, jj, it, m, n;
	int     N1, N2, N3, N4, N5, N6, iel=0, count;
	double  Ae, A1, A2, A3, B1, B2, B3;
	double  C1, C2, C3, LL1, LL2, LL3;
	XYCZdouble Node1, Node2, Node3, Node4, Node5, Node6;
	int     bpmNE, bpmNP;
	int     element[2000][3];
	FEM     *fem = &(data->fem);
	PortData *port;
	double  eps;
	double  le, a1, a2, b1, b2, L1, L2;
	double  le30, yy1, yy2, yy3, y31;
	double_complex cj(0.0, 1.0);
	CGtable *cg_table = data->fem.cg_table;

	for (i = 0; i < data->fem.np; i++) {
		data->fem.phi[i] = 0.0;
	}

	for (k = 0; k < data->port.number; k++) {
		if (data->par.inflag[k] == 1) {
			data->port.data[k]->Amp = (data->port.forward+data->port.backward);
		}
		else {
			data->port.data[k]->Amp = 0.0;
		}
	}

	// left hand side matrix
	Amatrix(data, cg_table, wl);

	// right hand side vector
	if (data->port.solver == POINTSOURCE) {
		for (i = 0; i < data->fem.np; i++) {
			if (i == data->port.centerkp) {
			fprintf(stderr, "centerkp = %d\n", data->port.centerkp);
			data->fem.b[i] = -cj*data->par.k0*Z0;
			}
			else {
				data->fem.b[i] = 0.0;
			}
		}

	} else if (data->port.solver == SHEETSOURCE) {
		for (i = 0; i < data->fem.np; i++) {
			if (data->port.xflag == 1) {
				if (fabs(data->port.x0-data->fem.node[i].x) < 1.0e-9) {
					data->fem.b[i] = -cj*data->par.k0*Z0;
				} else {
					data->fem.b[i] = 0.0;
				}
			}
			else {
				if (fabs(data->port.y0-data->fem.node[i].y) < 1.0e-9) {
					data->fem.b[i] = -cj*data->par.k0*Z0;
				}
				else {
					data->fem.b[i] = 0.0;
				}
			}
		}
	}
	else {
		for (i = 0; i < data->fem.np; i++) data->fem.b[i] = 0.0;

		for (k = 0; k < data->port.number; k++) {
			if (data->par.inflag[k] == 1) {
				Bmatrix(data->port.data[k], data, wl);            
				Fvector(data, wl, ll, k);    
			}
		}
	}


	// LU factorization and backward substitution
	callPARDISO(data, data->fem.A, data->fem.b, 
			data->fem.nr, data->fem.phi, data->fem.cg_table, 0);

	fprintf(stderr, "free -------------\n");
	callPARDISO(data, data->fem.A, data->fem.b, 
		data->fem.nr, data->fem.phi, data->fem.cg_table, -1);

	fprintf(stderr, "free -------------\n");
}

extern int rotation(Tensor *ep, Tensor *epR, double cosA, double sinA)
{
	Tensor rot1;

	rot1.xx = ep->xx*cosA-ep->xy*sinA;
	rot1.yx = ep->yx*cosA-ep->yy*sinA;
	rot1.zx = ep->zx*cosA-ep->zy*sinA;
	rot1.xy = ep->xx*sinA+ep->xy*cosA;
	rot1.yy = ep->yx*sinA+ep->yy*cosA;
	rot1.zy = ep->zx*sinA+ep->zy*cosA;
	rot1.xz = ep->xz;
	rot1.yz = ep->yz;
	rot1.zz = ep->zz;

	epR->xx = rot1.xx*cosA-rot1.yx*sinA;
	epR->xy = rot1.xy*cosA-rot1.yy*sinA;
	epR->xz = rot1.xz*cosA-rot1.yz*sinA;
	epR->yx = rot1.xx*sinA+rot1.yx*cosA;
	epR->yy = rot1.xy*sinA+rot1.yy*cosA;
	epR->yz = rot1.xz*sinA+rot1.yz*cosA;
	epR->zx = rot1.zx;
	epR->zy = rot1.zy;
	epR->zz = rot1.zz;
}
