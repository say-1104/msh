#include <stdio.h>
#include <math.h>

#define MAX(a, b) ((a >  b) ? a : b)
#define MIN(a, b) ((a <  b) ? a : b)

/* Table of constant values */

	static int c__9 = 9;
	static int c__1 = 1;
	static int c__3 = 3;

static int geigb_(int *, int *, int *, double *, double *, double *,
		double *, double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, int *, int *, int *, 
	    double *, double *, double *, double *, int *, int *, int *);
static int bisecg_(int *, int *, int *, double *, double *, double *,
		int *, double *, double *, double *, double *, double *, 
	    double *, int *, int *, double *, double *, 
	    int *, double *, double *);
static int mortho_(int *, int *, int *, double *, double *, double *,
		double *, double *, double *);
static int gsturm_(int *, int *, double *, double *, double *, double *, int *);
static int gsylv2_(double *, double *, int *, int *, double *, double *,
		int *, double *, double *, double *, int *);
static int gsturm_(int *, int *, double *, double *, double *, double *, int *);
static int sort_(int *, double *, double *, double *, int *, int *, int *,
		double *, double *, double *);
static int bglu2_(double *, int *, int *, int *, double *, double *,
		double *, int *, int *);
static int bslu2_(double *, int *, int *, double *, double *, double *, int *);
static int brylg_(int *, int *, double *, double *, double *, double *,
		double *, double *, int *, int *, double *, int *, 
	    int *, double *, double *, int *, int *, double *);
static int bgslv4_(double *, int *, int *, int *, double *, int *);
static int subs3g_(int *, int *, double *, double *, double *, double *,
		int *, int *, int *, double *, int *);
static int bsslv4_(double *, int *, int *, double *);
static int rylsfg_(int *, int *, double *, double *, double *, double *,
		double *, double *, double *, double *, int *, double *,
	     double *, int *, int *, int *, int *, 
	    double *, double *, double *, double *);
static int yax_(int *, int *, double *, double *, double *);
static int bglu2_(double *, int *, int *, int *, double *, double *,
		double *, int *, int *);
static int brylg_(int *, int *, double *, double *, double *, double *,
		double *, double *, int *, int *, double *, int *, 
	    int *, double *, double *, int *, int *, double *);
static int bgslv4_(double *, int *, int *, int *, double *, int *);
static int bgsv4c_(double *, int *, int *, int *, int *, double *, int *);
static int bssv4c_(double *, int *, int *, int *, double *);
static int jacobi_(double *, int *, double *, double *, double *, double *);

extern int eigv4p_(double *xa, double *xm, int *n, int *mu,
	int *ne0, int *ne, int *ialg, double *eval, double *eig, double *eigv,
	double *rn, double *wk, int *iwk, int *ier)
{
    /* System generated locals */
    int xa_dim1, xa_offset, xm_dim1, xm_offset, eigv_dim1, eigv_offset, 
	    wk_dim1, wk_offset, iwk_dim1, iwk_offset, i__1, i__2, i__3;
    double d__1;

    /* Local variables */
    static double gapt, epse, epsg;
    static int lmax;
    static double epsx, vtvl, xsup, epse0, epse1, epsg0;
    static int i__, j, k;
    static double s, epsbi;
    static int llmin;
    static double epsbi0, epsbi1, epsgl0;
    static int nn;
    static double epsglu;
    static double epssyl;
    static double xr1, xr2;
    static int nbi, nct;
    static double syl0;
    static int nxr1;

    /* Parameter adjustments */
    iwk_dim1 = *n + 1;
    iwk_offset = iwk_dim1 + 1;
    iwk -= iwk_offset;
    wk_dim1 = *n + 1;
    wk_offset = wk_dim1 + 1;
    wk -= wk_offset;
    eigv_dim1 = *n;
    eigv_offset = eigv_dim1 + 1;
    eigv -= eigv_offset;
    --eig;
    xm_dim1 = *mu + 1;
    xm_offset = xm_dim1;
    xm -= xm_offset;
    xa_dim1 = *mu + 1;
    xa_offset = xa_dim1;
    xa -= xa_offset;
    --rn;

    /* Function Body */
    epsx = 4.0;
    epse0 = 3.52e-15;
    epsbi0 = 2.2e-10;
    syl0 = epsbi0 * (double)4.0;
    epsg0 = 2.2e-10;
    epsgl0 = 2.2e-17;
    vtvl = 1.0e-4;
    nn = 4;
    nbi = 6;
    lmax = MAX(25,*mu);
    llmin = 6;

	/* computation of matrix norm  */
    xsup = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
		s = 0.;
		/* Computing MAX */
		i__2 = j - *mu;
		i__3 = j - 1;
		for (i__ = MAX(i__2,1); i__ <= i__3; ++i__) {
	    	s += (d__1 = xa[j - i__ + i__ * xa_dim1], fabs(d__1));
		}
		/* Computing MIN */
		i__2 = j + *mu;
		i__3 = MIN(i__2,*n);
		for (i__ = j; i__ <= i__3; ++i__) {
	    	s += (d__1 = xa[i__ - j + j * xa_dim1], fabs(d__1));
		}
		xsup = MAX(xsup,s);
    }

    epse = (double) (*mu + 1) * 4.0 * epse0 * xsup;
    epse1 = (double) (*mu + 1) * epse0;
    epsbi = xsup / ((double) (*n) * epsx);
    epsbi1 = epsbi0 * xsup;
    epssyl = syl0 * xsup;
    epsg = epsg0 * xsup;
    epsglu = epsgl0 * xsup;

	xr1 = 0.0;
    gsturm_(n, mu, &xa[xa_offset], &xm[xm_offset], &wk[wk_dim1 * 36 + 1],
			&xr1, &nct);
fprintf(stderr, "nct(%lf) = %d\n", xr1, nct);
    nxr1 = nct;
	xr2 = xsup;
    gsturm_(n, mu, &xa[xa_offset], &xm[xm_offset], &wk[wk_dim1 * 36 + 1],
			&xr2, &nct);
fprintf(stderr, "nct(%lf) = %d\n", xr2, nct);

    if (*ialg >= 10) {
		gsturm_(n, mu, &xa[xa_offset], &xm[xm_offset], &wk[wk_dim1 * 36 + 1], 
				eval, &nct);
		*ne = MIN(nct, *ne0);
    }

	/* Bisection */
    *ne = MIN(*ne, *n);
    bisecg_(n, mu, ne, &xa[xa_offset], &xm[xm_offset], &wk[wk_dim1 * 36 + 1], 
	    &iwk[iwk_dim1 + 1], &wk[wk_dim1 + 1], &wk[(wk_dim1 << 1) + 1], 
		&xsup, &epsbi, &epsbi1, &epssyl, &nn, &nbi, &wk[wk_dim1 * 3 + 1], 
		&wk[(wk_dim1 << 2) + 1], &nxr1, &xr1, &xr2);

	/* Subspace Iteration */
    iwk[(iwk_dim1 << 1) + 1] = 1;
    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
		if (iwk[k + iwk_dim1] == iwk[k - 1 + iwk_dim1]) {
	    	iwk[k + (iwk_dim1 << 1)] = iwk[k - 1 + (iwk_dim1 << 1)] + 1;
		} else {
	    	iwk[k + (iwk_dim1 << 1)] = 1;
		}
    }

    geigb_(n, mu, ne, &xa[xa_offset], &xm[xm_offset], &eig[1],
			&eigv[eigv_offset], &rn[1], &wk[wk_dim1 * 37 + 1],
			&wk[wk_dim1 * 6 + 1], &wk[wk_dim1 + 1], &wk[(wk_dim1 << 1) + 1],
			&epsbi1, &epse, &epse1, &epsg, &epsglu, &nn, &lmax, &llmin,
			&wk[wk_dim1 * 3 + 1], &wk[(wk_dim1 << 2) + 1],
			&wk[wk_dim1 * 5 + 1], &wk[wk_dim1 * 6 + 1], &iwk[iwk_dim1 + 1],
			&iwk[(iwk_dim1 << 1) + 1], &iwk[iwk_dim1 * 3 + 1]);

	/* Error Analysis */
    gapt = xsup * 1.0e-4;
    mortho_(n, mu, ne, &eigv[eigv_offset], &eig[1], &xm[xm_offset], &vtvl,
			&gapt, &wk[wk_dim1 * 3 + 1]);
    return 0;
}

static int bisecg_(int *n, int *m1, int *ne, double *xa, double *xm,
	double *b, int *nc, double *xi, double *xs, double *xsup,
	double *epsbi, double *epsbi1, double *epssyl, int *nn, int *nbie,
	double *wk, double *wk1, int *nxr1, double *xr1, double *xr2)
{
    /* System generated locals */
    int xa_dim1, xa_offset, xm_dim1, xm_offset, b_dim1, b_offset, i__1;
    double d__1, d__2;

    /* Local variables */
    static double alfa;
    static int nsyl, i__, k, kflag;
    static double w;
    static int index;
    static double epsbiw;
    static int nbi, nct;

    /* Parameter adjustments */
    --wk1;
    --wk;
    --xi;
    --nc;
    b_dim1 = (*m1 << 1) + 1 - (-(*m1) - 1) + 1;
    b_offset = -(*m1) - 1 + b_dim1;
    b -= b_offset;
    xm_dim1 = *m1 + 1;
    xm_offset = xm_dim1;
    xm -= xm_offset;
    xa_dim1 = *m1 + 1;
    xa_offset = xa_dim1;
    xa -= xa_offset;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		xi[i__] = *xr1;
		xs[i__] = *xsup;
		nc[i__] = *n - *nxr1;
    }
    xs[0] = *xr1;
    k = 1;

    alfa = *xr2;
    gsylv2_(&xa[xa_offset], &xm[xm_offset], n, m1, &alfa, epssyl, &nct,
			&b[b_offset], &wk[1], &wk1[1], &index);
    if (index != 0) {
		gsturm_(n, m1, &xa[xa_offset], &xm[xm_offset], &b[b_offset], &alfa,
				&nct);
    }
    nct -= *nxr1;
    i__1 = nct;
    for (i__ = k; i__ <= i__1; ++i__) {
		xs[i__] = alfa;
		nc[i__] = nct;
    }
    i__1 = *n;
    for (i__ = nct + 1; i__ <= i__1; ++i__) {
		/* Computing MAX */
		d__1 = alfa, d__2 = xi[i__];
		xi[i__] = MAX(d__1, d__2);
    }

	/* Bisection */
L1:
    kflag = 0;
    if (k <= 2) {
		epsbiw = *epsbi1;
    } else {
		epsbiw = *epsbi;
    }

L2:
    if (nc[k] > k || k != 1 && xi[k] == xs[k - 1]) {
		w = xs[k] - xi[k];
		if (w >= epsbiw || w >= *epsbi1 && nc[k] - k >= *nn) {
		    alfa = (xi[k] + xs[k]) / 2.;
		    if (kflag == 0) {
				nsyl = 0;
L3:
				gsylv2_(&xa[xa_offset], &xm[xm_offset], n, m1, &alfa, epssyl, 
					&nct, &b[b_offset], &wk[1], &wk1[1], &index);
fprintf(stderr, "%lf -- %d (index = %d)\n", alfa, nct, index);
				/* Retry to Count */
				if (index != 0) {
				    ++nsyl;
			    	if (nsyl == 1) {
						alfa = xi[k] * (double)0.75 + xs[k] * (double)0.25;
			    	}
			    	if (nsyl == 2) {
						alfa = xi[k] * (double)0.25 + xs[k] * (double)0.75;
			    	}
					/* write(*,*) ' *isyl=1*',k,nsyl,nc(k),xi (k),alfa,xs(k) */
			    	if (nsyl <= 2) {
						goto L3;
			    	}
				}
		    }
	
		    if (nsyl > 2 || kflag >= 1) {
				++kflag;
				alfa = (xi[k] + xs[k]) / 2.;
				gsturm_(n, m1, &xa[xa_offset], &xm[xm_offset], &b[b_offset],
						&alfa, &nct);
		    }
		    nct -= *nxr1;
	
		    i__1 = nct;
		    for (i__ = k; i__ <= i__1; ++i__) {
				xs[i__] = alfa;
				nc[i__] = nct;
		    }
		    i__1 = *n;
		    for (i__ = nct + 1; i__ <= i__1; ++i__) {
				/* Computing MAX */
				d__1 = alfa, d__2 = xi[i__];
				xi[i__] = MAX(d__1,d__2);
			}
		    goto L2;
		}
    } else { /* Isolated Eigenvalues */
		i__1 = *nbie;
		for (nbi = 1; nbi <= i__1; ++nbi) {
	    	if (xs[k] - xi[k] >= epsbiw) {
				alfa = (xi[k] + xs[k]) / 2.;
				gsylv2_(&xa[xa_offset], &xm[xm_offset], n, m1, &alfa, epssyl, 
					&nct, &b[b_offset], &wk[1], &wk1[1], &index);
				if (index != 0) {
				    gsturm_(n, m1, &xa[xa_offset], &xm[xm_offset], &b[b_offset],
					&alfa, &nct);
				}
				nct -= *nxr1;
				if (nct == k) {
				    xs[k] = alfa;
				} else {
				    xi[k] = alfa;
				}
	    	} else {
				goto L4;
	    	}
		}
L4:
	;
    }
    k = nc[k] + 1;
    if (k <= *ne) {
		goto L1;
    }
    return 0;
}

static int gsylv2_(double *xa, double *xm, int *n, int *mu, double *alfa,
	double *eps, int *nct, double *a, double *wk1, double *wk2, int *ier)
{
    /* System generated locals */
    int xa_dim1, xa_offset, xm_dim1, xm_offset, a_dim1, a_offset, i__1, 
	    i__2, i__3, i__4;
    double d__1;

    /* Local variables */
    static int i__, j, k;
    static double t;
    static int k1;
    static double t1;

    /* Parameter adjustments */
    a_dim1 = *mu + 2 + 1;
    a_offset = a_dim1;
    a -= a_offset;
    xm_dim1 = *mu + 1;
    xm_offset = xm_dim1;
    xm -= xm_offset;
    xa_dim1 = *mu + 1;
    xa_offset = xa_dim1;
    xa -= xa_offset;
    --wk1;
    --wk2;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* Computing MIN */
		i__3 = i__ + *mu;
		i__2 = MIN(i__3,*n);
		for (j = i__; j <= i__2; ++j) {
	    	a[j - i__ + i__ * a_dim1] = xa[j - i__ + i__ * xa_dim1]
				- *alfa *xm[j - i__ + i__ * xm_dim1];
		}
		/* Computing MIN */
		i__2 = i__ + *mu;
		i__3 = i__ + *mu + 2;
		for (j = MIN(i__2, *n) + 1; j <= i__3; ++j) {
		/* L22: */
	   		a[j - i__ + i__ * a_dim1] = 0.;
		}
    }
    if (*n % 2 != 0) {
		a[(*n + 1) * a_dim1] = 1.;
		i__1 = *mu + 2;
		for (j = 1; j <= i__1; ++j) {
	    	a[j + (*n + 1) * a_dim1] = 0.;
		}
    }

    *nct = 0;
    *ier = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; k += 2) {
		k1 = k + 1;

		if (a[k * a_dim1] < (double)0.) {
		    ++(*nct);
		}

		if ((d__1 = a[k * a_dim1], fabs(d__1)) <= *eps) {
		    *ier = k;
		    return 0;
		}

		t = -a[k * a_dim1 + 1] / a[k * a_dim1];
		/* Computing MIN */
		i__2 = k1 + *mu;
		i__3 = MIN(i__2, *n);
		for (j = k + 1; j <= i__3; ++j) {
		    wk1[j] = a[j - k + k * a_dim1];
		    a[j - k1 + k1 * a_dim1] += t * wk1[j];
		    wk2[j] = a[j - k1 + k1 * a_dim1];
		}

		if (a[k1 * a_dim1] < (double)0.) {
		    ++(*nct);
		}

		if ((d__1 = a[k1 * a_dim1], fabs(d__1)) <= *eps) {
		    *ier = k1;
	   		return 0;
		}

		/* Gaussian Elimination for K-th and K+1'th Process. */
		/* Computing MIN */
		i__2 = k1 + *mu;
		i__3 = MIN(i__2, *n);
		for (i__ = k1 + 1; i__ <= i__3; ++i__) {
		    t = -a[i__ - k + k * a_dim1] / a[k * a_dim1];
		    t1 = -a[i__ - k1 + k1 * a_dim1] / a[k1 * a_dim1];
			/* Computing MIN */
	   		i__4 = k1 + *mu;
	    	i__2 = MIN(i__4,*n);
	    	for (j = i__; j <= i__2; ++j) {
				a[j - i__ + i__ * a_dim1] = a[j - i__ + i__ * a_dim1]
						+ t * wk1[j] + t1 * wk2[j];
	 	  	}
		}
    }
    return 0;
}

static int gsturm_(int *n, int *m1, double *xa, double *xm,
		double *b, double *alfa, int *nct)
{
    /* System generated locals */
    int xa_dim1, xa_offset, xm_dim1, xm_offset, b_dim1, b_offset, i__1, 
	    i__2, i__3, i__4;
    double d__1, d__2;

    /* Local variables */
    static int jmin, jmax, i__, j, k;
    static double t, w;
    static int sgndt, psign, sgndt0, psign0, ip;

    /* Parameter adjustments */
    b_dim1 = (*m1 << 1) - (-(*m1)) + 1;
    b_offset = -(*m1) + b_dim1;
    b -= b_offset;
    xm_dim1 = *m1 + 1;
    xm_offset = xm_dim1;
    xm -= xm_offset;
    xa_dim1 = *m1 + 1;
    xa_offset = xa_dim1;
    xa -= xa_offset;

    /* Function Body */
    sgndt0 = 1;
    ip = 1;
    psign0 = 1;
    *nct = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* Computing MAX */
		i__2 = i__ - *m1;
		i__3 = i__ - 1;
		for (j = MAX(i__2,1); j <= i__3; ++j) {
	    	b[j - i__ + i__ * b_dim1] = xa[i__ - j + j * xa_dim1]
					- *alfa * xm[i__ - j + j * xm_dim1];
		}
		/* Computing MIN */
		i__2 = i__ + *m1;
		i__3 = MIN(i__2,*n);
		for (j = i__; j <= i__3; ++j) {
	    	b[j - i__ + i__ * b_dim1] = xa[j - i__ + i__ * xa_dim1]
				- *alfa * xm[j - i__ + i__ * xm_dim1];
		}
		/* Computing MIN */
		i__3 = i__ + *m1;
		/* Computing MIN */
		i__4 = i__ + (*m1 << 1);
		i__2 = MIN(i__4,*n);
		for (j = MIN(i__3,*n) + 1; j <= i__2; ++j) {
	    	b[j - i__ + i__ * b_dim1] = 0.;
		}
    }

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
		/* Martin's Special GAussian Elimination */
		/* Computing MAX */
		i__2 = 1, i__3 = k - *m1;
		jmin = MAX(i__2, i__3);
		/* Computing MIN */
		i__2 = k + *m1;
		jmax = MIN(i__2, *n);
		i__2 = k - 1;
		for (i__ = jmin; i__ <= i__2; ++i__) {
	    	if ((d__1 = b[i__ * b_dim1], fabs(d__1)) <
				(d__2 = b[i__ - k + k * b_dim1], fabs(d__2))) {
				ip = -ip;
				i__3 = jmax;
				for (j = i__; j <= i__3; ++j) {
		    		w = b[j - i__ + i__ * b_dim1];
		    		b[j - i__ + i__ * b_dim1] = b[j - k + k * b_dim1];
		    		b[j - k + k * b_dim1] = w;
				}
	    	}
	    	if ((d__1 = b[i__ * b_dim1], fabs(d__1)) > (double)0.0) {
				t = -b[i__ - k + k * b_dim1] / b[i__ * b_dim1];
				/* Voption Indep(B) */
				i__3 = jmax;
				for (j = i__; j <= i__3; ++j) {
		    		b[j - k + k * b_dim1] += t * b[j - i__ + i__ * b_dim1];
				}
	    	}
		}

		/* Count for Sign Change */
		if (k - *m1 <= 1) {
	    	psign = 1;
	    	i__2 = k;
	    	for (i__ = 1; i__ <= i__2; ++i__) {
				if (b[i__ * b_dim1] < (double)0.) {
		    		psign = -psign;
				}
	    	}
	    	sgndt = ip * psign;
		} else {
	    	if (b[(k - *m1 - 1) * b_dim1] < (double)0.) {
				psign0 = -psign0;
	    	}
	    	psign = 1;
	    	i__2 = k;
		    for (i__ = k - *m1; i__ <= i__2; ++i__) {
				if (b[i__ * b_dim1] < (float)0.) {
		    		psign = -psign;
				}
	    	}
	    	sgndt = ip * psign * psign0;
		}
		if (sgndt0 * sgndt < 0) {
	    	++(*nct);
		}
		sgndt0 = sgndt;
    }
    return 0;
}

static int geigb_(int *n, int *m1, int *ne, double *xa, double *xm,
	double *eig, double *v, double *r__, double *b, double *u,
	double *xi, double *xs, 
	double *epsbi1, double *epse, double *epse1, double *epsg,
	double *epsglu, int *nn, int *lmax, int *llmin, 
	double *w, double *aw, double *hv, double *wk1, 
	int *nc, int *ncli, int *ip)
{
    /* System generated locals */
    int xa_dim1, xa_offset, xm_dim1, xm_offset, v_dim1, v_offset, b_dim1, 
	    b_offset, u_dim1, u_offset, i__1, i__2, i__3;
    double d__1, d__2;

    /* Local variables */
    static double alfa;
    static int mbak, nbak;
    static double eigl;
    static int neig;
    static double eigr, eigx, gapx;
    static int lmin;
    static double rmax;
    static int isym;
    static int i__, j, k, l;
    static double s;
    static int irand;
    static double eigrm;
    static int kkmin, kxend, kkmax;
    static double ratio;
    static int lconv, k1;
    static double rnorm;
    static int jj, kk, ir, ll, kx, ix;
    static double ratmax;
    static int ibr, lle, llm, lls, num;
    static double eig0, gap0;


/*        GEIGB */
/*               COPYRIGHT : K.MURATA, H.HASEGAWA, OCT.  4 1991 V.1 */

/*               SUBSPACE ITERATION FOR SYMMETRIC BAND MATRIX */
/*              FROM THE SMALLEST EIGENVALUE TO THE LARGEST EIGENVALUE. */

    /* Parameter adjustments */
    --ip;
    --ncli;
    --nc;
    --wk1;
    --hv;
    --aw;
    --w;
    --xi;
    u_dim1 = *n + 1;
    u_offset = u_dim1 + 1;
    u -= u_offset;
    --r__;
    v_dim1 = *n;
    v_offset = v_dim1 + 1;
    v -= v_offset;
    --eig;
    b_dim1 = (*m1 << 1) + 1 - (-(*m1) - 1) + 1;
    b_offset = -(*m1) - 1 + b_dim1;
    b -= b_offset;
    xm_dim1 = *m1 + 1;
    xm_offset = xm_dim1;
    xm -= xm_offset;
    xa_dim1 = *m1 + 1;
    xa_offset = xa_dim1;
    xa -= xa_offset;

    /* Function Body */
    irand = 12345;
    ratmax = (double)16.0;
    k = 1;
	/* Computes the Range of Reorthogonalization */
L1:
    alfa = (xi[k] + xs[k]) / 2.;
	/* Computing MAX */
    i__1 = k - 1;
    mbak = MAX(i__1,1);
	/* Computing MAX */
    i__1 = mbak - ncli[mbak] + 1;
    nbak = MAX(i__1,1);
    gap0 = xs[k] - alfa;
L2:
    gapx = alfa - xs[nbak];
    ratio = gapx / gap0;
    if (ratio < ratmax && nbak > 1) {
		/* Computing MAX */
		i__1 = nbak - ncli[nbak];
		mbak = MAX(i__1,1);
		/* Computing MAX */
		i__1 = mbak - ncli[mbak] + 1;
		nbak = MAX(i__1,1);
		goto L2;
    }

    if (xs[k] == xi[nc[k] + 1] && nc[k] < *ne) {
		kkmax = nc[nc[k] + 1] - k;
    } else {
		kkmax = nc[k] - k;
    }
	/* Computing MAX */
    i__1 = 1, i__2 = k - 1;
    k1 = MAX(i__1, i__2);
    eigl = xi[k] - (xi[k] - xi[nc[k1]]) / 64.0;
    eigr = xs[k];
    eigrm = eigr - (xs[k] - xi[k]) / 16.0;
    if (k > 1 && xs[k - 1] == xi[k]) {
		kkmin = nc[k] - k + 1;
    } else {
		kkmin = 0;
    }
	/* M-Orthogonalize */
    i__1 = k + kkmin - 1;
    for (j = k; j <= i__1; ++j) {
		i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    w[i__] = v[i__ + j * v_dim1];
	}
	yax_(n, m1, &xm[xm_offset], &w[1], &hv[1]);
	i__2 = k - 1;
	for (jj = nbak; jj <= i__2; ++jj) {
	    s = 0.;
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
			s += hv[i__] * v[i__ + jj * v_dim1];
	    }
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
			w[i__] -= s * v[i__ + jj * v_dim1];
	    }
	}
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    v[i__ + j * v_dim1] = w[i__];
	}
    }

    i__1 = kkmax;
    for (kk = kkmin; kk <= i__1; ++kk) {
		/* Initial Vector Using Uniform Random Numbers */
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    irand = (irand * 1229 + 351750) % 1664501;
		    w[i__] = irand * 6.0078065438230435e-7;
		}
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    w[i__] = w[i__] * 2.0 - 1.0;
		}
		s = 0.;
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    s += w[i__] * w[i__];
		}
		s = 1.0 / sqrt(s);
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    w[i__] *= s;
		}
		yax_(n, m1, &xm[xm_offset], &w[1], &hv[1]);
		i__2 = k + kkmin - 1;
		for (j = nbak; j <= i__2; ++j) {
		    s = 0.;
		    i__3 = *n;
		    for (i__ = 1; i__ <= i__3; ++i__) {
				s += hv[i__] * v[i__ + j * v_dim1];
		    }
		    i__3 = *n;
		    for (i__ = 1; i__ <= i__3; ++i__) {
				w[i__] -= s * v[i__ + j * v_dim1];
		    }
		}
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    v[i__ + (k + kk) * v_dim1] = w[i__];
		}
    }

	/* LU Decomposition for A-ALFA*M */
    if (kkmax >= *nn << 1) {
		alfa = (xi[k] * (double)5.0 + xs[k] * (double)3.0) / (double)8.0;
    }
    isym = 0;
L3:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
		/* Computing MAX */
		i__3 = i__ - *m1;
		i__2 = MAX(i__3,1) - 1;
		for (j = i__ - *m1 - 1; j <= i__2; ++j) {
	    	b[j - i__ + i__ * b_dim1] = 0.;
		}
		/* Computing MAX */
		i__2 = i__ - *m1;
		i__3 = i__ - 1;
		for (j = MAX(i__2,1); j <= i__3; ++j) {
	    	b[j - i__ + i__ * b_dim1] = xa[i__ - j + j * xa_dim1]
				- alfa * xm[i__ - j + j * xm_dim1];
		}
		/* Computing MIN */
		i__2 = i__ + *m1;
		i__3 = MIN(i__2,*n);
		for (j = i__; j <= i__3; ++j) {
		    b[j - i__ + i__ * b_dim1] = xa[j - i__ + i__ * xa_dim1]
				- alfa * xm[j - i__ + i__ * xm_dim1];
		}
		/* Computing MIN */
		i__3 = i__ + *m1;
		i__2 = i__ + (*m1 << 1) + 1;
		for (j = MIN(i__3,*n) + 1; j <= i__2; ++j) {
		    b[j - i__ + i__ * b_dim1] = 0.;
		}
    }

    if (isym == 0) {
		bslu2_(&b[b_offset], n, m1, epsg, &aw[1], &wk1[1], &ir);
		if (ir != 0) {
		    isym = 1;
		    goto L3;
		}
    } else {
		bglu2_(&b[b_offset], n, m1, m1, epsglu, &wk1[1], &aw[1], &ip[1], &ir);
    }

	/* Subspace Iterations */
    ibr = 0;
    if (kkmax > 0 && kkmax <= 30) {
	ibr = 1;
	kxend = nc[k];
	num = kkmax + 1;
	lconv = nc[k] - k + 1;
	llm = *llmin - 1;

	i__1 = llm;
	for (ll = 1; ll <= i__1; ++ll) {
	    subs3g_(n, m1, &xm[xm_offset], &b[b_offset], &v[v_offset], &hv[1],
		     &ip[1], &k, &num, &u[u_offset], &isym);
/* L1010: */
	}

	lls = *llmin;
L4:
	i__1 = k + kkmax;
	for (j = k; j <= i__1; ++j) {
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* L1030: */
		w[i__] = v[i__ + j * v_dim1];
	    }
	    yax_(n, m1, &xm[xm_offset], &w[1], &hv[1]);
	    i__2 = k - 1;
	    for (jj = nbak; jj <= i__2; ++jj) {
		s = 0.;
		i__3 = *n;
		for (i__ = 1; i__ <= i__3; ++i__) {
/* L1050: */
		    s += hv[i__] * v[i__ + jj * v_dim1];
		}
		i__3 = *n;
		for (i__ = 1; i__ <= i__3; ++i__) {
/* L1060: */
		    w[i__] -= s * v[i__ + jj * v_dim1];
		}
/* L1040: */
	    }
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* L1070: */
		v[i__ + j * v_dim1] = w[i__];
	    }
/* L1020: */
	}

/* Computing MIN */
	i__1 = lls + 5;
	lle = MIN(i__1,*lmax);
	i__1 = lle;
	for (ll = lls; ll <= i__1; ++ll) {
	    subs3g_(n, m1, &xm[xm_offset], &b[b_offset], &v[v_offset], &hv[1],
		     &ip[1], &k, &num, &u[u_offset], &isym);
	    neig = 0;
	    rmax = 0.;
	    i__2 = k + num - 1;
	    for (kx = k; kx <= i__2; ++kx) {
		yax_(n, m1, &xa[xa_offset], &v[kx * v_dim1 + 1], &aw[1]);
		yax_(n, m1, &xm[xm_offset], &v[kx * v_dim1 + 1], &hv[1]);
		eigx = 0.;
		i__3 = *n;
		for (i__ = 1; i__ <= i__3; ++i__) {
/* L1130: */
		    eigx += v[i__ + kx * v_dim1] * aw[i__];
		}
		eig[kx] = eigx;
		s = 0.;
		i__3 = *n;
		for (ix = 1; ix <= i__3; ++ix) {
		    aw[ix] -= eigx * hv[ix];
/* L1140: */
/* Computing 2nd power */
		    d__1 = aw[ix];
		    s += d__1 * d__1;
		}
		r__[kx] = sqrt(s);
		if (eigx <= eigr + *epse1 && eigx >= eigl && r__[kx] < *epse *
			 .25) {
		    ++neig;
		}
/* Computing MAX */
		d__1 = rmax, d__2 = r__[kx];
		rmax = MAX(d__1,d__2);
/* L1110: */
	    }
	    if (neig >= lconv) {
		goto L5;
	    }
/* L1100: */
	}
	lls = lle + 1;
	if (lle < *lmax) {
	    goto L4;
	}

L5:
	sort_(n, &v[v_offset], &eig[1], &r__[1], &k, &num, &lconv, &eigl, &
		eigr, epse);
	i__1 = lconv - 1;
	for (kk = 0; kk <= i__1; ++kk) {
	    if (r__[k + kk] <= *epse && eig[k + kk] <= eigr + *epse1) {
		goto L1200;
	    }
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* L1210: */
		w[i__] = v[i__ + (k + kk) * v_dim1];
	    }
	    yax_(n, m1, &xm[xm_offset], &w[1], &hv[1]);
	    i__2 = k + kk - 1;
	    for (j = k; j <= i__2; ++j) {
		s = 0.;
		i__3 = *n;
		for (i__ = 1; i__ <= i__3; ++i__) {
/* L1230: */
		    s += hv[i__] * v[i__ + j * v_dim1];
		}
		i__3 = *n;
		for (i__ = 1; i__ <= i__3; ++i__) {
/* L1240: */
		    w[i__] -= s * v[i__ + j * v_dim1];
		}
/* L1220: */
	    }

	    yax_(n, m1, &xm[xm_offset], &w[1], &hv[1]);
	    i__2 = k + lconv - 1;
	    for (j = k + kk + 1; j <= i__2; ++j) {
		if (r__[j] <= *epse) {
		    s = 0.;
		    i__3 = *n;
		    for (i__ = 1; i__ <= i__3; ++i__) {
/* L1260: */
			s += hv[i__] * v[i__ + j * v_dim1];
		    }
		    i__3 = *n;
		    for (i__ = 1; i__ <= i__3; ++i__) {
/* L1270: */
			w[i__] -= s * v[i__ + j * v_dim1];
		    }
		}
/* L1250: */
	    }

/* Computing MIN */
	    d__1 = eig[k + kk];
	    eig[k + kk] = MIN(d__1,eigrm);
	    yax_(n, m1, &xm[xm_offset], &w[1], &hv[1]);
	    i__2 = k + kk;
	    rylsfg_(n, m1, &xa[xa_offset], &xm[xm_offset], &b[b_offset], &eig[
		    1], &v[v_offset], &r__[1], &w[1], &aw[1], &ip[1], epsglu, 
		    epse, &nbak, &i__2, &kxend, &ibr, &eigr, &hv[1], epse1, &
		    wk1[1]);
L1200:
	    ;
	}
/*             INVERSE ITERATIONS */
    } else {
	kxend = k;
	lmin = 4;
	if (xs[k] - xi[k] < *epsbi1 * (float)3.) {
	    lmin = 8;
	}

	yax_(n, m1, &xm[xm_offset], &w[1], &hv[1]);
	i__1 = nc[k] - k;
	for (kk = 0; kk <= i__1; ++kk) {
	    i__2 = *lmax;
	    for (l = 1; l <= i__2; ++l) {
		i__3 = *n;
		for (i__ = 1; i__ <= i__3; ++i__) {
/* L2020: */
		    w[i__] = hv[i__];
		}
		if (isym == 0) {
		    bsslv4_(&b[b_offset], n, m1, &w[1]);
		} else {
		    bgslv4_(&b[b_offset], n, m1, m1, &w[1], &ip[1]);
		}

		i__3 = k + kk;
		brylg_(n, m1, &xa[xa_offset], &xm[xm_offset], &v[v_offset], &
			w[1], &eig0, &rnorm, &nbak, &i__3, &aw[1], &kxend, &
			ibr, epse, &r__[1], &l, &lmin, &hv[1]);
		if (rnorm < *epse && l >= lmin && eig0 <= eigr + *epse1) {
		    eig[k + kk] = MAX(0.,eig0);
		    r__[k + kk] = rnorm;
		    i__3 = *n;
		    for (i__ = 1; i__ <= i__3; ++i__) {
/* L2030: */
			v[i__ + (k + kk) * v_dim1] = w[i__];
		    }
		    goto L6;
		}
/* L2010: */
	    }

	    r__[k + kk] = rnorm;
	    eig[k + kk] = MIN(eig0,eigrm);
	    i__2 = k + kk;
	    rylsfg_(n, m1, &xa[xa_offset], &xm[xm_offset], &b[b_offset], &eig[
		    1], &v[v_offset], &r__[1], &w[1], &aw[1], &ip[1], epsglu, 
		    epse, &nbak, &i__2, &kxend, &ibr, &eigr, &hv[1], epse1, &
		    wk1[1]);
L6:
/* L2000: */
	    ;
	}
    }

    k = nc[k] + 1;
    if (k <= *ne) {
	goto L1;
    }
    return 0;
} /* geigb_ */


static int rylsfg_(int *n, int *m1, double *xa, 
	double *xm, double *b, double *eig, double *v, 
	double *r__, double *w, double *aw, int *ip, 
	double *epsglu, double *epse, int *nbak, int *kx, 
	int *kxend, int *ibr, double *eigr, double *hv, 
	double *epse1, double *wk1)
{
    /* Initialized data */

    static int llmin = 2;

    /* System generated locals */
    int xa_dim1, xa_offset, xm_dim1, xm_offset, v_dim1, v_offset, b_dim1, 
	    b_offset, i__1, i__2, i__3;

    /* Local variables */
    static int i__, j;
    static double eigre;
    static double rnorm;
    static int ir, ll, lll;
    static double eig0;


/*        RYLSFG */
/*               COPYRIGHT : K.MURATA, H.HASEGAWA, OCT.  4 1991 V.1 */

/*               INVERSE ITERATION WITH RAYLEIGH SHIFT. */

    /* Parameter adjustments */
    --wk1;
    --hv;
    --ip;
    --aw;
    --w;
    --r__;
    v_dim1 = *n;
    v_offset = v_dim1 + 1;
    v -= v_offset;
    --eig;
    b_dim1 = (*m1 << 1) + 1 - (-(*m1) - 1) + 1;
    b_offset = -(*m1) - 1 + b_dim1;
    b -= b_offset;
    xm_dim1 = *m1 + 1;
    xm_offset = xm_dim1;
    xm -= xm_offset;
    xa_dim1 = *m1 + 1;
    xa_offset = xa_dim1;
    xa -= xa_offset;

    /* Function Body */

    eigre = *eigr + *epse1;
    for (lll = 1; lll <= 2; ++lll) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    i__3 = i__ - *m1;
	    i__2 = MAX(i__3,1) - 1;
	    for (j = i__ - *m1 - 1; j <= i__2; ++j) {
/* L120: */
		b[j - i__ + i__ * b_dim1] = 0.;
	    }
/* Computing MAX */
	    i__2 = i__ - *m1;
	    i__3 = i__ - 1;
	    for (j = MAX(i__2,1); j <= i__3; ++j) {
/* L130: */
		b[j - i__ + i__ * b_dim1] = xa[i__ - j + j * xa_dim1] - eig[*
			kx] * xm[i__ - j + j * xm_dim1];
	    }
/* Computing MIN */
	    i__2 = i__ + *m1;
	    i__3 = MIN(i__2,*n);
	    for (j = i__; j <= i__3; ++j) {
/* L140: */
		b[j - i__ + i__ * b_dim1] = xa[j - i__ + i__ * xa_dim1] - eig[
			*kx] * xm[j - i__ + i__ * xm_dim1];
	    }
/* Computing MIN */
	    i__3 = i__ + *m1;
	    i__2 = i__ + (*m1 << 1) + 1;
	    for (j = MIN(i__3,*n) + 1; j <= i__2; ++j) {
/* L150: */
		b[j - i__ + i__ * b_dim1] = 0.;
	    }
/* L110: */
	}
/*             LU DECOMPOSITION FOR A-EIG(KX)*I */
	bglu2_(&b[b_offset], n, m1, m1, epsglu, &wk1[1], &aw[1], &ip[1], &ir);
/*             INVERSE ITERATIONS( RETRY ) */
	for (ll = 1; ll <= 10; ++ll) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L210: */
		w[i__] = hv[i__];
	    }
	    bgslv4_(&b[b_offset], n, m1, m1, &w[1], &ip[1]);
	    brylg_(n, m1, &xa[xa_offset], &xm[xm_offset], &v[v_offset], &w[1],
		     &eig0, &rnorm, nbak, kx, &aw[1], kxend, ibr, epse, &r__[
		    1], &ll, &llmin, &hv[1]);
	    if (rnorm < *epse && ll >= llmin && eig0 < eigre) {
/* *          WRITE(*,600)  ' K  :', KX, ' LLL :', LLL, ' LL:'
,LL ,RNORM */
		goto L1;
	    }
/* L200: */
	}
/*             CONVERGENCE CHECK */
	if (rnorm < *epse * 2.) {
	    goto L1;
	}
	eig[*kx] = eig0;
/* L100: */
    }
/* L600: */
/*             STORE */
L1:
    eig[*kx] = MAX(0.,eig0);
    r__[*kx] = rnorm;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L300: */
	v[i__ + *kx * v_dim1] = w[i__];
    }
    return 0;
} /* rylsfg_ */


static int brylg_(int *n, int *m1, double *xa, 
	double *xm, double *v, double *w, double *eig, 
	double *rnorm, int *nbak, int *kx, double *aw, 
	int *kxend, int *ibr, double *epse, double *r__, 
	int *l, int *lmin, double *hv)
{
    /* System generated locals */
    int xa_dim1, xa_offset, xm_dim1, xm_offset, v_dim1, v_offset, i__1, 
	    i__2;
    double d__1, d__2, d__3;

    /* Local variables */
    static double wmax;
    static int j;
    static double s;
    static int ix;


/*        BRYLG */
/*               COPYRIGHT : K.MURATA, H.HASEGAWA, OCT.  4 1991 V.1 */

/*               COMPUTES RAYLEIGH QUOTIENT. */


    /* Parameter adjustments */
    --hv;
    --r__;
    --aw;
    --w;
    v_dim1 = *n;
    v_offset = v_dim1 + 1;
    v -= v_offset;
    xm_dim1 = *m1 + 1;
    xm_offset = xm_dim1;
    xm -= xm_offset;
    xa_dim1 = *m1 + 1;
    xa_offset = xa_dim1;
    xa -= xa_offset;

    /* Function Body */
    wmax = fabs(w[1]);
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
/* L125: */
/* Computing MAX */
	d__2 = wmax, d__3 = (d__1 = w[j], fabs(d__1));
	wmax = MAX(d__2,d__3);
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* L135: */
	w[j] /= wmax;
    }
/*             M-ORTHOGONALIZE */
    yax_(n, m1, &xm[xm_offset], &w[1], &hv[1]);
    i__1 = *kx - 1;
    for (j = *nbak; j <= i__1; ++j) {
	s = 0.;
	i__2 = *n;
	for (ix = 1; ix <= i__2; ++ix) {
/* L110: */
	    s += hv[ix] * v[ix + j * v_dim1];
	}
	i__2 = *n;
	for (ix = 1; ix <= i__2; ++ix) {
/* L120: */
	    w[ix] -= s * v[ix + j * v_dim1];
	}
/* L100: */
    }

    if (*ibr == 1) {
	i__1 = *kxend;
	for (j = *kx + 1; j <= i__1; ++j) {
	    if (r__[j] < *epse) {
		s = 0.;
		i__2 = *n;
		for (ix = 1; ix <= i__2; ++ix) {
/* L210: */
		    s += hv[ix] * v[ix + j * v_dim1];
		}
		i__2 = *n;
		for (ix = 1; ix <= i__2; ++ix) {
/* L220: */
		    w[ix] -= s * v[ix + j * v_dim1];
		}
	    }
/* L200: */
	}
    }
/*             M-NORMALIZE */
    yax_(n, m1, &xm[xm_offset], &w[1], &hv[1]);
    s = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* L130: */
	s += w[j] * hv[j];
    }
    s = sqrt(s);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* L140: */
	w[j] /= s;
    }
/*             COMPUTES AN EIGENVALUE */
    if (*l >= *lmin) {
	yax_(n, m1, &xa[xa_offset], &w[1], &aw[1]);
	yax_(n, m1, &xm[xm_offset], &w[1], &hv[1]);
	*eig = 0.;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* L170: */
	    *eig += w[j] * aw[j];
	}
/*             COMPUTES RESIDUAL NORM */
	s = 0.;
	i__1 = *n;
	for (ix = 1; ix <= i__1; ++ix) {
	    aw[ix] -= *eig * hv[ix];
/* L180: */
/* Computing 2nd power */
	    d__1 = aw[ix];
	    s += d__1 * d__1;
	}
	*rnorm = sqrt(s);
    }
    return 0;
} /* brylg_ */


static int sort_(int *n, double *v, double *eig, 
	double *r__, int *k, int *num, int *lconv, double 
	*eigl, double *eigr, double *epse)
{
    /* System generated locals */
    int v_dim1, v_offset, i__1, i__2, i__3;

    /* Local variables */
    static double work;
    static int i__;
    static double am;
    static int kx, kkk, ipk;


/*        SORT */
/*               COPYRIGHT : K.MURATA, H.HASEGAWA, OCT.  4 1991 V.1 */

/*               SORT EIGENVALUES, EIGENVECTORS AND RESIDUAL NORMS */
/*               WITH ASCENDING ORDER. */


    /* Parameter adjustments */
    --r__;
    --eig;
    v_dim1 = *n;
    v_offset = v_dim1 + 1;
    v -= v_offset;

    /* Function Body */
/* Computing MIN */
    i__2 = *lconv, i__3 = *num - 1;
    i__1 = *k + MIN(i__2,i__3);
    for (kx = *k; kx <= i__1; ++kx) {
	ipk = kx;
	am = *eigr;
	i__2 = *k + *num - 1;
	for (kkk = kx; kkk <= i__2; ++kkk) {
	    if (eig[kkk] < am && eig[kkk] >= *eigl && r__[kkk] < *epse * 16.) 
		    {
		am = eig[kkk];
		ipk = kkk;
	    }
/* L210: */
	}
	if (ipk != kx) {
	    work = eig[kx];
	    eig[kx] = eig[ipk];
	    eig[ipk] = work;
	    work = r__[kx];
	    r__[kx] = r__[ipk];
	    r__[ipk] = work;
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		work = v[i__ + kx * v_dim1];
		v[i__ + kx * v_dim1] = v[i__ + ipk * v_dim1];
/* L220: */
		v[i__ + ipk * v_dim1] = work;
	    }
	}
/* L200: */
    }
    return 0;
} /* sort_ */


static int subs3g_(int *n, int *m1, double *xm, 
	double *b, double *v, double *hv, int *ip, int *k,
	 int *num, double *u, int *isym)
{
    /* System generated locals */
    int xm_dim1, xm_offset, v_dim1, v_offset, b_dim1, b_offset, u_dim1, 
	    u_offset, i__1, i__2, i__3;

    /* Local variables */
    static double d__[50], h__[2500]	/* was [50][50] */;
    static int i__, j;
    static double s, w[2500]	/* was [50][50] */, b1[50];
    static int ix, kx;
    static double dinvrt;
    static int kk1;


/*        SUBS3G */
/*               COPYRIGHT : K.MURATA, H.HASEGAWA, OCT.  4 1991 V.1 */

/*               SUBSPACE ITERATION FOR EACH SUBSPACE. */


    /* Parameter adjustments */
    u_dim1 = *n + 1;
    u_offset = u_dim1 + 1;
    u -= u_offset;
    --ip;
    --hv;
    v_dim1 = *n;
    v_offset = v_dim1 + 1;
    v -= v_offset;
    b_dim1 = (*m1 << 1) + 1 - (-(*m1) - 1) + 1;
    b_offset = -(*m1) - 1 + b_dim1;
    b -= b_offset;
    xm_dim1 = *m1 + 1;
    xm_offset = xm_dim1;
    xm -= xm_offset;

    /* Function Body */
    i__1 = *k + *num - 1;
    for (kx = *k; kx <= i__1; ++kx) {
	kk1 = kx - *k + 1;
	yax_(n, m1, &xm[xm_offset], &v[kx * v_dim1 + 1], &hv[1]);
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
/* L220: */
	    u[j + kk1 * u_dim1] = hv[j];
	}
/* L200: */
    }

    if (*isym == 0) {
	bssv4c_(&b[b_offset], n, m1, num, &u[u_offset]);
    } else {
	bgsv4c_(&b[b_offset], n, m1, m1, num, &u[u_offset], &ip[1]);
    }

    i__1 = *num;
    for (j = 1; j <= i__1; ++j) {
	yax_(n, m1, &xm[xm_offset], &u[j * u_dim1 + 1], &hv[1]);
	i__2 = *num;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    s = 0.;
	    i__3 = *n;
	    for (ix = 1; ix <= i__3; ++ix) {
/* L330: */
		s += u[ix + i__ * u_dim1] * hv[ix];
	    }
	    h__[i__ + j * 50 - 51] = s;
/* L320: */
	}
/* L300: */
    }

    jacobi_(h__, num, d__, w, b1, &hv[1]);

    i__1 = *k + *num - 1;
    for (kx = *k; kx <= i__1; ++kx) {
	kk1 = kx - *k + 1;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L410: */
	    v[i__ + kx * v_dim1] = 0.;
	}
	i__2 = *num;
	for (ix = 1; ix <= i__2; ++ix) {
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/* L430: */
		v[i__ + kx * v_dim1] += u[i__ + ix * u_dim1] * w[ix + kk1 * 
			50 - 51];
	    }
/* L420: */
	}
	dinvrt = 1. / sqrt(d__[kk1 - 1]);
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L440: */
	    v[i__ + kx * v_dim1] *= dinvrt;
	}
/* L400: */
    }
    return 0;
} /* subs3g_ */


static int jacobi_(double *a, int *n, double *eig, 
	double *v, double *b, double *z__)
{
    /* System generated locals */
    int i__1, i__2, i__3;
    double d__1, d__2, d__3, d__4;

    /* Local variables */
    static int krot;
    static double g, h__;
    static int i__, j, k, p, q;
    static double tresh;
    static int kount2;
    static double th, sm, tu;
    static int kcount;
    static double tan__, cos__, cot, sin__;


/*        JACOBI */
/*               COPYRIGHT : K.MURATA, Y.OGAWA, OCT.  4 1991 V.1 */

/*               COMPUTES EIGENVALUES AND EIGENVECTORS BY RUTISHAUSER'S */
/*               JACOBI METHOD FOR GENERAL SYMMETRIC MATRICES. */
/*               TRANSLATED FROM CONTRIBUTION 2/1 OF LINEAR ALGEBRA */
/*               ( J.H.WILKINSON AND C.REINSCH ), SPRINGER, 1971. */

/*        INPUT - - */
/*             A(N,N)   R *8  : 2-DIM. ARRAY CONTAINING THE GENERAL */
/*                              SYMMETRIC MATRIX. */
/*             N        I *4  : ORDER OF MATRIX. */
/*        OUTPUT - - */
/*             EIG(N)   R *8  : 1-DIM. ARRAY. CONTAINING THE EIGENVALUE. 
*/
/*             V(N,N)   R *8  : 2-DIM. ARRAY. CONTAINING THE EIGENVECTOR. 
*/
/*        WORKING  - */
/*             B(N), Z(N) */
/*                      R *8  : 1-DIM. ARRAY. */

/* *    DIMENSION A(N,*), V(N,*), EIG(*), B(*), Z(*) */

    /* Parameter adjustments */
    --z__;
    --b;
    v -= 51;
    --eig;
    a -= 51;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L10: */
	    v[i__ + j * 50] = 0.;
	}
	v[j + j * 50] = 1.;
/* L20: */
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	eig[i__] = a[i__ + i__ * 50];
	b[i__] = eig[i__];
	z__[i__] = 0.;
/* L30: */
    }
    krot = 0;
    kcount = 0;
    kount2 = 0;

/* L1000: */
    for (k = 1; k <= 10; ++k) {
	sm = 0.;
	i__1 = *n - 1;
	for (p = 1; p <= i__1; ++p) {
	    i__2 = *n;
	    for (q = p + 1; q <= i__2; ++q) {
/* L35: */
		sm += (d__1 = a[p + q * 50], fabs(d__1));
	    }
/* L40: */
	}
	if (sm < 1e-72) {
	    ++kount2;
	    goto L2000;
	}

	if (k < 4) {
/* Computing MAX */
/* Computing 2nd power */
	    i__1 = *n;
	    d__1 = sm * .2 / (double) (i__1 * i__1);
	    tresh = MAX(d__1,1e-36);
	} else {
	    tresh = 1e-36;
	}

	i__1 = *n;
	for (p = 1; p <= i__1; ++p) {
	    i__2 = *n;
	    for (q = p + 1; q <= i__2; ++q) {
		g = (d__1 = a[p + q * 50], fabs(d__1)) * 100.;
		if (k > 4 && (d__1 = eig[p], fabs(d__1)) + g == (d__2 = eig[p],
			 fabs(d__2)) && (d__3 = eig[q], fabs(d__3)) + g == (
			d__4 = eig[q], fabs(d__4))) {
		    a[p + q * 50] = 0.;
		} else {
		    if ((d__1 = a[p + q * 50], fabs(d__1)) > tresh) {
			h__ = eig[q] - eig[p];
			if (fabs(h__) + g == fabs(h__)) {
			    tan__ = a[p + q * 50] / h__;
			} else {
			    th = h__ * .5 / a[p + q * 50];
/* Computing 2nd power */
			    d__1 = th;
			    cot = fabs(th) + sqrt(d__1 * d__1 + 1);
			    tan__ = 1. / cot;
			    if (th < 0.) {
				tan__ = -tan__;
			    }
			}
/* Computing 2nd power */
			d__1 = tan__;
			cos__ = 1. / sqrt(d__1 * d__1 + 1.);
			sin__ = tan__ * cos__;
			tu = sin__ / (cos__ + 1.);
			h__ = tan__ * a[p + q * 50];
			z__[p] -= h__;
			z__[q] += h__;
			eig[p] -= h__;
			eig[q] += h__;
			a[p + q * 50] = 0.;
			i__3 = p - 1;
			for (j = 1; j <= i__3; ++j) {
			    g = a[j + p * 50];
			    h__ = a[j + q * 50];
			    a[j + p * 50] = g - sin__ * (h__ + g * tu);
			    a[j + q * 50] = h__ + sin__ * (g - h__ * tu);
/* L50: */
			}
			i__3 = q - 1;
			for (j = p + 1; j <= i__3; ++j) {
			    g = a[p + j * 50];
			    h__ = a[j + q * 50];
			    a[p + j * 50] = g - sin__ * (h__ + g * tu);
			    a[j + q * 50] = h__ + sin__ * (g - h__ * tu);
/* L60: */
			}
			i__3 = *n;
			for (j = q + 1; j <= i__3; ++j) {
			    g = a[p + j * 50];
			    h__ = a[q + j * 50];
			    a[p + j * 50] = g - sin__ * (h__ + g * tu);
			    a[q + j * 50] = h__ + sin__ * (g - h__ * tu);
/* L70: */
			}
			i__3 = *n;
			for (j = 1; j <= i__3; ++j) {
			    g = v[j + p * 50];
			    h__ = v[j + q * 50];
			    v[j + p * 50] = g - sin__ * (h__ + g * tu);
			    v[j + q * 50] = h__ + sin__ * (g - h__ * tu);
/* L80: */
			}
			++krot;
		    }
		}
/* L100: */
	    }
/* L200: */
	}
	i__1 = *n;
	for (p = 1; p <= i__1; ++p) {
	    b[p] += z__[p];
	    eig[p] = b[p];
	    z__[p] = 0.;
/* L210: */
	}
/* L300: */
    }
    ++kcount;
L2000:
    return 0;
} /* jacobi_ */


static int bslu2_(double *a, int *n, int *mu, 
	double *eps, double *wk1, double *wk2, int *ier)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2, i__3, i__4;
    double d__1;

    /* Local variables */
    static int i__, j, k;
    static double t;
    static int k1;
    static double t1;


/*        BSLU2 */
/*               COPYRIGHT : H.HASEGAWA, OCT.  4 1991 V.1 */

/*               SOLVES SIMULTANEOUS LINEAR SYMMETRIC BAND EQUATIONS */
/*               BY SYMMETRIC GAUSSIAN ELIMINATION METHOD. */

/*        INPUT - - */
/*             A(0:MU+2,N+1) */
/*                      R *8  : 2-DIM. ARRAY FOR UPPER TRIANGULAR MATRIX. 
*/
/*             N        I *4  : ORDER OF MATRIX. */
/*             MU       I *4  : UPPER BAND WIDTH. */
/*             EPS      R *8  : PARAMETER TO CHECK SINGULARITY OF THE */
/*                              MATRIX. ( STANDARD VALUE 3.52D-15 ) */
/*        OUTPUT - - */
/*             A(0:MU+2,N+1)  : RESULT OF GAUSSIAN ELIMINATION. */
/*             IER      I *4  : = 0,  FOR NORMAL EXECUTION. */
/*                              = 1,  FOR SINGULAR MATRIX. */
/*                              = 3,  FOR INVALID ARGUEMENT. */
/*        WORKING  - */
/*             WK1(N), WK2(N) */
/*                      R *8  : 1-DIM. ARRAY. */

/* *    DIMENSION A(0:MU+2,*), WK1(*), WK2(*) */
/*             LEFT HAND SIDE */
    /* Parameter adjustments */
    a_dim1 = (*mu << 1) + 1 - (-(*mu) - 1) + 1;
    a_offset = -(*mu) - 1 + a_dim1;
    a -= a_offset;
    --wk1;
    --wk2;

    /* Function Body */
    if (*eps < 0.) {
	*eps = 3.52e-15;
    }
    if (*n <= 0 || *mu <= 0 || *mu >= *n) {
	*ier = 3;
	fprintf(stderr, " --- error ---\n");
/*
	s_wsle(&io___144);
	do_lio(&c__9, &c__1, "  (SUBR. BSLU2)  INVALID ARGUMENT.  MU, N =", 
		43L);
	do_lio(&c__3, &c__1, (char *)&(*mu), (ftnlen)sizeof(int));
	do_lio(&c__3, &c__1, (char *)&(*n), (ftnlen)sizeof(int));
	e_wsle();
*/
	return 0;
    }

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a[*mu + 1 + i__ * a_dim1] = 0.;
/* L10: */
	a[*mu + 2 + i__ * a_dim1] = 0.;
    }
    if (*n % 2 != 0) {
	a[(*n + 1) * a_dim1] = 1.;
	i__1 = *mu + 2;
	for (j = 1; j <= i__1; ++j) {
/* L20: */
	    a[j + (*n + 1) * a_dim1] = 0.;
	}
    }

    *ier = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; k += 2) {
	k1 = k + 1;

	if ((d__1 = a[k * a_dim1], fabs(d__1)) <= *eps) {
	    *ier = 1;
/* *         WRITE(*,*)  '  (SUBR. BSLU2)  MATRIX IS SINGULAR AT K
 =', K */
	    return 0;
	}

	t = -a[k * a_dim1 + 1] / a[k * a_dim1];
/* Computing MIN */
	i__3 = k1 + *mu;
	i__2 = MIN(i__3,*n);
	for (j = k + 1; j <= i__2; ++j) {
	    wk1[j] = a[j - k + k * a_dim1];
	    a[j - k1 + k1 * a_dim1] += t * wk1[j];
	    wk2[j] = a[j - k1 + k1 * a_dim1];
/* L110: */
	}

	if ((d__1 = a[k1 * a_dim1], fabs(d__1)) <= *eps) {
	    *ier = 1;
/* *         WRITE(*,*)  '  (SUBR. BSLU2)  MATRIX IS SINGULAR AT K
 =', K1 */
	    return 0;
	}
/*             GAUSSIAN ELIMINATION FOR K-TH AND K+1'TH PROCESS. */
/* Computing MIN */
	i__3 = k1 + *mu;
	i__2 = MIN(i__3,*n);
	for (i__ = k1 + 1; i__ <= i__2; ++i__) {
	    t = -a[i__ - k + k * a_dim1] / a[k * a_dim1];
	    t1 = -a[i__ - k1 + k1 * a_dim1] / a[k1 * a_dim1];
/* Computing MIN */
	    i__4 = k1 + *mu;
	    i__3 = MIN(i__4,*n);
	    for (j = i__; j <= i__3; ++j) {
/* L140: */
		a[j - i__ + i__ * a_dim1] = a[j - i__ + i__ * a_dim1] + t * 
			wk1[j] + t1 * wk2[j];
	    }
/* L130: */
	}
/* L100: */
    }
    return 0;
} /* bslu2_ */


static int bssv4c_(double *a, int *n, int *mu, int * m, double *b)
{
    /* System generated locals */
    int a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static int nend, i__, j, k, l;
    static double t;
    static int k1;
    static double s0, s1, bk, bk1;


/*        BSSV4C */
/*               COPYRIGHT : K.MURATA, H.HASEGAWA, OCT.  4 1991 V.1 */

/*               SOLVES SIMULTANEOUS LINEAR EQUATIONS WITH MULTIPLE */
/*               RIGHT HAND SIDE BY SYMMETRIC GAUSSIAN ELIMINATION. */

/*        INPUT - - */
/*             A(0:MU+2,N+1) */
/*                      R *8  : RESULT OF GAUSSIAN ELIMINATION. */
/*             N        I *4  : ORDER OF MATRIX. */
/*             MU       I *4  : UPPER BAND WIDTH. */
/*             M        I *4  : NUMBER OF RIGHT HAND SIDE. */
/*             B(N+1,M) R *8  : 2-DIM. ARRAY CONTAINING THE RIGHT HAND */
/*                              SIDE VECTORS. */
/*        OUTPUT - - */
/*             B(N+1,M)       : SOLUTION. */

/* *    DIMENSION A(0:MU+2,*), B(N+1,*) */
/*             FORWARD ELIMINATION PROCESS */
    /* Parameter adjustments */
    b_dim1 = *n + 1;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    a_dim1 = (*mu << 1) + 1 - (-(*mu) - 1) + 1;
    a_offset = -(*mu) - 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; k += 2) {
	k1 = k + 1;
	t = -a[k * a_dim1 + 1] / a[k * a_dim1];
	i__2 = *m;
	for (l = 1; l <= i__2; ++l) {
	    b[k1 + l * b_dim1] += t * b[k + l * b_dim1];
	    bk = -b[k + l * b_dim1] / a[k * a_dim1];
	    bk1 = -b[k1 + l * b_dim1] / a[k1 * a_dim1];
/* Computing MIN */
	    i__4 = k1 + *mu;
	    i__3 = MIN(i__4,*n);
	    for (i__ = k1 + 1; i__ <= i__3; ++i__) {
/* L110: */
		b[i__ + l * b_dim1] = b[i__ + l * b_dim1] + a[i__ - k + k * 
			a_dim1] * bk + a[i__ - k1 + k1 * a_dim1] * bk1;
	    }
/* L105: */
	}
/* L100: */
    }
/*             BACKWARD SUBSTITUTION PROCESS */
    if (*n % 2 == 0) {
	nend = *n;
    } else {
	nend = *n + 1;
	i__1 = *m;
	for (l = 1; l <= i__1; ++l) {
/* L201: */
	    b[*n + 1 + l * b_dim1] = 0.;
	}
    }
    for (k = nend; k >= 1; k += -2) {
	i__1 = *m;
	for (l = 1; l <= i__1; ++l) {
	    s1 = -b[k + l * b_dim1];
	    s0 = -b[k - 1 + l * b_dim1];
/* Computing MIN */
	    i__3 = k + *mu;
	    i__2 = MIN(i__3,*n);
	    for (j = k + 1; j <= i__2; ++j) {
		s1 += a[j - k + k * a_dim1] * b[j + l * b_dim1];
/* L210: */
		s0 += a[j - k + 1 + (k - 1) * a_dim1] * b[j + l * b_dim1];
	    }
	    b[k + l * b_dim1] = -s1 / a[k * a_dim1];
	    b[k - 1 + l * b_dim1] = (-s0 - a[(k - 1) * a_dim1 + 1] * b[k + l *
		     b_dim1]) / a[(k - 1) * a_dim1];
/* L300: */
	}
/* L200: */
    }
    return 0;
} /* bssv4c_ */


static int bsslv4_(double *a, int *n, int *mu, double *b)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static int nend, i__, j, k;
    static double t;
    static int k1;
    static double s0, s1, bk, bk1;


/*        BSSLV4 */
/*               COPYRIGHT : H.HASEGAWA, OCT.  4 1991 V.1 */

/*               SOLVES SIMULTANEOUS LINEAR SYMMETRIC BAND EQUATIONS */
/*               BY SYMMETRIC GAUSSIAN ELIMINATION METHOD. */

/*        INPUT - - */
/*             A(0:MU+2,N+1) */
/*                      R *8  : RESULT OF GAUSSIAN ELIMINATION. */
/*             N        I *4  : ORDER OF MATRIX. */
/*             MU       I *4  : UPPER BAND WIDTH. */
/*             B(N+1)   R *8  : 1-DIM. ARRAY CONTAINING THE RIGHT HAND */
/*                              SIDE VECTOR. */
/*        OUTPUT - - */
/*             B(N+1)         : SOLUTION. */

/* *    DIMENSION A(0:MU+2,*), B(*) */
/*             FORWARD ELIMINATION PROCESS */
    /* Parameter adjustments */
    a_dim1 = (*mu << 1) + 1 - (-(*mu) - 1) + 1;
    a_offset = -(*mu) - 1 + a_dim1;
    a -= a_offset;
    --b;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; k += 2) {
	k1 = k + 1;
	t = -a[k * a_dim1 + 1] / a[k * a_dim1];
	b[k1] += t * b[k];
	bk = -b[k] / a[k * a_dim1];
	bk1 = -b[k1] / a[k1 * a_dim1];
/* Computing MIN */
	i__3 = k1 + *mu;
	i__2 = MIN(i__3,*n);
	for (i__ = k1 + 1; i__ <= i__2; ++i__) {
/* L110: */
	    b[i__] = b[i__] + a[i__ - k + k * a_dim1] * bk + a[i__ - k1 + k1 *
		     a_dim1] * bk1;
	}
/* L100: */
    }
/*             BACKWARD SUBSTITUTION PROCESS */
    if (*n % 2 == 0) {
	nend = *n;
    } else {
	nend = *n + 1;
	b[*n + 1] = 0.;
    }
    for (k = nend; k >= 1; k += -2) {
	s1 = -b[k];
	s0 = -b[k - 1];
/* Computing MIN */
	i__2 = k + *mu;
	i__1 = MIN(i__2,*n);
	for (j = k + 1; j <= i__1; ++j) {
	    s1 += a[j - k + k * a_dim1] * b[j];
/* L210: */
	    s0 += a[j - k + 1 + (k - 1) * a_dim1] * b[j];
	}
	b[k] = -s1 / a[k * a_dim1];
	b[k - 1] = (-s0 - a[(k - 1) * a_dim1 + 1] * b[k]) / a[(k - 1) * 
		a_dim1];
/* L200: */
    }
    return 0;
} /* bsslv4_ */


static int bglu2_(double *a, int *n, int *ml, int *
	mu, double *eps, double *wk1, double *wk2, int *ip, 
	int *ier)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2, i__3;
    double d__1;

    /* Builtin functions */
/*
    int s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen), 
	    e_wsle();
*/

    /* Local variables */
    static double amax;
    static int imax, jmax;
    static double amax1;
    static int imax1, jmax1, i__, j, k;
    static double t, w;
    static int k1;
    static double t1, aik;
    static int ipk, ipk1;


/*        BGLU2 */
/*               COPYRIGHT : H.HASEGAWA, OCT.  4 1991 V.1 */

/*               SOLVES SIMULTANEOUS LINEAR EQUATIONS */
/*               BY GAUSSIAN ELIMINATION METHOD FOR GENERAL BAND MATRIX. 
*/

/*        INPUT - - */
/*             A(-ML-1:MU+ML+1,N+1) */
/*                      R *8  : 2-DIM. ARRAY CONTAINING REAL BAND MATRIX. 
*/
/*             N        I *4  : ORDER OF MATRIX. */
/*             ML       I *4  : LOWER BAND WIDTH. */
/*             MU       I *4  : UPPER BAND WIDTH. */
/*             EPS      R *8  : PARAMETER TO CHECK SINGULARITY OF THE */
/*                              MATRIX. ( STANDARD VALUE 3.52D-15 ) */
/*        OUTPUT - - */
/*             A(-ML-1:MU+ML+1,N+1) */
/*                            : RESULT OF GAUSSIAN ELIMINATION. */
/*             IP(N+1)  I *4  : PIVOT NUMBER. */
/*             IER      I *4  : = 0,  FOR NORMAL EXECUTION. */
/*                              = 1,  FOR SINGULAR MATRIX. */
/*                              = 3,  FOR INVALID ARGUEMENT. */
/*        WORKING  - */
/*             WK1(N), WK2(N) */
/*                      R *8  : 1-DIM. ARRAY. */

/*             LEFT HAND SIDE */
    /* Parameter adjustments */
    a_dim1 = *mu + *ml + 1 - (-(*ml) - 1) + 1;
    a_offset = -(*ml) - 1 + a_dim1;
    a -= a_offset;
    --wk1;
    --wk2;
    --ip;

    /* Function Body */
    if (*eps < 0.) {
	*eps = 3.52e-15;
    }
    if (*n <= 0 || *ml <= 0 || *mu <= 0 || *ml >= *n || *mu >= *n) {
	*ier = 3;
	fprintf(stderr, " --- error ---\n");
/*
	s_wsle(&io___172);
	do_lio(&c__9, &c__1, "  (SUBR. BGLU2)  INVALID ARGUMENT.  ML, MU, N ="
		, 47L);
	do_lio(&c__3, &c__1, (char *)&(*ml), (ftnlen)sizeof(int));
	do_lio(&c__3, &c__1, (char *)&(*mu), (ftnlen)sizeof(int));
	do_lio(&c__3, &c__1, (char *)&(*n), (ftnlen)sizeof(int));
	e_wsle();
*/
	return 0;
    }

    *ier = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a[-(*ml) - 1 + i__ * a_dim1] = 0.;
/* L10: */
	a[*mu + *ml + 1 + i__ * a_dim1] = 0.;
    }
    if (*n % 2 != 0) {
	i__1 = *mu + *ml + 1;
	for (j = -(*ml) - 1; j <= i__1; ++j) {
/* L20: */
	    a[j + (*n + 1) * a_dim1] = 0.;
	}
	a[(*n + 1) * a_dim1] = 1.;
    }

    i__1 = *n;
    for (k = 1; k <= i__1; k += 2) {

	k1 = k + 1;
/* Computing MIN */
	i__2 = k + *ml;
	imax = MIN(i__2,*n);
/* Computing MIN */
	i__2 = k + *mu + *ml;
	jmax = MIN(i__2,*n);
/* Computing MIN */
	i__2 = k1 + *ml;
	imax1 = MIN(i__2,*n);
/* Computing MIN */
	i__2 = k1 + *mu + *ml;
	jmax1 = MIN(i__2,*n);
/*             FIND MAXIMUM ELEMENT IN THE K-TH COLUMN. */
	amax = (d__1 = a[k * a_dim1], fabs(d__1));
	ipk = k;
	i__2 = imax;
	for (i__ = k + 1; i__ <= i__2; ++i__) {
	    aik = (d__1 = a[k - i__ + i__ * a_dim1], fabs(d__1));
	    if (aik > amax) {
		ipk = i__;
		amax = aik;
	    }
/* L110: */
	}
	ip[k] = ipk;
/*             EXCHANGE FOR K-TH ELIMINATION. */
	if (amax > *eps) {
	    if (ipk != k) {
		i__2 = jmax;
		for (j = k; j <= i__2; ++j) {
		    w = a[j - ipk + ipk * a_dim1];
		    a[j - ipk + ipk * a_dim1] = a[j - k + k * a_dim1];
		    a[j - k + k * a_dim1] = w;
/* L120: */
		}
	    }
/*             COMPUTE ALFA */
	    t = a[k * a_dim1];
	    i__2 = imax1;
	    for (i__ = k + 1; i__ <= i__2; ++i__) {
/* L130: */
		a[k - i__ + i__ * a_dim1] = -a[k - i__ + i__ * a_dim1] / t;
	    }
	    t = a[k * a_dim1 + 1];
	    i__2 = imax;
	    for (i__ = k + 1; i__ <= i__2; ++i__) {
/* L140: */
		a[k1 - i__ + i__ * a_dim1] += a[k - i__ + i__ * a_dim1] * t;
	    }
/*             MATRIX IS SINGULAR. */
	} else {
	    *ier = 1;
	    ip[k] = k;
	    i__2 = imax1;
	    for (i__ = k + 1; i__ <= i__2; ++i__) {
/* L150: */
		a[k - i__ + i__ * a_dim1] = 0.;
	    }
	    a[k * a_dim1] = *eps * .5;
/* *         WRITE(*,*)  '  (SUBR. BGLU2)  MATRIX IS SINGULAR AT K
 =', K */
/* *         RETURN */
	}
/*             FIND MAXIMUM ELEMENT IN THE K+1'TH COLUMN. */
	amax1 = (d__1 = a[k1 * a_dim1], fabs(d__1));
	ipk1 = k1;
	i__2 = imax1;
	for (i__ = k1 + 1; i__ <= i__2; ++i__) {
	    aik = (d__1 = a[k1 - i__ + i__ * a_dim1], fabs(d__1));
	    if (aik > amax1) {
		ipk1 = i__;
		amax1 = aik;
	    }
/* L200: */
	}
	ip[k1] = ipk1;
/*             EXCHANGE FOR K+1'TH ELIMINATION. */
	if (amax1 > *eps) {
	    if (ipk1 != k1) {
		i__2 = jmax1;
		for (j = k; j <= i__2; ++j) {
		    w = a[j - ipk1 + ipk1 * a_dim1];
		    a[j - ipk1 + ipk1 * a_dim1] = a[j - k1 + k1 * a_dim1];
		    a[j - k1 + k1 * a_dim1] = w;
/* L210: */
		}
	    }
/*             COMPUTE ALFA */
	    t = a[k1 * a_dim1];
	    i__2 = imax1;
	    for (i__ = k1 + 1; i__ <= i__2; ++i__) {
/* L220: */
		a[k1 - i__ + i__ * a_dim1] = -a[k1 - i__ + i__ * a_dim1] / t;
	    }
/*             MATRIX IS SINGULAR. */
	} else {
	    *ier = 1;
	    ip[k1] = k1;
	    i__2 = imax1;
	    for (i__ = k1 + 1; i__ <= i__2; ++i__) {
/* L230: */
		a[k1 - i__ + i__ * a_dim1] = 0.;
	    }
	    a[k1 * a_dim1] = *eps * .5;
/* *         WRITE(*,*)  '  (SUBR. BGLU2)  MATRIX IS SINGULAR AT K
 =', K1 */
/* *         RETURN */
	}

	t = a[k1 * a_dim1 - 1];
	i__2 = jmax1;
	for (j = k1 + 1; j <= i__2; ++j) {
	    wk1[j] = a[j - k + k * a_dim1];
	    a[j - k1 + k1 * a_dim1] += t * wk1[j];
/* L240: */
	    wk2[j] = a[j - k1 + k1 * a_dim1];
	}
/*             GAUSSIAN ELIMINATION FOR K-TH AND K+1'TH PROCESS. */
	i__2 = imax1;
	for (i__ = k1 + 1; i__ <= i__2; ++i__) {
	    t = a[k - i__ + i__ * a_dim1];
	    t1 = a[k1 - i__ + i__ * a_dim1];
	    i__3 = jmax1;
	    for (j = k1 + 1; j <= i__3; ++j) {
		a[j - i__ + i__ * a_dim1] = a[j - i__ + i__ * a_dim1] + t * 
			wk1[j] + t1 * wk2[j];
/* L310: */
	    }
/* L300: */
	}
/* L100: */
    }
    return 0;
} /* bglu2_ */


static int bgsv4c_(double *a, int *n, int *ml, int *
	mu, int *m, double *b, int *ip)
{
    /* System generated locals */
    int a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static int nend, i__, j, k, l;
    static double w;
    static int k1;
    static double s0, s1, bkl;
    static int ipk;
    static double bkl1;
    static int ipk1;


/*        BGSV4C */
/*               COPYRIGHT : K.MURATA, H.HASEGAWA, OCT.  4 1991 V.1 */

/*               SOLVES SIMULTANEOUS LINEAR EQUATIONS WITH MULTIPLE */
/*               RIGHT HAND SIDE BY GAUSSIAN ELIMINATION METHOD FOR */
/*               GENERAL BAND MATRIX. */

/*        INPUT - - */
/*             A(-ML-1:MU+ML+1,N+1) */
/*                      R *8  : RESULT OF GAUSSIAN ELIMINATION. */
/*             N        I *4  : ORDER OF MATRIX. */
/*             ML       I *4  : LOWER BAND WIDTH. */
/*             MU       I *4  : UPPER BAND WIDTH. */
/*             M        I *4  : NUMBER OF RIGHT HAND SIDE. */
/*             B(N+1,M) R *8  : 2-DIM. ARRAY CONTAINING THE RIGHT HAND */
/*                              SIDE VECTORS. */
/*             IP(N+1)  I *4  : PIVOT NUMBER. */
/*        OUTPUT - - */
/*             B(N+1,M)       : SOLUTION. */

/*             FORWARD ELIMINATION PROCESS */
    /* Parameter adjustments */
    b_dim1 = *n + 1;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    a_dim1 = *mu + *ml + 1 - (-(*ml) - 1) + 1;
    a_offset = -(*ml) - 1 + a_dim1;
    a -= a_offset;
    --ip;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; k += 2) {
/*             EXCHANGE FOR K-TH ELIMINATION. */
	ipk = ip[k];
	if (ipk != k) {
	    i__2 = *m;
	    for (l = 1; l <= i__2; ++l) {
		w = b[ipk + l * b_dim1];
		b[ipk + l * b_dim1] = b[k + l * b_dim1];
/* L110: */
		b[k + l * b_dim1] = w;
	    }
	}
/*             EXCHANGE FOR K+1'TH ELIMINATION. */
	k1 = k + 1;
	ipk1 = ip[k1];
	if (ipk1 != k1) {
	    i__2 = *m;
	    for (l = 1; l <= i__2; ++l) {
		w = b[ipk1 + l * b_dim1];
		b[ipk1 + l * b_dim1] = b[k1 + l * b_dim1];
/* L120: */
		b[k1 + l * b_dim1] = w;
	    }
	}
/*             GAUSSIAN ELIMINATION FOR K-TH AND K+1'TH PROCESS. */
	i__2 = *m;
	for (l = 1; l <= i__2; ++l) {
	    b[k1 + l * b_dim1] += a[k1 * a_dim1 - 1] * b[k + l * b_dim1];
	    bkl = b[k + l * b_dim1];
	    bkl1 = b[k1 + l * b_dim1];
/* Computing MIN */
	    i__4 = k1 + *ml;
	    i__3 = MIN(i__4,*n);
	    for (i__ = k1 + 1; i__ <= i__3; ++i__) {
/* L140: */
		b[i__ + l * b_dim1] = b[i__ + l * b_dim1] + a[k - i__ + i__ * 
			a_dim1] * bkl + a[k1 - i__ + i__ * a_dim1] * bkl1;
	    }
/* L130: */
	}
/* L100: */
    }
/*             BACKWARD SUBSTITUTION PROCESS */
    if (*n % 2 == 0) {
	nend = *n;
    } else {
	nend = *n + 1;
	i__1 = *m;
	for (l = 1; l <= i__1; ++l) {
/* L199: */
	    b[*n + 1 + l * b_dim1] = 0.;
	}
    }
    for (k = nend; k >= 1; k += -2) {
	i__1 = *m;
	for (l = 1; l <= i__1; ++l) {
	    s1 = -b[k + l * b_dim1];
	    s0 = -b[k - 1 + l * b_dim1];
/* Computing MIN */
	    i__3 = k + *ml + *mu;
	    i__2 = MIN(i__3,*n);
	    for (j = k + 1; j <= i__2; ++j) {
		s1 += a[j - k + k * a_dim1] * b[j + l * b_dim1];
/* L210: */
		s0 += a[j - k + 1 + (k - 1) * a_dim1] * b[j + l * b_dim1];
	    }
	    b[k + l * b_dim1] = -s1 / a[k * a_dim1];
	    b[k - 1 + l * b_dim1] = (-s0 - a[(k - 1) * a_dim1 + 1] * b[k + l *
		     b_dim1]) / a[(k - 1) * a_dim1];
/* L205: */
	}
/* L200: */
    }
    return 0;
} /* bgsv4c_ */


static int bgslv4_(double *a, int *n, int *ml, int *
	mu, double *b, int *ip)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static int nend, i__, j, k;
    static double t, w;
    static int k1;
    static double s0, s1, t1;
    static int ipk, ipk1;


/*        BGSLV4 */
/*               COPYRIGHT : H.HASEGAWA, OCT.  4 1991 V.1 */

/*               SOLVES SIMULTANEOUS LINEAR EQUATIONS */
/*               BY GAUSSIAN ELIMINATION METHOD FOR GENERAL BAND MATRIX. 
*/

/*        INPUT - - */
/*             A(-ML-1:MU+ML+1,N+1) */
/*                      R *8  : RESULT OF GAUSSIAN ELIMINATION. */
/*             N        I *4  : ORDER OF MATRIX. */
/*             ML       I *4  : LOWER BADN WIDTH. */
/*             MU       I *4  : UPPER BAND WIDTH. */
/*             B(N+1)   R *8  : 1-DIM. ARRAY CONTAINING THE RIGHT HAND */
/*                              SIDE VECTOR. */
/*             IP(N+1)  I *4  : PIVOT NUMBER. */
/*        OUTPUT - - */
/*             B(N+1)         : SOLUTION. */

/*             FORWARD ELIMINATION PROCESS */
    /* Parameter adjustments */
    a_dim1 = *mu + *ml + 1 - (-(*ml) - 1) + 1;
    a_offset = -(*ml) - 1 + a_dim1;
    a -= a_offset;
    --b;
    --ip;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; k += 2) {
/*             EXCHANGE FOR K-TH ELIMINATION. */
	ipk = ip[k];
	if (ipk != k) {
	    w = b[ipk];
	    b[ipk] = b[k];
	    b[k] = w;
	}
/*             EXCHANGE FOR K+1'TH ELIMINATION. */
	k1 = k + 1;
	ipk1 = ip[k1];
	if (ipk1 != k1) {
	    w = b[ipk1];
	    b[ipk1] = b[k1];
	    b[k1] = w;
	}

	b[k1] += a[k1 * a_dim1 - 1] * b[k];
/*             GAUSSIAN ELIMINATION FOR K-TH AND K+1'TH PROCESS. */
	t = b[k];
	t1 = b[k1];
/* Computing MIN */
	i__3 = k1 + *ml;
	i__2 = MIN(i__3,*n);
	for (i__ = k1 + 1; i__ <= i__2; ++i__) {
/* L110: */
	    b[i__] = b[i__] + a[k - i__ + i__ * a_dim1] * t + a[k1 - i__ + 
		    i__ * a_dim1] * t1;
	}
/* L100: */
    }
/*             BACKWARD SUBSTITUTION PROCESS */
    if (*n % 2 == 0) {
	nend = *n;
    } else {
	nend = *n + 1;
	b[*n + 1] = 0.;
    }
    for (k = nend; k >= 1; k += -2) {
	s1 = -b[k];
	s0 = -b[k - 1];
/* Computing MIN */
	i__2 = k + *mu + *ml;
	i__1 = MIN(i__2,*n);
	for (j = k + 1; j <= i__1; ++j) {
	    s1 += a[j - k + k * a_dim1] * b[j];
/* L210: */
	    s0 += a[j - k + 1 + (k - 1) * a_dim1] * b[j];
	}
	b[k] = -s1 / a[k * a_dim1];
	b[k - 1] = (-s0 - a[(k - 1) * a_dim1 + 1] * b[k]) / a[(k - 1) * 
		a_dim1];
/* L200: */
    }
    return 0;
} /* bgslv4_ */


static int yax_(int *n, int *m1, double *a, double *x, double *ax)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static int i__, j;
    static double s, t;


/*        YAX */
/*               COPYRIGHT : H.HASEGAWA, OCT.  4 1991 V.1 */

/*               COMPUTE A*X FOR UPPER TRIANGULAR BAND MATRIX A. */

/*        INPUT - - */
/*             A(0:M1,N) */
/*                      R *8  : 2-DIM. ARRAY CONTAINING UPPER POSITION OF 
*/
/*                              SYMMETRIC BAND MATRIX. */
/*             X(N)     R *8  : 1-DIM. ARRAY CONTAINING X. */
/*             N        I *4  : ORDER OF MATRIX. */
/*             M1       I *4  : UPPER BAND WIDTH. */
/*        OUTPUT - - */
/*             AX(N)    R *8  : 1-DIM. ARRAY CONTAINING A*X. */


    /* Parameter adjustments */
    a_dim1 = *m1 + 1;
    a_offset = a_dim1;
    a -= a_offset;
    --x;
    --ax;

    /* Function Body */
    for (i__ = *n; i__ >= 1; --i__) {
	s = a[i__ * a_dim1] * x[i__];
	t = x[i__];
/* Computing MIN */
	i__2 = i__ + *m1;
	i__1 = MIN(i__2,*n);
	for (j = i__ + 1; j <= i__1; ++j) {
	    ax[j] += a[j - i__ + i__ * a_dim1] * t;
/* L120: */
	    s += a[j - i__ + i__ * a_dim1] * x[j];
	}
	ax[i__] = s;
/* L100: */
    }
    return 0;
} /* yax_ */


static int mortho_(int *n, int *m1, int *ne, double 
	*v, double *eig, double *xm, double *vtvl2, double *
	gapt, double *w)
{
    /* Format strings */
    static char fmt_600[] = "(\0021\002/12x,\002*** \002,a25,\002 ***\002/)";
    static char fmt_610[] = "(\002 \002,2i4,2(5x,d10.2,a1))";

    /* System generated locals */
    int v_dim1, v_offset, xm_dim1, xm_offset, i__1, i__2, i__3;
    double d__1;

    /* Builtin functions */
/*
    int s_wsfe(cilist *), do_fio(int *, char *, ftnlen), e_wsfe(), 
	    s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen), 
	    e_wsle();
*/

    /* Local variables */
    static int nvtv, i__, j, k;
    static double gap;
    static double vtv;


/*        MORTHO */
/*               COPYRIGHT : K.MURATA, H.HASEGAWA, OCT.  4 1991 V.1 */

/*               CHECKS M-ORTHOGONALITY OF EIGENVECTORS */


    /* Parameter adjustments */
    v_dim1 = *n;
    v_offset = v_dim1 + 1;
    v -= v_offset;
    xm_dim1 = *m1 + 1;
    xm_offset = xm_dim1;
    xm -= xm_offset;
    --eig;
    --w;

    /* Function Body */
    nvtv = 0;
    i__1 = *ne;
    for (i__ = 1; i__ <= i__1; ++i__) {
	yax_(n, m1, &xm[xm_offset], &v[i__ * v_dim1 + 1], &w[1]);
	i__2 = *ne;
	for (j = i__ + 1; j <= i__2; ++j) {
	    gap = (d__1 = eig[j] - eig[i__], fabs(d__1));
	    if (gap <= *gapt) {
		vtv = 0.;
		i__3 = *n;
		for (k = 1; k <= i__3; ++k) {
/* L130: */
		    vtv += w[k] * v[k + j * v_dim1];
		}
		if (fabs(vtv) > *vtvl2) {
		    if (nvtv == 0) {
			fprintf(stderr, "--- error ---\n");
/*
			s_wsfe(&io___224);
			do_fio(&c__1, "    M-ORTHOGONALITY     ", 24L);
			e_wsfe();
*/
		    }
		    ++nvtv;
			fprintf(stderr, "--- error ---\n");
/*
		    s_wsfe(&io___225);
		    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(int));
		    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(int));
		    do_fio(&c__1, (char *)&vtv, (ftnlen)sizeof(double));
		    do_fio(&c__1, "*", 1L);
		    e_wsfe();
*/
		}
	    }
/* L120: */
	}
/* L110: */
    }
    if (nvtv != 0) {
		fprintf(stderr, " --- error ---\n");
/*
	s_wsle(&io___226);
	do_lio(&c__9, &c__1, " ", 1L);
	e_wsle();
	s_wsle(&io___227);
	do_lio(&c__9, &c__1, " *** TOTAL     : ", 17L);
	do_lio(&c__3, &c__1, (char *)&nvtv, (ftnlen)sizeof(int));
	e_wsle();
*/
    }
    return 0;
} /* mortho_ */

