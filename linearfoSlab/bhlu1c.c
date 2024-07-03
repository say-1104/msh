#include "fem.h"

	static int c__9 = 9;
	static int c__1 = 1;
	static int c__3 = 3;

extern int bhlu1c_(double_complex *a, double_complex *ac, int *n, 
	int *ml, int *mu, double *eps, double_complex *wk, int *ip, int *ier)
{
    /* System generated locals */
    int a_dim1, a_offset, ac_dim1, ac_offset,
		i__1, i__2, i__3, i__4, i__5, i__6;
    double d__1, d__2;
    double_complex z__1, z__2;

    /* Local variables */
    static double amax;
    static int jmax, i__, j, k;
    static double_complex t, w;
    static double aik;
    static int ipk;

/*             LEFT HAND SIDE */
    /* Parameter adjustments */
    ac_dim1 = *ml;
    ac_offset = ac_dim1 + 1;
    ac -= ac_offset;
    a_dim1 = *mu + *ml + 1;
    a_offset = a_dim1;
    a -= a_offset;
    --wk;
    --ip;

    /* Function Body */
    if (*eps < 0.0) {
      *eps = 3.52e-15;
    }
    if (*n <= 0 || *ml <= 0 || *mu <= 0 || *ml >= *n || *mu >= *n) {
      *ier = 3;
      fprintf(stderr, "  (SUBR. BHLU1)  INVALID ARGUMENT.  ML, MU, N =",
	      *ml, *mu, *n);
      return 0;
    }
    
    *ier = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
      /* Computing MIN */
      i__2 = k + *mu + *ml;
      jmax = (i__2 < *n)? i__2: *n;
      /* FIND MAXIMUM ELEMENT IN THE K-TH COLUMN. */
      amax = (d__1 = real(a[k * a_dim1]), fabs(d__1))
	+(d__2 = imag(a[k * a_dim1]), fabs(d__2));
      ipk = k;
      /* Computing MIN */
      i__3 = k + *ml;
      i__2 = (i__3 < *n) ? i__3 : *n;
      for (i__ = k + 1; i__ <= i__2; ++i__) {
	aik = (d__1 = real(a[i__ * a_dim1]), fabs(d__1))
	  + (d__2 = imag(a[i__ * a_dim1]), fabs(d__2));
	if (aik > amax) {
	  ipk = i__;
	  amax = aik;
	}
      }
      ip[k] = ipk;
      
      if (amax > *eps) {
	if (ipk != k) {
	  /* Computing MIN */
	  i__3 = k + *mu + *ml;
	  i__2 = (i__3 < *n) ? i__3 : *n;
	  for (j = k; j <= i__2; ++j) {
	    i__3 = j - k + ipk * a_dim1;
	    w = a[i__3];
	    i__3 = j - k + ipk * a_dim1;
	    i__4 = j - k + k * a_dim1;
	    a[i__3] = a[i__4];
	    i__3 = j - k + k * a_dim1;
	    a[i__3] = w;
	  }
	}
	
	/* Computing MIN */
	i__3 = k + *mu + *ml;
	i__2 = (i__3 < *n) ? i__3 : *n;
	for (j = k + 1; j <= i__2; ++j) {
	  i__3 = j;
	  i__4 = j - k + k * a_dim1;
	  wk[i__3] = a[i__4];
	}
	/* COMPUTE ALFA AND PERFORM GAUSSIAN ELIMINATION. */
	/* Computing MIN */
	i__4 = k + *ml;
	i__3 = (i__4 < *n) ? i__4 : *n;
	for (i__ = k + 1; i__ <= i__3; ++i__) {
	  i__4 = i__ - k + k * ac_dim1;
	  i__2 = i__ * a_dim1;
	  z__2 = -a[i__2];
	  z__1 = z__2/a[k*a_dim1];
	  /* z_div(&z__1, &z__2, &a[k * a_dim1]); */
	  ac[i__4] = z__1;
	  i__4 = i__ * a_dim1;
	  z__2 = -a[i__4];
	  z__1 = z__2/a[k*a_dim1];
	  /* z_div(&z__1, &z__2, &a[k * a_dim1]); */
	  t = z__1;
	  /* Computing MIN */
	  i__2 = k + *mu + *ml;
	  i__4 = (i__2 < *n) ? i__2 : *n;
	  for (j = k + 1; j <= i__4; ++j) {
	    i__2 = j - k - 1 + i__ * a_dim1;
	    i__5 = j - k + i__ * a_dim1;
	    i__6 = j;
	    z__2 = t*wk[i__6];
	    z__1 = a[i__5] + z__2;
	    a[i__2] = z__1;
	  }
	  i__2 = jmax - k + i__ * a_dim1;
	  a[i__2] = 0.0;
	}
	/* MATRIX IS SINGULAR. */
      } else {
	*ier = 1;
	i__3 = k * a_dim1;
	a[i__3] = *eps;
	/* Computing MIN */
	i__2 = k + *ml;
	i__3 = (i__2 < *n) ? i__2 : *n;
	for (i__ = k + 1; i__ <= i__3; ++i__) {
	  i__2 = i__ - k + k * ac_dim1;
	  ac[i__2] = 0.0;
	  /* Computing MIN */
	  i__5 = k + *mu + *ml;
	  i__2 = (i__5 < *n) ? i__5 : *n;
	  for (j = k + 1; j <= i__2; ++j) {
	    i__5 = j - k - 1 + i__ * a_dim1;
	    i__6 = j - k + i__ * a_dim1;
	    a[i__5] = a[i__6];
	  }
	  i__5 = jmax - k + i__ * a_dim1;
	  a[i__5] = 0.0;
	}
	fprintf(stderr, "  (SUBR. BHLU1)  MATRIX IS SINGULAR AT K = %d\n", k);
	return 0;
      }
    }
    return 0;
}


extern int bhslv1c_(double_complex *a, double_complex *ac, int *n,
	 int *ml, int *mu, double_complex *b, int *ip)
{
    /* System generated locals */
    int a_dim1, a_offset, ac_dim1, ac_offset, i__1, i__2, i__3, i__4, i__5;
    double_complex z__1, z__2;

    /* Local variables */
    static int jmax, i__, j, k;
    static double_complex s, t, w;
    static int ipk;

	/* FORWARD ELIMINATION PROCESS */
    /* Parameter adjustments */
    ac_dim1 = *ml;
    ac_offset = ac_dim1 + 1;
    ac -= ac_offset;
    a_dim1 = *mu + *ml + 1;
    a_offset = a_dim1;
    a -= a_offset;
    --b;
    --ip;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
		ipk = ip[k];
		if (ipk != k) {
		    i__2 = ipk;
		    w = b[i__2];
		    i__2 = ipk;
		    i__3 = k;
		    b[i__2] = b[i__3];
		    i__2 = k;
		    b[i__2] = w;
		}
		/* GAUSSIAN ELIMINATION */
		i__2 = k;
		t = b[i__2];
		/* Computing MIN */
		i__3 = k + *ml;
		i__2 = (i__3 < *n) ? i__3 : *n;
		for (i__ = k + 1; i__ <= i__2; ++i__) {
	    	i__3 = i__;
		    i__4 = i__;
		    i__5 = i__ - k + k * ac_dim1;
			z__2 = ac[i__5]*t;
	    	z__1 = b[i__4] + z__2;
		    b[i__3] = z__1;
		}
    }
	/* BACKWARD SUBSTITUTION PROCESS */
    for (k = *n; k >= 1; --k) {
		i__1 = k;
		z__1 = -b[i__1];
		s = z__1;
		/* Computing MIN */
		i__1 = k + *mu + *ml;
		jmax = (i__1 < *n) ? i__1 : *n;
		i__1 = jmax;
		for (j = k + 1; j <= i__1; ++j) {
	    	i__3 = j - k + k * a_dim1;
		    i__4 = j;
			z__2 = a[i__3]*b[i__4];
		    z__1 = s + z__2;
		    s = z__1;
		}
		i__3 = k;
		z__2 = -s;
		/* z_div(&z__1, &z__2, &a[k * a_dim1]); */
		z__1 = z__2/a[k * a_dim1];
		b[i__3] = z__1;
    }
    return 0;
}
