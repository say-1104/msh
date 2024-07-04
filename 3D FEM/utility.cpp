#include "fem3d.h"

#define TYPE2

void createS3d(DataTable *data) {
  long long int i, kk, it;
  double aa[6], bb[6], cc[6];
  long long int n1, n2;
  static long long int edge[6][2] = { { 0, 1 }, { 0, 2 }, { 0, 3 }, 
				      { 1, 2 }, { 3, 1 },{ 2, 3 } };
#ifdef TYPE1
  double x[10], y[10], z[10];
#endif

  FEM *fem = &(data->fem);

  ALLOCATION(fem->s3d, long long int, fem->n_edge);
  if (fem->order != 1) {
    for (i = 0; i < fem->n_edge; i++)
      fem->s3d[i] = 1;
  } else {
    for (i = 0; i < fem->n_edge; i++)
      fem->s3d[i] = 0;
    for (kk = 0; kk < fem->ne; kk++) {
#ifdef TYPE1
      for (i = 0; i < fem->n_node; i++) {
	it = fem->element[kk].kn[i];
	x[i] = fem->node[it].x;
	y[i] = fem->node[it].y;
	z[i] = fem->node[it].z;
      }
      for (i = 0; i < 6; i++) {
	aa[i] = x[edge[i][1]]-x[edge[i][0]];
	bb[i] = y[edge[i][1]]-y[edge[i][0]];
	cc[i] = z[edge[i][1]]-z[edge[i][0]];
      }
#endif
      for (i = 0; i < 6; i++) {
#ifdef TYPE1
	it = fem->element[kk].kk[i];
	n1 = fem->element[kk].kn[edge[i][0]];
	n2 = fem->element[kk].kn[edge[i][1]];
#endif
#ifdef TYPE2
	it = fem->element[kk].kk[i];
	n1 = fem->element[kk].kn[edge[i][0]];
	n2 = fem->element[kk].kn[edge[i][1]];
	aa[i] = fem->node[n2].x - fem->node[n1].x;
	bb[i] = fem->node[n2].y - fem->node[n1].y;
	cc[i] = fem->node[n2].z - fem->node[n1].z;
#endif
	if (cc[i] > ZERO || (std::abs(cc[i]) <= ZERO && bb[i] > ZERO)
	    || (std::abs(cc[i]) <= ZERO && std::abs(bb[i]) <= ZERO
		&& aa[i] > ZERO)) {
	  if (n2 > n1) {
	    if (fem->s3d[it] == -1)
	      fprintf(stderr, "s3d\n");
	    fem->s3d[it] = 1;
	  } else {
	    if (fem->s3d[it] == 1)
	      fprintf(stderr, "s3d\n");
	    fem->s3d[it] = -1;
	  }
	} else {
	  if (n2 > n1) {
	    if (fem->s3d[it] == 1)
	      fprintf(stderr, "s3d\n");
	    fem->s3d[it] = -1;
	  } else {
	    if (fem->s3d[it] == -1)
	      fprintf(stderr, "s3d\n");
	    fem->s3d[it] = 1;
	  }
	}
      }
    }
  }
  
  fprintf(stderr, "- create s3d -\n");
  
  return;
}

void createS2d(DataTable *data, PortInfo *port, long long int n_port) {
  long long int kk, i, it, n1, n2;
  double b[N_NODE], c[N_NODE];
  double cosA, sinA;
#ifdef TYPE1
  double x[N_NODE], y[N_NODE], z[N_NODE], xr[N_NODE];
#endif
#ifdef TYPE2
  double x1, x2, y1, y2, z1, z2, xr1, xr2;
#endif

  FEM *fem = &(data->fem);
  
  ALLOCATION(port->s2d, long long int, port->n_edge+port->np);
  ALLOCATION(port->s23d, long long int, port->n_edge+port->np);
  if (fem->order != 1) {
    for (i = 0; i < port->n_edge + port->np; i++)
      port->s2d[i] = 1;
  } else {
    for (i = 0; i < port->n_edge + port->np; i++) {
      port->s2d[i] = 0;
      // port->s2d[i] = 1;
    }
    
    cosA = cos(port->angle * M_PI / 180.0);
    sinA = sin(port->angle * M_PI / 180.0);
    
    for (kk = 0; kk < port->ne; kk++) {
#ifdef TYPE1
      for (i = 0; i < port->l_node; i++) {
	it = port->element[kk].kn[i];
	x[i] = fem->node[it].x - port->x;
	y[i] = fem->node[it].y - port->y;
	z[i] = fem->node[it].z - port->z;
	xr[i] = (x[i] * cosA - z[i] * sinA);
      }
      b[0] = y[1] - y[2], b[1] = y[2] - y[0], b[2] = y[0] - y[1];
      c[0] = xr[2] - xr[1], c[1] = xr[0] - xr[2], c[2] = xr[1] - xr[0];
#endif
      for (i = 0; i < 3; i++) {
#ifdef TYPE1
	it = port->toLe[port->element[kk].ke[(i+1)%3]];
	n1 = port->element[kk].kn[(i+1)%3];
	n2 = port->element[kk].kn[(i+2)%3];
#endif
#ifdef TYPE2
	it = port->toLe[port->element[kk].ke[i]];
	n1 = port->element[kk].kn[(i) % 3];
	n2 = port->element[kk].kn[(i + 1) % 3];
	x1 = fem->node[n1].x - port->x;
	y1 = fem->node[n1].y - port->y;
	z1 = fem->node[n1].z - port->z;
	xr1 = (x1 * cosA - z1 * sinA);
	x2 = fem->node[n2].x - port->x;
	y2 = fem->node[n2].y - port->y;
	z2 = fem->node[n2].z - port->z;
	xr2 = (x2 * cosA - z2 * sinA);
	b[i] = y1 - y2;
	c[i] = xr2 - xr1;
#endif
	if (b[i] < -1.0e-4 || (std::abs(b[i]) < 1.0e-4 && c[i] > 0.0)) {
	  if (n2 > n1) {
	    if (port->s2d[it] == -1)
	      fprintf(stderr, "s2d\n");
	    port->s2d[it] = 1;
	  } else {
	    if (port->s2d[it] == 1)
	      fprintf(stderr, "s2d\n");
	    port->s2d[it] = -1;
	  }
	} else {
	  if (n2 > n1) {
	    if (port->s2d[it] == 1)
	      fprintf(stderr, "s2d\n");
	    port->s2d[it] = -1;
	  } else {
	    if (port->s2d[it] == -1)
	      fprintf(stderr, "s2d\n");
	    port->s2d[it] = 1;
	  }
	}
      }
    }
    for (i = 0; i < fem->n_edge; i++) {
      if (port->toLe[i] >= 0 && port->toLe[i] < port->nr) {
	port->s23d[port->toLe[i]] = fem->s3d[i]
	  * port->s2d[port->toLe[i]];
      }
    }
  }
  
  fprintf(stderr, "- create s2d -\n");
  
  return;
}
