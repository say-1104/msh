#include "optfem.h"

extern int allocMatrix(DataTable *data)
{
  int i, j, k, ll, mm;
  FEM *fem = &(data->fem);
  Port *port = &(data->port);
  CGtable *cg_table;

  ALLOCATION(fem->cg_table, CGtable, fem->nr);
  cg_table = fem->cg_table;
  makeCGmatrixTable(fem->element, fem->ne, fem->nr, cg_table, 6);
  ALLOCATION(fem->A, double_complex *, fem->nr)
  for (i = 0; i < fem->nr; i++) 
    ALLOCATION(fem->A[i], double_complex, cg_table[i].n);
 
  if (data->port.inputDirection == SINGLE_DIRECTION) {
    ALLOCATION(fem->A12, double_complex *, fem->nr);
    for (i = 0; i < fem->nr; i++) 
      ALLOCATION(fem->A12[i], double_complex, cg_table[i].n);

    ALLOCATION(fem->B12, double_complex *, fem->nr);
    for (i = 0; i < fem->nr; i++) 
      ALLOCATION(fem->B12[i], double_complex, cg_table[i].n);

  }

  ALLOCATION(fem->phi, double_complex, fem->np);
  ALLOCATION(fem->b, double_complex, fem->np);

  for (ll = 0; ll < data->par.numoflambda; ll++) {
    for (mm = 0; mm < data->par.num_mode; mm++) {
      for (k = 0; k < port->number; k++) {
        ALLOCATION(port->dataMemory[ll][mm][k].NN, double_complex *, port->dataMemory[ll][mm][k].np);
        for (i = 0; i < port->dataMemory[ll][mm][k].np; i++) {
          ALLOCATION(port->dataMemory[ll][mm][k].NN[i], double_complex, port->dataMemory[ll][mm][k].np);
        }
        ALLOCATION(port->dataMemory[ll][mm][k].NNy, double_complex *, port->dataMemory[ll][mm][k].np);
        for (i = 0; i < port->dataMemory[ll][mm][k].np; i++) {
          ALLOCATION(port->dataMemory[ll][mm][k].NNy[i], double_complex, port->dataMemory[ll][mm][k].np);
        }
      }

      for (k = 0; k < port->number; k++) {
        ALLOCATION(port->dataMemory[ll][mm][k].phi, double_complex, port->dataMemory[ll][mm][k].np);
      }
    }
  }
}
