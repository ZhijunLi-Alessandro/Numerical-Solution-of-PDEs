#ifndef PARABOLIC_H
#define PARABOLIC_H
#include "grid.h"
#include "csr.h"

typedef double (*parabolic_source_term)(double x, double y, double t, double hx, double hy);
typedef double (*parabolic_Dirichlet_boundary)(double x, double y, double t, int boundary_type);

SparseCSR* assemble_Matrix_Parabolic_Explicit(Grid2D* grid, double tau);
void assemble_RHS_Parabolic(Grid2D* grid, parabolic_source_term f, parabolic_Dirichlet_boundary compute_boundary_value, double *b, double t, double tau);
SparseCSR** assemble_Matrix_Parabolic_ADI(Grid2D* grid, double tau);

#endif