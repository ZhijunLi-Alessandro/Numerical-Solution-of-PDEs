#ifndef PARABOLIC_H
#define PARABOLIC_H
#include "grid.h"
#include "csr.h"

SparseCSR* assemble_Matrix_Parabolic_Explicit(Grid2D* grid, double t);

#endif