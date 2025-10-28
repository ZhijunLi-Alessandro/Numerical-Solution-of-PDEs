/**
 * @file poisson2d.h
 * @brief Header file for assembling 2D Poisson equation matrix and RHS vector.
 * 
 * This header file declares functions to assemble the sparse matrix
 * and right-hand side vector for the 2D Poisson equation on a given grid.
 * @see poisson2d.c
 * @author Li Zhijun
 * @date 2025-10-22
 */
# ifndef POISSON2D_H
# define POISSON2D_H
# include "grid.h"
# include "csr.h"

typedef double (*f_func)(double, double);
typedef double (*boundary_func)(double, double, int);
typedef double (*normal_func)(int);

SparseCSR* assemble_Matrix_Dirichlet(Grid2D* grid);
double* assemble_RHS_Dirichlet(Grid2D* grid, f_func f, boundary_func compute_boundary_value);
SparseCSR* assemble_Matrix_Neumann(Grid2D* grid, normal_func get_normal);
double* assemble_RHS_Neumann(Grid2D* grid, f_func f, boundary_func compute_boundary_value, f_func get_exact);

# endif