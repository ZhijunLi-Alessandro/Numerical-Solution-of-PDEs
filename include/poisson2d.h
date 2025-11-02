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

/**
 * @brief Assemble the system matrix for the numerical solution of the Dirichlet problem.
 * @param grid Pointer to the grid structure which defines the solution region and boundary conditions.
 * 
 * @note User should also call assemble_RHS_Dirichlet() to assemble the right-hand side vector.
 * @see assmble_RHS_Dirichlet()
 */
SparseCSR* assemble_Matrix_Dirichlet(Grid2D* grid);

/**
 * @brief Assemble the right-hand side (RHS) vector for the numerical solution of the Dirichlet problem.
 * @param grid Pointer to the grid structure which defines the solution region and boundary conditions.
 * @param f Function pointer to calculate the source term of the Poisson function.
 * @param compute_boundary_value Function pointer to assigns values to boundaries.
 * @par Function signature requirements:
 * @code
 * double f(double x, double y)
 * @endcode
 *      - param x: X-coordinate of the grid point
 *      - param y: Y-coordinate of the grid point
 *      - return value: the source term value at point (x, y)
 * @code
 * double compute_boundary_value(double x, double y, int boundary_type)
 * @endcode
 *      - param x: X-coordinate of the grid point
 *      - param y: Y-coordinate of the grid point
 *      - param boundary_type: the region type. (defined by Grid2d->region)
 *      - return value: Interpolation of the boundary point's value
 * 
 * @note User should also call assemble_Matrix_Dirichlet() to assemble the system matrix
 * @see assemble_Matrix_Dirichlet()
 * @see Grid2D
 */
double* assemble_RHS_Dirichlet(Grid2D* grid, f_func f, boundary_func compute_boundary_value);

/**
 * @brief 
 */
SparseCSR* assemble_Matrix_Neumann(Grid2D* grid, normal_func get_normal);
double* assemble_RHS_Neumann(Grid2D* grid, f_func f, boundary_func compute_boundary_value, f_func get_exact);

# endif