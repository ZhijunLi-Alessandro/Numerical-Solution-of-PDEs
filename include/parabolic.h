/**
 * @file parabolic.h
 * @brief API for parabolic/heat equation matrix and RHS assembly helpers.
 *
 * This header exposes typedefs and assembly routines used by explicit and
 * ADI finite-difference solvers for the parabolic PDE examples in the
 * repository.
 * 
 * @author Li Zhijun
 * @date 2025-12-02
 */
#ifndef PARABOLIC_H
#define PARABOLIC_H
#include "grid.h"
#include "csr.h"

/**
 * @brief Function type for parabolic source terms f(x,y,t;hx,hy).
 *
 * The source term receives the spatial coordinates, current time, and the
 * local cell sizes `hx` and `hy` in case the source depends on the mesh.
 */
typedef double (*parabolic_source_term)(double x, double y, double t, double hx, double hy);

/**
 * @brief Function type for Dirichlet boundary values used by parabolic assembly.
 *
 * The callback receives (x,y,t) and an integer `boundary_type` which allows
 * callers to distinguish between multiple boundary segments.
 */
typedef double (*parabolic_Dirichlet_boundary)(double x, double y, double t, int boundary_type);

/**
 * @brief Assemble the system matrix for an explicit parabolic time-step.
 *
 * Produces a CSR matrix representing the update operator for the explicit
 * time integration (I - tau*L) or similar form depending on implementation.
 *
 * @param grid Pointer to the Grid2D structure describing the mesh and indexing.
 * @param tau Time-step size.
 * @return Pointer to a newly allocated SparseCSR matrix. Caller owns and must
 *         free the returned matrix using `freeSparseCSR`.
 */
SparseCSR* assemble_Matrix_Parabolic_Explicit(Grid2D* grid, double tau);

/**
 * @brief Assemble the right-hand side vector for a parabolic time step.
 *
 * This fills the provided array `b` with contributions from the source term
 * and Dirichlet boundary values for the current time `t` and step length
 * `tau`.
 *
 * @param grid Pointer to grid structure.
 * @param f Source term callback.
 * @param compute_boundary_value Dirichlet boundary callback.
 * @param b Pre-allocated array of length `grid->n_active` to receive RHS values.
 * @param t Current time.
 * @param tau Time-step size.
 */
void assemble_RHS_Parabolic(Grid2D* grid, parabolic_source_term f, parabolic_Dirichlet_boundary compute_boundary_value, double *b, double t, double tau);

/**
 * @brief Assemble ADI operator matrices used by the Peaceman–Rachford ADI method.
 *
 * Returns an array of four `SparseCSR*` pointers representing the split
 * operators; caller owns the returned array and the matrices and must free
 * them when no longer needed.
 *
 * @param grid Pointer to the mesh/grid structure.
 * @param tau Time-step size.
 * @return Array of four pointers to SparseCSR matrices.
 */
SparseCSR** assemble_Matrix_Parabolic_ADI(Grid2D* grid, double tau);

#endif