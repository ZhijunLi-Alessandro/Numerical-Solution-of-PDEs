/**
 * @file poisson2d.c
 * @brief Implementation of assembling the matrix and right-hand side (RHS) for the 
 *        Poisson equation with Dirichlet and Neumann boundary conditions.
 * 
 * This file includes functions to construct the finite difference discretization of 
 * the Poisson equation, handling both Dirichlet (essential) and Neumann (natural) 
 * boundary conditions. The system matrix and RHS vector are assembled accordingly for 
 * numerical solution.
 * 
 * @author Li Zhijun
 * @date 2025-10-22
 */
#include <stdlib.h>
#include <math.h>
#include "poisson2d.h"

SparseCSR* assemble_Matrix_Dirichlet(Grid2D* grid) {
    SparseCSR *matrix = createSparseCSR(grid->n_active, grid->n_active, 5 * grid->n_active);
    int *row_ptr = matrix->row_ptr;
    int *col_ind = matrix->col_ind;
    double *values = matrix->values;

    int idx = 0;
    row_ptr[0] = 0;

    for (int i = 0; i < grid->n_active; i++) {
        int gi = grid->id_i[i];
        int gj = grid->id_j[i];
        double xi = grid->x[gi];
        double yj = grid->y[gj];

        // Interior point
        if (grid->region[gi][gj] == 1) {
            // Center
            int row = i;
            int col = i;
            col_ind[idx] = col;
            values[idx] = 4.0;
            idx++;

            // Left
            col = grid->id_map[gi - 1][gj];
            col_ind[idx] = col;
            values[idx] = -1.0;
            idx++;

            // Right
            col = grid->id_map[gi + 1][gj];
            col_ind[idx] = col;
            values[idx] = -1.0;
            idx++;

            // Down
            col = grid->id_map[gi][gj - 1];
            col_ind[idx] = col;
            values[idx] = -1.0;
            idx++;

            // Up
            col = grid->id_map[gi][gj + 1];
            col_ind[idx] = col;
            values[idx] = -1.0;
            idx++;
        }

        // Boundary point
        else {
            int row = i;
            int col = i;
            col_ind[idx] = col;
            values[idx] = 1.0;
            idx++;
        }

        row_ptr[i + 1] = idx;
    }
    matrix->nnz = idx; // Update nnz
    return matrix;
}

double* assemble_RHS_Dirichlet(Grid2D* grid, f_func f, boundary_func compute_boundary_value) {
    double *b = (double *)malloc(grid->n_active * sizeof(double));
    double h = grid->hx; // Assuming hx = hy
    for (int i = 0; i < grid->n_active; i++) {
        int gi = grid->id_i[i];
        int gj = grid->id_j[i];
        double xi = grid->x[gi];
        double yj = grid->y[gj];

        // Interior point
        if (grid->region[gi][gj] == 1) {
            b[i] = f(xi, yj) * h * h; // Source term
        }
        // Boundary point
        else {
            // Dirichlet conditions based on boundary type
            b[i] = compute_boundary_value(xi, yj, grid->region[gi][gj]);
        }
    }
    return b;
}

SparseCSR* assemble_Matrix_Neumann(Grid2D* grid, normal_func get_normal) {
    SparseCSR *matrix = createSparseCSR(grid->n_active, grid->n_active, 5 * grid->n_active);
    int *row_ptr = matrix->row_ptr;
    int *col_ind = matrix->col_ind;
    double *values = matrix->values;

    int idx = 0, flag = 1;
    row_ptr[0] = 0;

    for (int i = 0; i < grid->n_active; i++) {
        int gi = grid->id_i[i];
        int gj = grid->id_j[i];
        double xi = grid->x[gi];
        double yj = grid->y[gj];

        // Interior point
        if (grid->region[gi][gj] == 1) {
            if (flag) {
                col_ind[idx] = i;
                values[idx] = 1.0;
                idx++;
                flag = 0;
            } else {
                // Center
                int row = i;
                int col = i;
                col_ind[idx] = col;
                values[idx] = 4.0;
                idx++;

                // Left
                col = grid->id_map[gi - 1][gj];
                col_ind[idx] = col;
                values[idx] = -1.0;
                idx++;

                // Right
                col = grid->id_map[gi + 1][gj];
                col_ind[idx] = col;
                values[idx] = -1.0;
                idx++;

                // Down
                col = grid->id_map[gi][gj - 1];
                col_ind[idx] = col;
                values[idx] = -1.0;
                idx++;

                // Up
                col = grid->id_map[gi][gj + 1];
                col_ind[idx] = col;
                values[idx] = -1.0;
                idx++;
            }
        }

        // Boundary point
        else {
            double alpha = get_normal(grid->region[gi][gj]);
            int row = i;
            int col = i;
            col_ind[idx] = col;
            values[idx] = fabs(sin(alpha)) + fabs(cos(alpha));
            idx++;

            if (fabs(sin(alpha)) > 1e-12) {
                if (sin(alpha) > 0) { // Bottom boundary
                    col = grid->id_map[gi][gj - 1];
                    col_ind[idx] = col;
                    values[idx] = -sin(alpha);
                    idx++;
                } else { // Top boundary
                    col = grid->id_map[gi][gj + 1];
                    col_ind[idx] = col;
                    values[idx] = sin(alpha);
                    idx++;
                }
            }
            if (fabs(cos(alpha)) > 1e-12) {
                if (cos(alpha) > 0) { // Left boundary
                    col = grid->id_map[gi - 1][gj];
                    col_ind[idx] = col;
                    values[idx] = -cos(alpha);
                    idx++;
                } else { // Right boundary
                    col = grid->id_map[gi + 1][gj];
                    col_ind[idx] = col;
                    values[idx] = cos(alpha);
                    idx++;
                }
            }
        }
        row_ptr[i + 1] = idx;
    }
    matrix->nnz = idx; // Update nnz
    return matrix;
}

double* assemble_RHS_Neumann(Grid2D* grid, f_func f, boundary_func compute_boundary_value, f_func get_exact) {
    double *b = (double *)malloc(grid->n_active * sizeof(double));
    double h = grid->hx; // Assuming hx = hy
    int flag = 1;
    for (int i = 0; i < grid->n_active; i++) {
        int gi = grid->id_i[i];
        int gj = grid->id_j[i];
        double xi = grid->x[gi];
        double yj = grid->y[gj];

        // Interior point
        if (grid->region[gi][gj] == 1) {
            if (flag) {
                b[i] = get_exact(xi, yj);
                flag = 0;
            } else {
                b[i] = f(xi, yj) * h * h; // Source term
            }
        }
        // Boundary point
        else {
            // Neumann conditions based on boundary type
            b[i] = compute_boundary_value(xi, yj, grid->region[gi][gj]) * h;
        }
    }
    return b;
}