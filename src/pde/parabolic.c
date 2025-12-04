/**
 * @file parabolic.c
 * @brief Assembly helpers for finite-difference parabolic (heat) solvers.
 *
 * This file implements routines to assemble sparse matrices and right-hand-side
 * vectors used by time-stepping schemes for the parabolic (heat) equation on the
 * repository's structured grids. Implementations include:
 *  - `assemble_Matrix_Parabolic_Explicit`: build the operator for an explicit
 *    time-step update on the active grid points.
 *  - `assemble_RHS_Parabolic`: fill the RHS vector using a provided source-term
 *    callback and Dirichlet boundary value callback.
 *  - `assemble_Matrix_Parabolic_ADI`: construct directional ADI split operators
 *    (arrays of CSR matrices) for alternating-direction implicit methods.
 *
 * The functions operate on the `Grid2D` structure and return newly allocated
 * `SparseCSR` objects (or arrays of them) which the caller must free when no
 * longer needed.
 *
 * @author Li Zhijun
 * @date 2025-12-02
 */
#include <stdlib.h>
#include <math.h>
#include "parabolic.h"

/**
 * @brief Assemble the sparse matrix for an explicit parabolic time-step.
 *
 * Constructs a CSR matrix containing the coefficients appropriate for an
 * explicit finite-difference update on the active points of the grid. The
 * returned matrix is newly allocated and must be freed using `freeSparseCSR`.
 *
 * @param grid Pointer to Grid2D describing the mesh and active indices.
 * @param tau Time-step size.
 * @return Pointer to the assembled SparseCSR matrix.
 */
SparseCSR* assemble_Matrix_Parabolic_Explicit(Grid2D* grid, double tau) {
    SparseCSR *matrix = createSparseCSR(grid->n_active, grid->n_active, 5 * grid->n_active);
    int *row_ptr = matrix->row_ptr;
    int *col_ind = matrix->col_ind;
    double *values = matrix->values;

    int idx = 0;
    row_ptr[0] = 0;

    double mu_x = tau / grid->hx / grid->hx;
    double mu_y = tau / grid->hy / grid->hy;

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
            values[idx] = 1.0 - 2 * (mu_x + mu_y);
            idx++;

            // Left
            col = grid->id_map[gi - 1][gj];
            col_ind[idx] = col;
            values[idx] = mu_x;
            idx++;

            // Right
            col = grid->id_map[gi + 1][gj];
            col_ind[idx] = col;
            values[idx] = mu_x;
            idx++;

            // Down
            col = grid->id_map[gi][gj - 1];
            col_ind[idx] = col;
            values[idx] = mu_y;
            idx++;

            // Up
            col = grid->id_map[gi][gj + 1];
            col_ind[idx] = col;
            values[idx] = mu_y;
            idx++;
        }

        // Boundary point
        else {
            int row = i;
            int col = i;
            col_ind[idx] = col;
            values[idx] = 0.0;
            idx++;
        }

        row_ptr[i + 1] = idx;
    }
    matrix->nnz = idx; // Update nnz
    return matrix;
}

/**
 * @brief Assemble the right-hand side vector for a parabolic step.
 *
 * For interior (active) points this fills the contribution from the source
 * term integrated over the time-step; for boundary points it fills the
 * Dirichlet boundary value via the provided callback.
 *
 * @param grid Pointer to Grid2D describing the mesh and active indices.
 * @param f Source term callback with signature `f(x,y,t,hx,hy)`.
 * @param compute_boundary_value Dirichlet boundary callback.
 * @param b Preallocated array of length grid->n_active to be filled with RHS.
 * @param t Current time.
 * @param tau Time-step size.
 */
void assemble_RHS_Parabolic(Grid2D* grid, parabolic_source_term f, parabolic_Dirichlet_boundary compute_boundary_value, double *b, double t, double tau) {
    // double *b = (double *)malloc(grid->n_active * sizeof(double));
    double hx = grid->hx;
    double hy = grid->hy;
    for (int i = 0; i < grid->n_active; i++) {
        int gi = grid->id_i[i];
        int gj = grid->id_j[i];
        double xi = grid->x[gi];
        double yj = grid->y[gj];

        // Interior point
        if (grid->region[gi][gj] == 1) {
            b[i] = f(xi, yj, (t - tau / 2), hx, hy) * tau; // Intergrated Source term!!
        }
        // Boundary point
        else {
            // Dirichlet conditions based on boundary type
            b[i] = compute_boundary_value(xi, yj, t, grid->region[gi][gj]);
        }
    }
    // return b;
}

/**
 * @brief Assemble the split ADI matrices used by an ADI parabolic solver.
 *
 * Returns four CSR matrices (allocated) representing directional operators
 * used in a Peaceman–Rachford / Douglas ADI step. The caller receives an
 * array of four `SparseCSR*` pointers and is responsible for freeing them.
 *
 * @param grid Pointer to the Grid2D structure.
 * @param tau Time-step size.
 * @return Array of 4 pointers to newly allocated SparseCSR matrices.
 */
SparseCSR** assemble_Matrix_Parabolic_ADI(Grid2D* grid, double tau) {
    SparseCSR **ADI_matrixs = (SparseCSR **)malloc(4 * sizeof(SparseCSR*));

    double mu_x = tau / grid->hx / grid->hx;
    double mu_y = tau / grid->hy / grid->hy;

    SparseCSR *plus_delta_y = createSparseCSR(grid->n_active, grid->n_active, 3 * grid->n_active);
    int *row_ptr = plus_delta_y->row_ptr;
    int *col_ind = plus_delta_y->col_ind;
    double *values = plus_delta_y->values;

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
            values[idx] = 1.0 - mu_y;
            idx++;

            // Down
            col = grid->id_map[gi][gj - 1];
            col_ind[idx] = col;
            values[idx] = mu_y / 2;
            idx++;

            // Up
            col = grid->id_map[gi][gj + 1];
            col_ind[idx] = col;
            values[idx] = mu_y / 2;
            idx++;
        }
        // Boundary point
        else {
            int row = i;
            int col = i;
            col_ind[idx] = col;
            values[idx] = 0.0;
            idx++;
        }
        row_ptr[i + 1] = idx;
    }
    plus_delta_y->nnz = idx; // Update nnz
    ADI_matrixs[0] = plus_delta_y;

    SparseCSR *minus_delta_x = createSparseCSR(grid->n_active, grid->n_active, 3 * grid->n_active);
    row_ptr = minus_delta_x->row_ptr;
    col_ind = minus_delta_x->col_ind;
    values = minus_delta_x->values;

    idx = 0;
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
            values[idx] = 1.0 + mu_x;
            idx++;

            // Left
            col = grid->id_map[gi - 1][gj];
            col_ind[idx] = col;
            values[idx] = - mu_x / 2;
            idx++;

            // Right
            col = grid->id_map[gi + 1][gj];
            col_ind[idx] = col;
            values[idx] = - mu_x / 2;
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
    minus_delta_x->nnz = idx; // Update nnz
    ADI_matrixs[1] = minus_delta_x;
    
    SparseCSR *plus_delta_x = createSparseCSR(grid->n_active, grid->n_active, 3 * grid->n_active);
    row_ptr = plus_delta_x->row_ptr;
    col_ind = plus_delta_x->col_ind;
    values = plus_delta_x->values;

    idx = 0;
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
            values[idx] = 1.0 - mu_x;
            idx++;

            // Left
            col = grid->id_map[gi - 1][gj];
            col_ind[idx] = col;
            values[idx] = mu_x / 2;
            idx++;

            // Right
            col = grid->id_map[gi + 1][gj];
            col_ind[idx] = col;
            values[idx] = mu_x / 2;
            idx++;
        }
        // Boundary point
        else {
            int row = i;
            int col = i;
            col_ind[idx] = col;
            values[idx] = 0.0;
            idx++;
        }
        row_ptr[i + 1] = idx;
    }
    plus_delta_x->nnz = idx; // Update nnz
    ADI_matrixs[2] = plus_delta_x;

    SparseCSR *minus_delta_y = createSparseCSR(grid->n_active, grid->n_active, 3 * grid->n_active);
    row_ptr = minus_delta_y->row_ptr;
    col_ind = minus_delta_y->col_ind;
    values = minus_delta_y->values;

    idx = 0;
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
            values[idx] = 1.0 + mu_y;
            idx++;

            // Down
            col = grid->id_map[gi][gj - 1];
            col_ind[idx] = col;
            values[idx] = - mu_y / 2;
            idx++;

            // Up
            col = grid->id_map[gi][gj + 1];
            col_ind[idx] = col;
            values[idx] = - mu_y / 2;
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
    minus_delta_y->nnz = idx; // Update nnz
    ADI_matrixs[3] = minus_delta_y;

    return ADI_matrixs;
}

void* solve_Parabolic_ADI(Grid2D* grid, SparseCSR** ADI_matrix, double t, double tau, double* b, double) {
    SparseCSR *plus_delta_y = ADI_matrix[0], *minus_delta_x = ADI_matrix[1],
              *plus_delta_x = ADI_matrix[2], *minus_delta_y = ADI_matrix[3];
    
}

// double* assemble_RHS_Neumann(Grid2D* grid, f_func f, boundary_func compute_boundary_value, f_func get_exact) {
//     double *b = (double *)malloc(grid->n_active * sizeof(double));
//     double h = grid->hx; // Assuming hx = hy
//     int flag = 1;
//     for (int i = 0; i < grid->n_active; i++) {
//         int gi = grid->id_i[i];
//         int gj = grid->id_j[i];
//         double xi = grid->x[gi];
//         double yj = grid->y[gj];

//         // Interior point
//         if (grid->region[gi][gj] == 1) {
//             if (flag) {
//                 b[i] = get_exact(xi, yj);
//                 flag = 0;
//             } else {
//                 b[i] = f(xi, yj) * h * h; // Source term
//             }
//         }
//         // Boundary point
//         else {
//             // Neumann conditions based on boundary type
//             b[i] = compute_boundary_value(xi, yj, grid->region[gi][gj]) * h;
//         }
//     }
//     return b;
// }