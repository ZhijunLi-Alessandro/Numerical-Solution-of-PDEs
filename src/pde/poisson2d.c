#include <stdlib.h>
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
            // b[i] = 0.0;
        }
    }
    return b;
}