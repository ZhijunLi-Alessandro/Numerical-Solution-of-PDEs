# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <csr.h>
# include <grid.h>
# include <utils.h>
# define Pi 3.14159265358979323846

int is_activate(double x, double y, double hx, double hy) {
    double eps = 1e-12;
    if (y > 1.0 && y <= (2.0 + eps)) {
        if (x >= (y - 1.0 - eps) && x <= (3.0 - y + eps)) {
            if (x <= y - 1.0 + hx - 2 * eps) {
                return 2; // Top left slant boundary
            } else if (x >= 3.0 - y - hx + 2 * eps) {
                return 3; // Top right slant boundary
            } else {
                return 1; // Active interior point
            }
        } else {
            return 0;
        }
    }
    else if (y > -1.0 && y <= 1.0) {
        if (x >= -eps && x <= (0.5 * y + 1.5 + eps)) {
            if (x <= hx - 2 *eps) {
                return 4; // Left boundary
            } else if (x >= 0.5 * y + 1.5 - hx + 2 * eps) {
                return 5; // Upper right boundary
            } else {
                return 1; // Active interior point
            }
        } else {
            return 0;
        }
    }
    else if (y >= (-2.0 - eps) && y <= -1.0) {
        if (x >= -eps && x <= (-y + eps)) {
            if (x <= hx - 2 * eps) {
                return 4; // Left boundary
            } else if (x >= -y - hx + 2 * eps) {
                return 6; // Lower right boundary
            } else if (y <= -2.0 + hy - 2 * eps) {
                return 7; // Bottom boundary
            } else {
                return 1; // Active interior point
            }
        } else {
            return 0;
        }
    }
};

double compute_f(double x, double y) {
    return sin(Pi * x) * cos(2 * Pi * y);
}

double compute_u_exact(double x, double y) {
    return 1.0 / (5 * Pi * Pi) * sin(Pi * x) * cos(2 * Pi * y);
}

double compute_boundary_value(double x, double y, int boundary_type) {
    double x_b, y_b;
    switch (boundary_type) {
        case 2: // Top left slant boundary
            x_b = (x + y - 1.0) / 2.0;
            y_b = (x + y + 1.0) / 2.0;
        case 3: // Top right slant boundary
            x_b = (x - y + 3.0) / 2.0;
            y_b = (-x + y + 3.0) / 2.0;
        case 4: // Left boundary
            x_b = 0.0;
            y_b = y;
        case 5: // Upper right boundary
            x_b = (x + 2.0 * y + 6.0) / 5.0;
            y_b = (2.0 * x + 4.0 * y -3.0) / 5.0;
        case 6: // Lower right boundary
            x_b = (x - y) / 2.0;
            y_b = (-x + y) / 2.0;
        case 7: // Bottom boundary
            x_b = x;
            y_b = -2.0;
        default:
            return 0.0;
        return compute_u_exact(x_b, y_b);
    }
}

Grid2D* initial_Grid(int nx, int ny) {
    double x0 = 0.0, x1 = 2.0, y0 = -2.0, y1 = 2.0;
    Grid2D *grid = create_uniform_grid(nx, ny, x0, x1, y0, y1);
    grid->n_active = 0;
    for (int i = 0; i < grid->nx; i++) {
        for (int j = 0; j < grid->ny; j++) {
            int region_value = is_activate(grid->x[i], grid->y[j], grid->hx, grid->hy);
            if (region_value > 0) {
                grid->region[i][j] = region_value;
                grid->id_map[i][j] = grid->n_active;
                grid->n_active++;
            } else {
                grid->region[i][j] = 0;
                grid->id_map[i][j] = -1;
            }
        }
    }
    grid->id_i = (int *)malloc(grid->n_active * sizeof(int));
    grid->id_j = (int *)malloc(grid->n_active * sizeof(int));
    for (int i = 0; i < grid->nx; i++) {
        for (int j = 0; j < grid->ny; j++) {
            if (grid->region[i][j] == 1) {
                grid->id_i[grid->id_map[i][j]] = i;
                grid->id_j[grid->id_map[i][j]] = j;
            }
        }
    }
    return grid;
}

SparseCSR* assemble_Matrix(Grid2D* grid) {
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

double* assemble_RHS(Grid2D* grid) {
    double *b = (double *)malloc(grid->n_active * sizeof(double));
    double h = grid->hx; // Assuming hx = hy
    for (int i = 0; i < grid->n_active; i++) {
        int gi = grid->id_i[i];
        int gj = grid->id_j[i];
        double xi = grid->x[gi];
        double yj = grid->y[gj];

        // Interior point
        if (grid->region[gi][gj] == 1) {
            b[i] = compute_f(xi, yj) * h * h; // Source term
        }
        // Boundary point
        else {
            // Dirichlet conditions based on boundary type
            b[i] = compute_boundary_value(xi, yj, grid->region[gi][gj]);
        }
    }
    return b;
}

int main() {
    int nx = 41;
    int ny = 81;
    Grid2D* grid = initial_Grid(nx, ny);
    printf("Number of active grid points: %d\n", grid->n_active);
    print_int_matrix((const int **)grid->region, grid->nx, grid->ny);
    SparseCSR* matrix = assemble_Matrix(grid);

    // printf("Assembled sparse matrix in CSR format:\n");
    // if (grid->n_active <= 36) {
    //     print_SparseCSR(matrix, 1);
    // } else {
    //     printf("Matrix too large to display fully.\n");
    //     print_SparseCSR_simple(matrix, 1);
    // }

    double* rhs = assemble_RHS(grid);
    // Solve the linear system using Conjugate Gradient method
    double* solution = (double *)malloc(grid->n_active * sizeof(double));
    for (int i = 0; i < grid->n_active; i++) {
        solution[i] = 0.0; // Initial guess
    }
    CG_csr(matrix, rhs, solution, 1000, 1e-6);
    
    // Output the solution at active grid points
    printf("Computed solution at active grid points:\n");
    for (int i = 0; i < grid->n_active; i++) {
        int gi = grid->id_i[i];
        int gj = grid->id_j[i];
        double xi = grid->x[gi];
        double yj = grid->y[gj];
        double u_exact = compute_u_exact(xi, yj);
        printf("u(%.4f, %.4f) = %.6f, Exact = %.6f, Error = %.6e\n", xi, yj, solution[i], u_exact, fabs(solution[i] - u_exact));
    }
    // Clean up
    free_grid(grid);
    freeSparseCSR(matrix);
    free(rhs);
    free(solution);
    return 0;
}