# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <csr.h>
# include <poisson2d.h>
# include <utils.h>
# define Pi 3.14159265358979323846

int region_divider(double x, double y, double hx, double hy) {
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
    double x_b, y_b, result;
    switch (boundary_type) {
        case 2: // Top left slant boundary
            x_b = (x + y - 1.0) / 2.0;
            y_b = (x + y + 1.0) / 2.0;
            result = compute_u_exact(x_b, y_b);
            break;
        case 3: // Top right slant boundary
            x_b = (x - y + 3.0) / 2.0;
            y_b = (-x + y + 3.0) / 2.0;
            result = compute_u_exact(x_b, y_b);
            break;
        case 4: // Left boundary
            x_b = 0.0;
            y_b = y;
            result = compute_u_exact(x_b, y_b);
            break;
        case 5: // Upper right boundary
            x_b = (x + 2.0 * y + 6.0) / 5.0;
            y_b = (2.0 * x + 4.0 * y -3.0) / 5.0;
            result = compute_u_exact(x_b, y_b);
            break;
        case 6: // Lower right boundary
            x_b = (x - y) / 2.0;
            y_b = (-x + y) / 2.0;
            result = compute_u_exact(x_b, y_b);
            break;
        case 7: // Bottom boundary
            x_b = x;
            y_b = -2.0;
            result = compute_u_exact(x_b, y_b);
            // printf("Boundary point (%.4f, %.4f) with value %.6f\n", x_b, y_b, result);
            break;
        default:
            result = 0.0;
    }
    return result;
}

int main() {
    int nx = 41;
    int ny = 81;
    Grid2D* grid = initialize_Grid(nx, ny, 0.0, 2.0, -2.0, 2.0, region_divider);
    printf("Number of active grid points: %d\n", grid->n_active);
    // printf("%.6f\n", compute_u_exact(grid->x[10], grid->y[40]));
    print_int_matrix((const int **)grid->region, grid->nx, grid->ny);
    SparseCSR* matrix = assemble_Matrix_Dirichlet(grid);

    // printf("Assembled sparse matrix in CSR format:\n");
    // if (grid->n_active <= 36) {
    //     print_SparseCSR(matrix, 1);
    // } else {
    //     printf("Matrix too large to display fully.\n");
    //     print_SparseCSR_simple(matrix, 1);
    // }

    double* rhs = assemble_RHS_Dirichlet(grid, compute_f, compute_boundary_value);
    // Solve the linear system using Conjugate Gradient method
    double* solution = (double *)malloc(grid->n_active * sizeof(double));
    for (int i = 0; i < grid->n_active; i++) {
        solution[i] = 0.0; // Initial guess
    }
    GaussSeidel_csr(matrix, rhs, solution, 1000, 1e-6);

    double *exact = (double *)malloc(grid->n_active * sizeof(double));
    for (int i = 0; i < grid->n_active; i++) {
        int gi = grid->id_i[i];
        int gj = grid->id_j[i];
        double xi = grid->x[gi];
        double yj = grid->y[gj];
        exact[i] = compute_u_exact(xi, yj);
        // if (fabs(yj + 2.0) < 1e-12) {
        //     printf("Point (%.4f, %.4f): Exact = %.6f\n", xi, yj, exact[i]);
        //     printf("Computed = %.6f\n", solution[i]);
        // }
    }
    
    // Output the solution at active grid points
    // printf("Computed solution at active grid points:\n");
    // for (int i = 0; i < grid->n_active; i++) {
    //     int gi = grid->id_i[i];
    //     int gj = grid->id_j[i];
    //     double xi = grid->x[gi];
    //     double yj = grid->y[gj];
    //     double u_exact = compute_u_exact(xi, yj);
    //     printf("u(%.4f, %.4f) = %.6f, Exact = %.6f, Error = %.6e\n", xi, yj, solution[i], u_exact, fabs(solution[i] - u_exact));
    // }

    double **data_points = read_indices_to_points(grid, solution);
    // print_matrix((const double **)data_points, grid->nx, grid->ny, 6);
    double **exact_points = read_indices_to_points(grid, exact);
    // print_matrix((const double **)exact_points, grid->nx, grid->ny, 6);

    // Optionally, write the solution to a CSV file for visualization
    write_csv_matrix("results/Dirichlet_solution.csv", data_points, grid->nx, grid->ny);
    write_csv_matrix("results/Dirichlet_exact.csv", exact_points, grid->nx, grid->ny);

    // Clean up
    // Save grid dimensions because free_grid(grid) will deallocate the structure
    int nx_grid = grid->nx;
    int ny_grid = grid->ny;

    freeSparseCSR(matrix);
    free(rhs);
    free(solution);

    for (int i = 0; i < nx_grid; i++) {
        free(data_points[i]);
        free(exact_points[i]);
    }
    free(data_points);
    free(exact_points);

    // Now free the grid structure itself
    free_grid(grid);
    return 0;
}