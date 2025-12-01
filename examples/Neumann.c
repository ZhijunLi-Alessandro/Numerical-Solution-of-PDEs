/**
 * @file Neumann.c
 * @brief Example program to solve a 2D Poisson problem with Neumann boundary condition.
 * 
 * @details
 * This program define a Poisson problem of Neumann boundary conditions on the two-
 * dimensional irregular computational region of periodic solutions and constructs a 
 * finite difference method to numerically solve and calculate it. 
 * The solution is set to be periodic. (u=(1/(5*Pi*Pi))*sin*(Pi*x)*cos(2*Pi*x)).
 * See details in ReadMe.md
 * 
 * @see csr.h, poisson2d.h, grid.h
 * @author Li Zhijun
 * @date 2025-10-28
 * @example Neumann.c
 */
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
                // if (y <= -2.0 + hy - 2 * eps) {
                //     return 0; // Exclude Corner
                // }
                return 4; // Left boundary
            } else if (x >= -y - hx + 2 * eps) {
                if (y <= -2.0 + hy - 2 * eps) {
                    return 0; // Exclude Corner
                }
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

double get_normal(int boundary_type) {
    switch (boundary_type) {
        case 2: // Top left slant boundary
            return 3.0 * Pi / 4.0;
        case 3: // Top right slant boundary
            return Pi / 4.0;
        case 4: // Left boundary
            return Pi;
        case 5: // Upper right boundary
            return - atan(0.5);
        case 6: // Lower right boundary
            return Pi / 4.0;
        case 7: // Bottom boundary
            return -Pi / 2.0;
        default:
            return 0.0; // Not a boundary
    }
}

double compute_f(double x, double y) {
    return sin(Pi * x) * cos(2 * Pi * y);
}

double compute_u_exact(double x, double y) {
    return 1.0 / (5 * Pi * Pi) * sin(Pi * x) * cos(2 * Pi * y);
}

double compute_derivative_exact(double x, double y, int direction) {
    if (direction == 0) { // derivative w.r.t x
        return (1.0 / (5 * Pi * Pi)) * Pi * cos(Pi * x) * cos(2 * Pi * y);
    } else if (direction == 1) { // derivative w.r.t y
        return (1.0 / (5 * Pi * Pi)) * (-2.0 * Pi) * sin(Pi * x) * sin(2 * Pi * y);
    } else {
        return 0.0;
    }
}

double compute_boundary_value(double x, double y, int boundary_type) {
    double x_b, y_b, result, derivx, derivy, alpha;
    alpha = get_normal(boundary_type);
    switch (boundary_type) {
        case 2: // Top left slant boundary
            x_b = (x + y - 1.0) / 2.0;
            y_b = (x + y + 1.0) / 2.0;
            break;
        case 3: // Top right slant boundary
            x_b = (x - y + 3.0) / 2.0;
            y_b = (-x + y + 3.0) / 2.0;
            break;
        case 4: // Left boundary
            x_b = 0.0;
            y_b = y;
            break;
        case 5: // Upper right boundary
            x_b = (x + 2.0 * y + 6.0) / 5.0;
            y_b = (2.0 * x + 4.0 * y -3.0) / 5.0;
            break;
        case 6: // Lower right boundary
            x_b = (x - y) / 2.0;
            y_b = (-x + y) / 2.0;
            break;
        case 7: // Bottom boundary
            x_b = x;
            y_b = -2.0;
            break;
        default:
            return 0.0;
    }
    derivx = compute_derivative_exact(x_b, y_b, 0);
    derivy = compute_derivative_exact(x_b, y_b, 1);
    result = (derivx * cos(alpha) + derivy * sin(alpha));
    return result;
}

double *numerical_deriv(Grid2D* grid, double*data_indices, int direction) {
    double *derivative = (double *)malloc(grid->n_active * sizeof(double));
    int **region = grid->region, **id_map = grid->id_map;
    double hx = grid->hx, hy = grid->hy;
    for (int i = 0; i < grid->n_active; i++) {
        int gi = grid->id_i[i];
        int gj = grid->id_j[i];
        if (direction == 0) { // numerical derivative w.r.t x
            if ((gi == 0) || (region[gi - 1][gj] == 0)) {
                derivative[i] = (data_indices[id_map[gi + 1][gj]] - data_indices[i]) / hx;
            }
            else if ((gi == grid->nx-1) || (region[gi + 1][gj] == 0)) {
                derivative[i] = (data_indices[i] - data_indices[id_map[gi - 1][gj]]) / hx;
            }
            else {
                derivative[i] = (data_indices[id_map[gi + 1][gj]] - data_indices[id_map[gi - 1][gj]]) / hx / 2;
            }
        }
        else if (direction == 1) { // numerical derivative w.r.t y
            if ((gj == 0) || (region[gi][gj - 1] == 0)) {
                derivative[i] = (data_indices[id_map[gi][gj + 1]] - data_indices[i]) / hy;
            }
            else if ((gj == grid->ny-1) || (region[gi][gj + 1] == 0)) {
                derivative[i] = (data_indices[i] - data_indices[id_map[gi][gj - 1]]) / hy;
            }
            else {
                derivative[i] = (data_indices[id_map[gi][gj + 1]] - data_indices[id_map[gi][gj - 1]]) / hy / 2;
            }
        }
    }
    return derivative;
}

int main() {
    int nx = 41;
    int ny = 81;
    Grid2D* grid = initialize_Grid(nx, ny, 0.0, 2.0, -2.0, 2.0, region_divider);
    // printf("Number of active grid points: %d\n", grid->n_active);
    // printf("%.6f\n", compute_u_exact(grid->x[10], grid->y[40]));
    printf("Grid region layout (0: exterior, 1: interior, others: boundary types):\n");
    print_int_matrix((const int **)grid->region, grid->nx, grid->ny);
    SparseCSR* matrix = assemble_Matrix_Neumann(grid, get_normal);

    // printf("Assembled sparse matrix in CSR format:\n");
    // if (grid->n_active <= 36) {
    //     print_SparseCSR(matrix, 1);
    // } else {
    //     printf("Matrix too large to display fully.\n");
    //     print_SparseCSR_simple(matrix, 1);
    // }

    double* rhs = assemble_RHS_Neumann(grid, compute_f, compute_boundary_value, compute_u_exact);
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

    double *solution_deriv = numerical_deriv(grid, solution, 0);
    double *exact_deriv = numerical_deriv(grid, exact, 0);

    double **data_points = create_grid_2D_array(grid);
    read_indices_to_points(grid, solution, data_points);
    // double **data_points = read_indices_to_points(grid, solution_deriv);
    // print_matrix((const double **)data_points, grid->nx, grid->ny, 6);
    double **exact_points = create_grid_2D_array(grid);
    read_indices_to_points(grid, exact, exact_points);
    // double **exact_points = read_indices_to_points(grid, exact_deriv);
    // print_matrix((const double **)exact_points, grid->nx, grid->ny, 6);

    // Optionally, write the solution to a CSV file for visualization
    write_csv_matrix("results/Poisson/data/Neumann_solution.csv", data_points, grid->nx, grid->ny);
    write_csv_matrix("results/Poisson/data/Neumann_exact.csv", exact_points, grid->nx, grid->ny);
    write_csv_int_matrix("results/Poisson/data/grid_data.csv", grid->region, grid->nx, grid->ny);

    // Clean up
    // Save grid dimensions because free_grid(grid) will deallocate the structure
    int nx_grid = grid->nx;
    int ny_grid = grid->ny;

    freeSparseCSR(matrix);
    free(rhs);
    free(exact);
    free(solution);

    free_grid_2D_array(data_points, grid);
    free_grid_2D_array(exact_points, grid);

    // Now free the grid structure itself
    free_grid(grid);
    return 0;
}