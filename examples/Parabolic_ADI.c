/**
 * @file Parabolic_ADI.c
 * @brief Example: solve the 2D parabolic PDE with Dirichlet boundary conditions using ADI.
 *
 * @details
 * This example demonstrates constructing the mesh, assembling ADI split
 * operators and time-stepping a manufactured solution using an alternating
 * direction implicit (ADI) method. Output CSVs are produced for visualization
 * and verification; see `ReadMe.md` for usage notes and expected outputs.
 *
 * @see csr.h, parabolic.h, bessel.h
 * @author Li Zhijun
 * @date 2025-12-03
 * @example Parabolic_ADI.c
 */
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <vec.h>
# include <utils.h>
# include <bessel.h>
# include <parabolic.h>

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

double compute_u_exact(double x, double y, double t, double hx, double hy) {
    double r = sqrt((x - 1) * (x - 1) + (y - 1) * (y - 1));
    if (r > (sqrt(hx * hx + hy * hy) / 2)) {
        return -plane_solution_function(r, t) / 4;
    }
    else {
        return -average_cell(hx, hy, t) / 4;
    }
}

double compute_u_boundary(double x, double y, double t) {
    double r = sqrt((x - 1) * (x - 1) + (y - 1) * (y - 1));
    return -plane_solution_function(r, t) / 4;
}

double integrated_source_term(double x, double y, double t, double hx, double hy) {
    if ((fabs(x - 1) < (hx / 2)) && (fabs(y - 1) < (hy / 2))) {
        return sin(t) / hx / hy;
    }
    else {
        return 0;
    }
}

double compute_boundary_value(double x, double y, double t, int boundary_type) {
    double x_b, y_b, result;
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
    result = compute_u_boundary(x_b, y_b, t);
    return result;
}

int main(){
    double T_max = 2 * M_PI;
    int nx = 41;
    int ny = 81;
    Grid2D* grid = initialize_Grid(nx, ny, 0.0, 2.0, -2.0, 2.0, region_divider);
    double tau = grid->hx * grid->hx * grid->hy * grid->hy / (grid->hx * grid->hx + grid->hy * grid->hy) / 2 * 5;
    SparseCSR** ADI_matrixs = assemble_Matrix_Parabolic_ADI(grid, tau);
    SparseCSR *plus_delta_y = ADI_matrixs[0], *minus_delta_x = ADI_matrixs[1],
              *plus_delta_x = ADI_matrixs[2], *minus_delta_y = ADI_matrixs[3];
    
    // printf("Number of active grid points: %d\n", grid->n_active);
    // printf("%.6f\n", compute_u_exact(grid->x[10], grid->y[40]));
    // printf("Grid region layout (0: exterior, 1: interior, others: boundary types):\n");
    // print_int_matrix((const int **)grid->region, grid->nx, grid->ny);

    double t_now = 0.0;
    int step = 0;
    int output_interval = 20;

    double *exact = (double *)malloc(grid->n_active * sizeof(double));
    double *solution = (double *)malloc(grid->n_active * sizeof(double));
    double *rhs = (double *)malloc(grid->n_active * sizeof(double));
    double *temp = (double *)malloc(grid->n_active * sizeof(double));
    // double *u_star = (double *)malloc(grid->n_active * sizeof(double));
    
    double **exact_points = create_grid_2D_array(grid);
    double **solution_points = create_grid_2D_array(grid);

    for (int i = 0; i < grid->n_active; i++) {
        int gi = grid->id_i[i];
        int gj = grid->id_j[i];
        double xi = grid->x[gi];
        double yj = grid->y[gj];
        solution[i] = compute_u_exact(xi, yj, t_now, grid->hx, grid->hy);
    }

    write_csv_int_matrix("results/Parabolic/data/ADI/grid_data.csv", grid->region, grid->nx, grid->ny);

    while (t_now < T_max) {
        t_now += tau;
        step ++;
        
        spmv_csr(plus_delta_y, solution, temp);
        assemble_RHS_Parabolic(grid, integrated_source_term, compute_boundary_value, rhs, t_now - tau / 2, tau / 2);
        vec_add(rhs, temp, grid->n_active);
        GaussSeidel_csr(minus_delta_x, rhs, solution, 20, 1e-6);

        spmv_csr(plus_delta_x, solution, temp);
        assemble_RHS_Parabolic(grid, integrated_source_term, compute_boundary_value, rhs, t_now, tau / 2);
        vec_add(rhs, temp, grid->n_active);
        GaussSeidel_csr(minus_delta_y, rhs, solution, 20, 1e-6);

        for (int i = 0; i < grid->n_active; i++) {
            int gi = grid->id_i[i];
            int gj = grid->id_j[i];
            double xi = grid->x[gi];
            double yj = grid->y[gj];
            exact[i] = compute_u_exact(xi, yj, t_now, grid->hx, grid->hy);
        }
        if ((step % output_interval) == 0) {
            printf("Current Step: %06d, Writing Output\n", step);
            read_indices_to_points(grid, exact, exact_points);
            read_indices_to_points(grid, solution, solution_points);
            
            char fname_exact[256];
            sprintf(fname_exact, "results/Parabolic/data/ADI/exact_%06d.csv", step);
            char fname_rhs[256];
            sprintf(fname_rhs, "results/Parabolic/data/ADI/solution_%06d.csv", step);
            write_csv_matrix(fname_exact, exact_points, grid->nx, grid->ny);
            write_csv_matrix(fname_rhs, solution_points, grid->nx, grid->ny);
        }
    }

    free(exact);
    free(solution);
    free(rhs);
    free(temp);
    free_grid_2D_array(exact_points, grid);
    free_grid_2D_array(solution_points, grid);
    free_grid(grid);
    return 0;
}