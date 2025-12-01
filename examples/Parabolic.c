# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <csr.h>
# include <poisson2d.h>
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

int main(){
    double T_max = 6 * M_PI;
    int nx = 41;
    int ny = 81;
    Grid2D* grid = initialize_Grid(nx, ny, 0.0, 2.0, -2.0, 2.0, region_divider);
    double t = grid->hx * grid->hx * grid->hy * grid->hy / (grid->hx * grid->hx + grid->hy * grid->hy) / 4 * 12;
    SparseCSR* iteration_matrix = assemble_Matrix_Parabolic_Explicit(grid, t);
    // printf("Number of active grid points: %d\n", grid->n_active);
    // printf("%.6f\n", compute_u_exact(grid->x[10], grid->y[40]));
    // printf("Grid region layout (0: exterior, 1: interior, others: boundary types):\n");
    // print_int_matrix((const int **)grid->region, grid->nx, grid->ny);
    double t_now = 0.0;
    int step = 0;
    int output_interval = 100;

    double *exact = (double *)malloc(grid->n_active * sizeof(double));
    double *compute = (double *)malloc(grid->n_active * sizeof(double));
    double *rhs = (double *)malloc(grid->n_active * sizeof(double));
    
    double **exact_points = create_grid_2D_array(grid);
    double **rhs_points = create_grid_2D_array(grid);

    for (int i = 0; i < grid->n_active; i++) {
        int gi = grid->id_i[i];
        int gj = grid->id_j[i];
        double xi = grid->x[gi];
        double yj = grid->y[gj];
        exact[i] = compute_u_exact(xi, yj, t_now, grid->hx, grid->hy);
    }

    write_csv_int_matrix("results/Parabolic/data/grid_data.csv", grid->region, grid->nx, grid->ny);

    while (t_now < T_max) {
        t_now += t;
        step ++;
        spmv_csr(iteration_matrix, exact, compute);
        for (int i = 0; i < grid->n_active; i++) {
            int gi = grid->id_i[i];
            int gj = grid->id_j[i];
            double xi = grid->x[gi];
            double yj = grid->y[gj];
            exact[i] = compute_u_exact(xi, yj, t_now, grid->hx, grid->hy);
        }
        for (int i = 0; i < grid->n_active; i++) {
            int gi = grid->id_i[i];
            int gj = grid->id_j[i];
            if (grid->region[gi][gj] == 1) {
                rhs[i] = exact[i] - compute[i];
            }
            else {
                rhs[i] = 0;
            }
        }
        if ((step % output_interval) == 0) {
            printf("Current Step: %06d, Writing Output\n", step);
            read_indices_to_points(grid, exact, exact_points);
            read_indices_to_points(grid, rhs, rhs_points);
            
            char fname_exact[256];
            sprintf(fname_exact, "results/Parabolic/data/exact_%06d.csv", step);
            char fname_rhs[256];
            sprintf(fname_rhs, "results/Parabolic/data/rhs_%06d.csv", step);
            write_csv_matrix(fname_exact, exact_points, grid->nx, grid->ny);
            write_csv_matrix(fname_rhs, rhs_points, grid->nx, grid->ny);
        }
    }

    free(exact);
    free(compute);
    free(rhs);
    free_grid_2D_array(exact_points, grid);
    free_grid_2D_array(rhs_points, grid);
    free_grid(grid);
    return 0;
}