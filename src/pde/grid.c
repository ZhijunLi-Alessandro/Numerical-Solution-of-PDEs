/**
 * @file grid.c
 * @brief Implementation of uniform grid operations.
 * 
 * This file includes functions for creating, freeing, initializing, and
 * remapping data in form of column vectors to data in form of matrix.
 * 
 * @author Li Zhijun
 * @date 2025-10-21
 */
#include <stdlib.h>
#include "grid.h"

Grid2D* create_uniform_grid(int nx, int ny, double x0, double x1, double y0, double y1) {
    Grid2D *grid = (Grid2D *)malloc(sizeof(Grid2D));
    grid->nx = nx;
    grid->ny = ny;
    grid->x0 = x0;
    grid->x1 = x1;
    grid->y0 = y0;
    grid->y1 = y1;
    grid->hx = (x1 - x0) / (nx - 1);
    grid->hy = (y1 - y0) / (ny - 1);

    grid->x = (double *)malloc(nx * sizeof(double));
    grid->y = (double *)malloc(ny * sizeof(double));

    for (int i = 0; i < nx; i++) {
        grid->x[i] = x0 + i * grid->hx;
    }
    for (int j = 0; j < ny; j++) {
        grid->y[j] = y0 + j * grid->hy;
    }

    grid->region = (int **)malloc(nx * sizeof(int *));
    for (int i = 0; i < nx; i++) {
        grid->region[i] = (int *)malloc(ny * sizeof(int));
        for (int j = 0; j < ny; j++) {
            grid->region[i][j] = 0; // Default: all points are exterior
        }
    }
    grid->id_map = (int **)malloc(nx * sizeof(int *));
    for (int i = 0; i < nx; i++) {
        grid->id_map[i] = (int *)malloc(ny * sizeof(int));
        for (int j = 0; j < ny; j++) {
            grid->id_map[i][j] = -1; // Initialize to -1
        }
    }

    return grid;
}

Grid2D* initialize_Grid(int nx, int ny, double x0, double x1, double y0, double y1, region_divider_func region_divider) {
    Grid2D *grid = create_uniform_grid(nx, ny, x0, x1, y0, y1);
    grid->n_active = 0;
    grid->n_interior = 0;
    for (int i = 0; i < grid->nx; i++) {
        for (int j = 0; j < grid->ny; j++) {
            int region_value = region_divider(grid->x[i], grid->y[j], grid->hx, grid->hy);
            if (region_value > 0) {
                grid->region[i][j] = region_value;
                grid->id_map[i][j] = grid->n_active;
                grid->n_active++;
                if (region_value == 1) {
                    grid->n_interior++;
                }
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
            if (grid->region[i][j] > 0) {
                grid->id_i[grid->id_map[i][j]] = i;
                grid->id_j[grid->id_map[i][j]] = j;
            }
        }
    }
    return grid;
}

double **read_indices_to_points(Grid2D *grid, double* data_indices) {
    double **data_points = (double **)malloc(grid->nx * sizeof(double *));
    for (int i = 0; i < grid->nx; i++) {
        data_points[i] = (double *)malloc(grid->ny * sizeof(double));
        for (int j = 0; j < grid->ny; j++) {
            if (grid->region[i][j] == 0) {
                data_points[i][j] = 0.0; // or some sentinel value for inactive points
            } else {
                data_points[i][j] = data_indices[grid->id_map[i][j]];
            }
        }
    }
    return data_points;
}

void* free_grid(Grid2D *grid) {
    if (grid) {
        free(grid->x);
        free(grid->y);
        if (grid->region) {
            for (int i = 0; i < grid->nx; i++) {
                free(grid->region[i]);
            }
            free(grid->region);
        }
        if (grid->id_map) {
            for (int i = 0; i < grid->nx; i++) {
                free(grid->id_map[i]);
            }
            free(grid->id_map);
        }
        if (grid->id_i) {
            free(grid->id_i);
        }
        if (grid->id_j) {
            free(grid->id_j);
        }
        free(grid);
    }
}