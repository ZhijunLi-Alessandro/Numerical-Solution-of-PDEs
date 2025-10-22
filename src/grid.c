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