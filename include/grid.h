/**
 * @file grid.h
 * @brief Header file for a uniform grid structure and related functions.
 * 
 * This header file declares the Grid structure and functions to create and free a grid.
 * The Grid structure represents a 2D grid with specified dimensions and spacing.
 * @see grid.c
 * @author Li Zhijun
 * @date 2025-10-10
 */
# ifndef GRID_H
# define GRID_H

typedef struct {
    int nx;         /**< Number of grid points in the x direction */
    int ny;         /**< Number of grid points in the y direction */
    double hx;      /**< Grid spacing in the x direction */
    double hy;      /**< Grid spacing in the y direction */
    double x0, x1;  /**< Domain boundaries in the x direction */
    double y0, y1;  /**< Domain boundaries in the y direction */

    double *x;      /**< Array of x-coordinates of grid points */
    double *y;      /**< Array of y-coordinates of grid points */
    int **region;   /**< Region array indicating active grid points and boundaries */

    int **id_map;   /**< Mapping from 2D grid points to 1D indices */
    int n_active;   /**< Number of active grid points */

    int *id_i;      /**< Array of i-coordinates of active points */
    int *id_j;      /**< Array of j-coordinates of active points */
} Grid2D;

Grid2D* create_uniform_grid(int nx, int ny, double x0, double x1, double y0, double y1);
void* free_grid(Grid2D *grid);

# endif