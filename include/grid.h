/**
 * @file grid.h
 * @brief Header file for a uniform grid structure and related functions.
 * 
 * This header file declares the Grid structure and functions to create and free a grid.
 * The Grid structure represents a 2D grid with specified dimensions and spacing.
 * This file also declares functions for creating, initializing, and data remapping.
 * @see grid.c
 * @author Li Zhijun
 * @date 2025-10-21
 */
# ifndef GRID_H
# define GRID_H

typedef int (*region_divider_func)(double, double, double, double);

/**
 * @struct Grid2D
 * @brief Structure to represent a 2D orthogonal grid.
 * 
 * 
 */
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
    int n_interior; /**< Number of interior grid points */

    int *id_i;      /**< Array of i-coordinates of active points */
    int *id_j;      /**< Array of j-coordinates of active points */
} Grid2D;

/**
 * @brief Create new uniform 2D grid points.
 * @param nx Number of points in x direction.
 * @param ny Number of points in y direction.
 * @param x0 X-coordinate on the left of the domain.
 * @param x1 X-coordinate on the right of the domain.
 * @param y0 Y-coordinate on the bottom of the domain.
 * @param y1 Y-coordinate on the top of the domain.
 * 
 * @note This function is a built-in method and it is not recommended to call it directly from the outside.
 * @note Use initialize_Grid to create a new grid structure
 * @see initialize_Grid()
 */
Grid2D* create_uniform_grid(int nx, int ny, double x0, double x1, double y0, double y1);

/**
 * @brief Create a new 2D grid structure with the division of the computing domain.
 * @param nx Number of points in x direction.
 * @param ny Number of points in y direction.
 * @param x0 X-coordinate on the left of the domain.
 * @param x1 X-coordinate on the right of the domain.
 * @param y0 Y-coordinate on the bottom of the domain.
 * @param y1 Y-coordinate on the top of the domain.
 * @param region_divider Function pointer to define which domain a grid point belongs to.
 * @par Function signature requirements:
 * @code
 * int region_divider(double x, double y, double hx, double hy)
 * @endcode
 *      - param x: X-coordinate of the grid point
 *      - param y: Y-coordinate of the grid point
 *      - param hx: Grid spacing in the x direction
 *      - param hy: Grid spacing in the y direction
 *      - return value: An integer representing the area to which the grid point belongs.
 *      - return value: ==0 -> not calculated; ==1 -> interior point; >1 -> boundarys
 * 
 * @note The caller is responsible for freeing the allocated memory using free_grid().
 * @see free_grid()
 */
Grid2D* initialize_Grid(int nx, int ny, double x0, double x1, double y0, double y1, region_divider_func region_divider);

double **create_grid_2D_array(Grid2D *grid);

void* free_grid_2D_array(double** array, Grid2D *grid);

/**
 * @brief Remap the data in the form of column vectors to the grid points.
 * @param grid Pointer to the grid structure with mapping relationships established by initalize_Grid().
 * @param data_indices The data in the form of column vectors.
 * 
 * @see initialize_Grid()
 */
void read_indices_to_points(Grid2D *grid, double* data_indices, double** data_points);

/**
 * @brief Free the memory allocated for a grid structure.
 * @param grid Pointer to the grid structure to free.
 */
void* free_grid(Grid2D *grid);

# endif