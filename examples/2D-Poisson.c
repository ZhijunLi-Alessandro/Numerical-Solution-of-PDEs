/**
 * @file 2D-Poisson.c
 * @brief Example program to generate a 2D Poisson matrix in CSR format.
 * 
 * @details
 * This program generates a sparse matrix representing the 2D Poisson equation
 * using finite difference discretization on a grid of size (nx+1) x (ny+1).
 * The matrix is stored in Compressed Sparse Row (CSR) format. The program also
 * demonstrates the use of the get_diag_L_U_csr function to extract the diagonal,
 * lower, and upper triangular parts of the matrix.
 * 
 * Usage:
 * Compile the program and run it with two command-line arguments specifying
 * the number of interior grid points in the x and y directions (nx and ny).
 * 
 * Example:
 * \verbatim
   mkdir build && cd build
   cmake ..
   make
   ../bin/2D-Poisson_2d 5 5 \endverbatim
 * 
 * @see csr.h, utils.h
 * @author Li Zhijun
 * @date 2025-10-10
 * @example 2D-Poisson.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <csr.h>
#include <utils.h>

/**
 * @brief Map 2D grid coordinates to 1D index.
 * @param i x-coordinate (0 <= i <= nx)
 * @param j y-coordinate (0 <= j <= ny)
 * @param nx Number of interior grid points in x direction.
 */
int mapping(int i, int j, int nx) {
    return j * (nx + 1) + i;
}

/**
 * @brief Create a 2D Poisson matrix in CSR format.
 * @param nx Number of interior grid points in x direction.
 * @param ny Number of interior grid points in y direction.
 * @return Pointer to the generated SparseCSR matrix.
 */
SparseCSR* createMatrix(int nx, int ny) {
    int rows = (nx + 1) * (ny + 1);
    int cols = (nx + 1) * (ny + 1);
    int nnz = 5 * rows; // Approximate number of non-zeros

    SparseCSR *matrix = createSparseCSR(rows, cols, nnz);
    int *row_ptr = matrix->row_ptr;
    int *col_ind = matrix->col_ind;
    double *values = matrix->values;

    int idx = 0;
    row_ptr[0] = 0;

    for (int j = 0; j <= ny; j++) {
        for (int i = 0; i <= nx; i++) {
            // Boundary
            if (i == 0 || i == nx || j == 0 || j == ny) {
                col_ind[idx] = mapping(i, j, nx);
                values[idx] = 1.0;
                idx++;
                row_ptr[mapping(i, j, nx) + 1] = idx;
                continue;
            }

            // Interior
            else{
                // Down
                col_ind[idx] = mapping(i, j - 1, nx);
                values[idx] = -1.0;
                idx++;

                // Left
                col_ind[idx] = mapping(i - 1, j, nx);
                values[idx] = -1.0;
                idx++;

                // Center
                col_ind[idx] = mapping(i, j, nx);
                values[idx] = 4.0;
                idx++;

                // Right
                col_ind[idx] = mapping(i + 1, j, nx);
                values[idx] = -1.0;
                idx++;

                // Up
                col_ind[idx] = mapping(i, j + 1, nx);
                values[idx] = -1.0;
                idx++;
            }

            row_ptr[mapping(i, j, nx) + 1] = idx;
        }
    }

    matrix->nnz = idx; // Update nnz in case of boundary conditions
    return matrix;
}

/**
 * @brief Main function to generate and display a 2D Poisson matrix in CSR format.
 * @param argc Argument count.
 * @param argv Argument vector. Expects two arguments: nx and ny.
 * @return int Exit status of the program.
 * @details
 * The main function generates a 2D Poisson matrix based on the provided grid size,
 * displays the matrix, and extracts its diagonal, lower, and upper triangular parts.
 * The results are printed to the console.
 * @note This program requires two command-line arguments: nx and ny.
 * @note nx > 6 or ny > 6 will result in a simplified matrix display.
 * Usage:
 * ```
 * mkdir build && cd build
 * cmake ..
 * make
 * ../bin/2D-Poisson_2d <nx> <ny>
 * ```
 */
int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Usage: %s <nx> <ny>\n", argv[0]);
        return -1;
    }
    int nx = atoi(argv[1]);
    int ny = atoi(argv[2]);

    SparseCSR *matrix = createMatrix(nx, ny);
    printf("Generated 2D Poisson matrix in CSR format:\n");
    if (nx <= 6 && ny <= 6) {
        print_SparseCSR(matrix, 1);
    }
    else {
        printf("Matrix too large to display fully.\n");
        print_SparseCSR_simple(matrix, 1);
    }
    SparseCSR** diag_matrices = get_D_L_U_csr(matrix); // Example usage of get_diag_L_U_csr
    SparseCSR* D = diag_matrices[0];
    SparseCSR* L = diag_matrices[1];
    SparseCSR* U = diag_matrices[2];
    printf("Diagonal matrix D:\n");
    if (nx <= 6 && ny <= 6) {
        print_SparseCSR(D, 1);
    } else {
        printf("Matrix too large to display fully.\n");
        print_SparseCSR_simple(D, 1);
    }
    printf("Lower triangular matrix L:\n");
    if (nx <= 6 && ny <= 6) {
        print_SparseCSR(L, 1);
    } else {
        printf("Matrix too large to display fully.\n");
        print_SparseCSR_simple(L, 1);
    }
    printf("Upper triangular matrix U:\n");
    if (nx <= 6 && ny <= 6) {
        print_SparseCSR(U, 1);
    } else {
        printf("Matrix too large to display fully.\n");
        print_SparseCSR_simple(U, 1);
    }
    // Free allocated memory for D, L, U
    freeSparseCSR(D);
    freeSparseCSR(L);
    freeSparseCSR(U);
    free(diag_matrices);

    // Clean up
    freeSparseCSR(matrix);
    return 0;
}