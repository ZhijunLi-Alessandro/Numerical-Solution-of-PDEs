/**
 * @file test_csr_5x5.c
 * @brief Test program for 5x5 CSR matrix operations and iterative solvers.
 * 
 * @details
 * This test program creates a 5x5 sparse matrix in CSR format, defines a right-hand side vector,
 * and solves the linear system using Jacobi, Gauss-Seidel, and Conjugate Gradient methods.
 * The results are printed to the console.
 * 
 * Usage:
 * Compile the program and run it.
 * 
 * Example:
 * \verbatim
   mkdir build && cd build
   cmake ..
   make
   ../bin/test_csr_5x5 \endverbatim
 * @note This file will be regarded as a test file in doxygen documentation in next edition.
 * @see csr.h, utils.h
 * @author Li Zhijun
 * @date 2025-10-10
 * @test test_csr_5x5.c
 */
# include <stdio.h>
# include "csr.h"
# include "utils.h"

/**
 * @brief Main function to test 5x5 CSR matrix operations and iterative solvers.
 * 
 * @details 
 * The main function creates a 5x5 sparse matrix in CSR format, defines a right-hand side vector,
 * and solves the linear system using Jacobi, Gauss-Seidel, and Conjugate Gradient methods.
 * The results are printed to the console.
 * 
 * example matrix:
 * \verbatim
 * [ 2 -1  0  0  0]
 * [-1  2 -1  0  0]
 * [ 0 -1  2 -1  0]
 * [ 0  0 -1  2 -1]
 * [ 0  0  0 -1  2]
 * b = [1 2 3 4 5]^T
 * x = [35/6, 32/3, 27/2, 40/3, 55/6]^T \endverbatim
 * 
 * @return int Exit status of the program.
 */
int main() {

    int rows = 5, cols = 5, nnz = 13;
    SparseCSR *matrix = createSparseCSR(rows, cols, nnz);
    int row_ptr[] = {0, 2, 5, 8, 11, 13};
    int col_ind[] = {0, 1, 0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4, 5, 4, 5};
    double values[] = {2.0, -1.0,  -1.0, 2.0, -1.0,  -1.0, 2.0, -1.0,  -1.0, 2.0, -1.0,  -1.0, 2.0};
    for (int i = 0; i < rows + 1; i++) {
        matrix->row_ptr[i] = row_ptr[i];
    }
    for (int i = 0; i <nnz; i++) {
        matrix->col_ind[i] = col_ind[i];
        matrix->values[i] = values[i];
    }
    double b[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    printf("Matrix A in CSR format:\n");
    print_SparseCSR(matrix, 6);
    printf("Right-hand side vector b:\n");
    print_vector(b, rows, 6);
    double x_Jacobi[] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double x_GS[] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double x_CG[] = {0.0, 0.0, 0.0, 0.0, 0.0};
    int max_iter = 50;
    double tol = 1e-6;
    Jacobi_csr_debug(matrix, b, x_Jacobi, max_iter, tol);
    printf("Jacobi Solution x:\n");
    print_vector(x_Jacobi, rows, 6);
    GaussSeidel_csr_debug(matrix, b, x_GS, max_iter, tol);
    printf("Gauss-Seidel Solution x:\n");
    print_vector(x_GS, rows, 6);
    CG_csr_debug(matrix, b, x_CG, max_iter, tol);
    printf("Conjugate Gradient Solution x:\n");
    print_vector(x_CG, rows, 6);
    freeSparseCSR(matrix);
    return 0;
}