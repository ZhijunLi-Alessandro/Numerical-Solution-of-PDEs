/**
 * @file test_csr_3x3.c
 * @brief Test program for 3x3 CSR matrix operations and iterative solvers.
 * 
 * @details
 * This test program creates a 3x3 sparse matrix in CSR format, defines a right-hand side vector,
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
   ../bin/test_csr_3x3 \endverbatim
 * @note This file will be regarded as a test file in doxygen documentation in next edition.
 * @see csr.h, utils.h
 * @author Li Zhijun
 * @date 2025-10-10
 * @test test_csr_3x3.c
 */
# include <stdio.h>
# include "csr.h"
# include "utils.h"

/**
 * @brief Main function to test 3x3 CSR matrix operations and iterative solvers.
 * 
 * @details 
 * The main function creates a 3x3 sparse matrix in CSR format, defines a right-hand side vector,
 * and solves the linear system using Jacobi, Gauss-Seidel, and Conjugate Gradient methods.
 * The results are printed to the console.
 * 
 * example matrix:
 * \verbatim
   [ 4 -1  0]
   [-1  4 -1]
   [ 0 -1  3]
   b = [15, 10 10]^T
   x = [5 5 5]^T \endverbatim
 * 
 * @return int Exit status of the program.
 */
int main() {

    int rows = 3, cols = 3, nnz = 7;
    SparseCSR *matrix = createSparseCSR(rows, cols, nnz);
    int row_ptr[] = {0, 2, 5, 7};
    int col_ind[] = {0, 1, 0, 1, 2, 1, 2};
    double values[] = {4.0, -1.0, -1.0, 4.0, -1.0, -1.0, 3.0};
    for (int i = 0; i < rows + 1; i++) {
        matrix->row_ptr[i] = row_ptr[i];
    }
    for (int i = 0; i <nnz; i++) {
        matrix->col_ind[i] = col_ind[i];
        matrix->values[i] = values[i];
    }
    double b[] = {15.0, 10.0, 10.0};
    printf("Matrix A in CSR format:\n");
    print_SparseCSR(matrix, 6);
    printf("Right-hand side vector b:\n");
    print_vector(b, rows, 6);
    double x_Jacobi[] = {0.0, 0.0, 0.0};
    double x_GS[] = {0.0, 0.0, 0.0};
    double x_CG[] = {0.0, 0.0, 0.0};
    int max_iter = 50;
    double tol = 1e-6;
    Jacobi_csr(matrix, b, x_Jacobi, max_iter, tol);
    printf("Jacobi Solution x:\n");
    print_vector(x_Jacobi, rows, 6);
    GaussSeidel_csr(matrix, b, x_GS, max_iter, tol);
    printf("Gauss-Seidel Solution x:\n");
    print_vector(x_GS, rows, 6);
    CG_csr(matrix, b, x_CG, max_iter, tol);
    printf("Conjugate Gradient Solution x:\n");
    print_vector(x_CG, rows, 6);
    freeSparseCSR(matrix);
    return 0;
}