/**
 * @file csr.h
 * @brief Header file for CSR matrix operations and iterative solvers.
 * 
 * This header file declares the SparseCSR structure and functions for
 * creating, freeing, decomposing, performing matrix-vector multiplication,
 * and solving linear systems using iterative methods (Jacobi, Gauss-Seidel,
 * Conjugate Gradient) for sparse matrices stored in Compressed Sparse Row (CSR)
 * format.
 * @see csr.c
 * @author Li Zhijun
 * @date 2025-10-10
 */
# ifndef CSR_H
# define CSR_H
# include "vec.h"

/**
 * @struct SparseCSR
 * @brief Structure to represent a sparse matrix in Compressed Sparse Row (CSR) format.
 * 
 * CSR format stores a sparse matrix using three arrays:
 * - row_ptr: Array of size (rows + 1) where row_ptr[i] indicates the start of row i in col_ind and values.
 * - col_ind: Array of size nnz (number of non-zero elements) storing the column indices of non-zero elements.
 * - values: Array of size nnz storing the non-zero values of the matrix.
 */
typedef struct {
    int rows;       /**< Number of rows in the matrix. */
    int cols;       /**< Number of columns in the matrix. */
    int nnz;        /**< Number of non-zero elements in the matrix. */
    int *row_ptr;   /**< Row pointer array of size 'rows + 1'. */
    int *col_ind;   /**< Column index array of size 'nnz'. */
    double *values; /**< Non-zero values array of size 'nnz'. */
} SparseCSR;

/**
 * @brief Create a new SparseCSR matrix structure.
 * @param rows Number of rows.
 * @param cols Number of columns.
 * @param nnz Number of non-zero elements.
 * @return Pointer to the newly allocated SparseCSR structure.
 * 
 * @note The caller is responsible for freeing the allocated memory using freeSparseCSR().
 * @see freeSparseCSR()
 */
SparseCSR* createSparseCSR(int rows, int cols, int nnz);

/**
 * @brief Free the memory allocated for a SparseCSR matrix.
 * @param matrix Pointer to the SparseCSR structure to free.
 */
void freeSparseCSR(SparseCSR *matrix);

/**
 * @brief Sparse matrix-vector multiplication (y = A*x) for CSR format.
 * @param matrix Pointer to the SparseCSR matrix.
 * @param x Input vector.
 * @param y Output vector (result).
 */
void spmv_csr(const SparseCSR *matrix, const double *x, double *y);

/**
 * @brief Decompose a CSR matrix into diagonal, lower, and upper matrices.
 *
 * Returns an array of three SparseCSR pointers:
 *   - result[0]: Diagonal matrix (D)
 *   - result[1]: Lower triangular matrix (L)
 *   - result[2]: Upper triangular matrix (U)
 *
 * @param matrix Pointer to the input SparseCSR matrix.
 * @return Array of three SparseCSR* (D, L, U)
 */
SparseCSR** get_D_L_U_csr(const SparseCSR *matrix);

/**
 * @brief Solve Ax = b using the Jacobi iterative method for CSR matrices.
 * @param matrix Pointer to the SparseCSR matrix (A).
 * @param b Right-hand side vector.
 * @param x Solution vector (input: initial guess, output: result).
 * @param max_iter Maximum number of iterations.
 * @param tol Tolerance for convergence.
 */
void Jacobi_csr(const SparseCSR *matrix, const double *b, double *x, int max_iter, double tol);

/**
 * @brief Solve Ax = b using the Gauss-Seidel iterative method for CSR matrices.
 * @param matrix Pointer to the SparseCSR matrix (A).
 * @param b Right-hand side vector.
 * @param x Solution vector (input: initial guess, output: result).
 * @param max_iter Maximum number of iterations.
 * @param tol Tolerance for convergence.
 */
void GaussSeidel_csr(const SparseCSR *matrix, const double *b, double *x, int max_iter, double tol);

/**
 * @brief Solve Ax = b using the Conjugate Gradient method for CSR matrices.
 * @param matrix Pointer to the SparseCSR matrix (A).
 * @param b Right-hand side vector.
 * @param x Solution vector (input: initial guess, output: result).
 * @param max_iter Maximum number of iterations.
 * @param tol Tolerance for convergence.
 */
void CG_csr(const SparseCSR *matrix, const double *b, double *x, int max_iter, double tol);

# endif