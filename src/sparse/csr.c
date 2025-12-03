/**
 * @file csr.c
 * @brief Implementation of sparse matrix operations in CSR format.
 * 
 * This file includes functions for creating, freeing, decomposing, matrix-vector
 * multiplication, and solving linear systems using iterative methods (Jacobi,
 * Gauss-Seidel, Conjugate Gradient) for sparse matrices stored in Compressed
 * Sparse Row (CSR) format.
 * 
 * @author Li Zhijun
 * @date 2025-10-10
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "csr.h"

SparseCSR* createSparseCSR(int rows, int cols, int nnz) {
    SparseCSR *matrix = (SparseCSR *)malloc(sizeof(SparseCSR));
    matrix->rows = rows;
    matrix->cols = cols;
    matrix->nnz = nnz;
    matrix->row_ptr = (int *)malloc((rows + 1) * sizeof(int));
    matrix->col_ind = (int *)malloc(nnz * sizeof(int));
    matrix->values = (double *)malloc(nnz * sizeof(double));
    return matrix;
}

void freeSparseCSR(SparseCSR *matrix) {
    if (matrix) {
        free(matrix->row_ptr);
        free(matrix->col_ind);
        free(matrix->values);
        free(matrix);
    }
}

void spmv_csr(const SparseCSR *matrix, const double *x, double *y) {
    for (int i = 0; i < matrix->rows; i++) {
        y[i] = 0.0;
        for (int j = matrix->row_ptr[i]; j < matrix->row_ptr[i + 1]; j++) {
            y[i] += matrix->values[j] * x[matrix->col_ind[j]];
        }
    }
}

SparseCSR** get_D_L_U_csr(const SparseCSR *matrix) {
    int rows = matrix->rows;
    int cols = matrix->cols;
    int nnz_diag = rows;
    int nnz_L = 0, nnz_U = 0;
    // First pass: count nnz for L and U
    for (int i = 0; i < rows; i++) {
        for (int j = matrix->row_ptr[i]; j < matrix->row_ptr[i + 1]; j++) {
            int col = matrix->col_ind[j];
            if (col < i) nnz_L++;
            else if (col > i) nnz_U++;
        }
    }

    SparseCSR *diag = createSparseCSR(rows, cols, nnz_diag);
    SparseCSR *L = createSparseCSR(rows, cols, nnz_L);
    SparseCSR *U = createSparseCSR(rows, cols, nnz_U);

    // Initialize row_ptr
    diag->row_ptr[0] = 0;
    L->row_ptr[0] = 0;
    U->row_ptr[0] = 0;
    int idx_diag = 0, idx_L = 0, idx_U = 0;
    for (int i = 0; i < rows; i++) {
        // Diagonal always has one entry per row
        diag->col_ind[idx_diag] = i;
        diag->values[idx_diag] = 0.0; // will be set below
        idx_diag++;

        for (int j = matrix->row_ptr[i]; j < matrix->row_ptr[i + 1]; j++) {
            int col = matrix->col_ind[j];
            double val = matrix->values[j];
            if (col < i) {
                L->col_ind[idx_L] = col;
                L->values[idx_L] = val;
                idx_L++;
            } else if (col == i) {
                diag->values[i] = val;
            } else if (col > i) {
                U->col_ind[idx_U] = col;
                U->values[idx_U] = val;
                idx_U++;
            }
        }
        diag->row_ptr[i + 1] = idx_diag;
        L->row_ptr[i + 1] = idx_L;
        U->row_ptr[i + 1] = idx_U;
    }

    SparseCSR **result = (SparseCSR **)malloc(3 * sizeof(SparseCSR*));
    result[0] = diag;
    result[1] = L;
    result[2] = U;
    return result;
}

void Jacobi_csr_debug(const SparseCSR *matrix, const double *b, double *x, int max_iter, double tol) {
    double *x_new = (double *)malloc(matrix->rows * sizeof(double));
    for (int iter = 0; iter < max_iter; iter++) {
        for (int i = 0; i < matrix->rows; i++) {
            double sum = 0.0;
            double diag = 0.0;
            for (int j = matrix->row_ptr[i]; j < matrix->row_ptr[i + 1]; j++) {
                if (matrix->col_ind[j] == i) {
                    diag = matrix->values[j];
                } else {
                    sum += matrix->values[j] * x[matrix->col_ind[j]];
                }
            }
            x_new[i] = (b[i] - sum) / diag;
        }
        double norm = 0.0;
        for (int i = 0; i < matrix->rows; i++) {
            norm += (x_new[i] - x[i]) * (x_new[i] - x[i]);
            x[i] = x_new[i];
        }
        norm = sqrt(norm);
        printf("Jacobi Iteration %d: Residual = %e\n", iter + 1, norm);
        if (norm < tol) break;
    }
    free(x_new);
}

void Jacobi_csr(const SparseCSR *matrix, const double *b, double *x, int max_iter, double tol) {
    double *x_new = (double *)malloc(matrix->rows * sizeof(double));
    for (int iter = 0; iter < max_iter; iter++) {
        for (int i = 0; i < matrix->rows; i++) {
            double sum = 0.0;
            double diag = 0.0;
            for (int j = matrix->row_ptr[i]; j < matrix->row_ptr[i + 1]; j++) {
                if (matrix->col_ind[j] == i) {
                    diag = matrix->values[j];
                } else {
                    sum += matrix->values[j] * x[matrix->col_ind[j]];
                }
            }
            x_new[i] = (b[i] - sum) / diag;
        }
        double norm = 0.0;
        for (int i = 0; i < matrix->rows; i++) {
            norm += (x_new[i] - x[i]) * (x_new[i] - x[i]);
            x[i] = x_new[i];
        }
        norm = sqrt(norm);
        if (norm < tol) break;
    }
    free(x_new);
}

void GaussSeidel_csr_debug(const SparseCSR *matrix, const double *b, double *x, int max_iter, double tol) {
    for (int iter = 0; iter < max_iter; iter++) {
        double norm = 0.0;
        for (int i = 0; i < matrix->rows; i++) {
            double sum = 0.0;
            double diag = 0.0;
            for (int j = matrix->row_ptr[i]; j < matrix->row_ptr[i + 1]; j++) {
                if (matrix->col_ind[j] == i) {
                    diag = matrix->values[j];
                } else {
                    sum += matrix->values[j] * x[matrix->col_ind[j]];
                }
            }
            double x_old = x[i];
            x[i] = (b[i] - sum) / diag;
            norm += (x[i] - x_old) * (x[i] - x_old);
        }
        norm = sqrt(norm);
        printf("GS Iteration %d: Residual = %e\n", iter + 1, norm);
        if (norm < tol) break;
    }
}

void GaussSeidel_csr(const SparseCSR *matrix, const double *b, double *x, int max_iter, double tol) {
    for (int iter = 0; iter < max_iter; iter++) {
        double norm = 0.0;
        for (int i = 0; i < matrix->rows; i++) {
            double sum = 0.0;
            double diag = 0.0;
            for (int j = matrix->row_ptr[i]; j < matrix->row_ptr[i + 1]; j++) {
                if (matrix->col_ind[j] == i) {
                    diag = matrix->values[j];
                } else {
                    sum += matrix->values[j] * x[matrix->col_ind[j]];
                }
            }
            double x_old = x[i];
            x[i] = (b[i] - sum) / diag;
            norm += (x[i] - x_old) * (x[i] - x_old);
        }
        norm = sqrt(norm);
        if (norm < tol) break;
    }
}

void CG_csr_debug(const SparseCSR *matrix, const double *b, double *x, int max_iter, double tol) {
    double *r = (double *)malloc(matrix->rows * sizeof(double));
    double *p = (double *)malloc(matrix->rows * sizeof(double));
    double *Ap = (double *)malloc(matrix->rows * sizeof(double));

    // r = b - A*x
    spmv_csr(matrix, x, r);
    for (int i = 0; i < matrix->rows; i++) {
        r[i] = b[i] - r[i];
        p[i] = r[i];
    }

    double rsold = vec_dot(r, r, matrix->rows);

    for (int iter = 0; iter < max_iter; iter++) {
        spmv_csr(matrix, p, Ap);
        double alpha = rsold;
        double pAp = vec_dot(p, Ap, matrix->rows);
        alpha /= pAp;

        for (int i = 0; i < matrix->rows; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        double rsnew = vec_dot(r, r, matrix->rows);

        printf("CG Iteration %d: Residual = %e\n", iter + 1, sqrt(rsnew));

        if (sqrt(rsnew) < tol) break;

        for (int i = 0; i < matrix->rows; i++) {
            p[i] = r[i] + (rsnew / rsold) * p[i];
        }
        rsold = rsnew;
    }

    free(r);
    free(p);
    free(Ap);
}

void CG_csr(const SparseCSR *matrix, const double *b, double *x, int max_iter, double tol) {
    double *r = (double *)malloc(matrix->rows * sizeof(double));
    double *p = (double *)malloc(matrix->rows * sizeof(double));
    double *Ap = (double *)malloc(matrix->rows * sizeof(double));

    // r = b - A*x
    spmv_csr(matrix, x, r);
    for (int i = 0; i < matrix->rows; i++) {
        r[i] = b[i] - r[i];
        p[i] = r[i];
    }

    double rsold = vec_dot(r, r, matrix->rows);

    for (int iter = 0; iter < max_iter; iter++) {
        spmv_csr(matrix, p, Ap);
        double alpha = rsold;
        double pAp = vec_dot(p, Ap, matrix->rows);
        alpha /= pAp;

        for (int i = 0; i < matrix->rows; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        double rsnew = vec_dot(r, r, matrix->rows);

        if (sqrt(rsnew) < tol) break;

        for (int i = 0; i < matrix->rows; i++) {
            p[i] = r[i] + (rsnew / rsold) * p[i];
        }
        rsold = rsnew;
    }

    free(r);
    free(p);
    free(Ap);
}