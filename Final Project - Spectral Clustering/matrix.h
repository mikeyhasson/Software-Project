
#ifndef MATRIX_H
#define MATRIX_H

#include "vector.h"

#include <stdlib.h>
#include <stdio.h>

#define INITIALIZE 1

typedef struct matrix{
    int d1;
    int d2;
    double *block;
} matrix;

matrix* alloc_matrix(int m, int n, int init);

void free_matrix(matrix *mat);

matrix* vec_to_diag(vector *vals);

int mat_get_nrows(matrix *mat);

int mat_get_ncolumns(matrix *mat);

double mat_get(matrix *mat, int i, int j);

vector* mat_get_row(matrix *mat, int i);

void mat_set(matrix *mat, int i, int j, double val);

void mat_set_row(matrix *mat, int i, vector *vec);

matrix* diag(int n, double val);

matrix* mat_mul(matrix *A, matrix *B);

matrix* mat_add(matrix *A, matrix *B);

matrix* mat_sub(matrix *A, matrix *B);

void add_vector_to_row(matrix *A, int i, vector *vec);

void divide_rows_by_vec(matrix *A, vector *vec);

int matrix_eq(matrix *A, matrix *B);

matrix *transpose(matrix *A);

void set_row_to_vec(matrix *A, int i, vector *vec);

void copy_to_matrix(matrix *to, matrix *from, int rows, int cols);

void mat_print(matrix *A);

/* 
receives a matrix A, 2 pointers for storing the result (eig_values and eig_vectors) 
of the calculation of the matrix A eigen values and their corresponding eigen vectors.
epsilon, jacobi_max_iter are parameters controlling on the accuracy of the result.
*/
void Find_eig(matrix *A, vector **eig_values, matrix **eig_vectors,
             double epsilon, int jacobi_max_iter);

/* 
receives a matrix A, 2 pointers for storing the result (eig_values and eig_vectors) 
of the calculation of the matrix A eigen values and their corresponding eigen vectors,
sorted in a descending order w.r.t the eigen values.
epsilon, jacobi_max_iter are parameters controlling on the accuracy of the result.
*/
void Find_eig_Sorted(matrix *A, vector **eig_values, matrix **eig_vectors,
             double epsilon, int jacobi_max_iter);

/* sorts eigen values and eigen vectors */
void sort_eigen(vector *eigen_values, matrix *eigen_vectors);

/* receives a filename to the data, and returns a matrix of the scanned data */
matrix* filename_to_matrix(char *filename);

#endif
