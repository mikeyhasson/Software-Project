#ifndef NSC_H
#define NSC_H

#include "matrix.h"

/* calculates the Weigthed Adjacency matrix given the data points. */
matrix* wam(matrix *data);

/* calculates the Diagonal Degree matrix given the data points. */
matrix* ddg(matrix *data);

/* calculates the Normalized Graph Laplacian matrix given the data points. */
matrix* lnorm(matrix *data);

/* 
calculates the eig vals and vecs of a given matrix an puts them in the given pointers.
*/
void jacobi(matrix *mat, vector **eig_Values, matrix **eig_Vectors, 
            double epsilon, int jacobi_max_iter);

/* 
calculates the suitable centroids given the data points according to the full 
spectral clustering algorithm.
The variable "initial_centroids_indices" used as the indices of the initial_centorids 
of the matrix T in the spk algorithm.
*/
matrix* spk(matrix *data, int k, int *initial_centroids_indices, 
            double epsilon, int jacboi_max_iter, int kmeans_max_iter);

/*
just like the previous function, but returns a vector of indices correspond to the 
cluster of each of the data points, after clustered by the spk algorithm.
*/
vector* spk_cluster(matrix *data, int k, int *initial_centroids_indices, 
                    double epsilon, int jacobi_max_iter, int kmeans_max_iter);

#endif
