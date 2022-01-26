#ifndef KMEANS_H
#define KMEANS_H

#include "matrix.h"
#include <stdlib.h>

/* returns: j at {0, ..., rows_offset-1} such that (*dist)(centroids[j], x) is minimal; */
int find_closest_centroid(matrix *centroids, vector *x, int rows_offset);

matrix* kmeans(matrix *data_points, int k, long max_iterations, int *initial_centroids_indices);

#endif
