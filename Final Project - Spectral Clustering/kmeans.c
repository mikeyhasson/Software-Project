#include "kmeans.h"
#include "error_msg.h"
#include "vector.h"

matrix* init_initial_centroids(char *init_cent_indices_filename, matrix *data_points){
    FILE *fp;
    int curr_value;
    char curr_char;
    int buffer_size, *indices, *tmp_buff;
    int k, i, d;
    matrix *init_centroids;
    vector *row;

    d = mat_get_ncolumns(data_points);

    if ((fp = fopen(init_cent_indices_filename, "r+")) == NULL){
        print_input_err_exit();
    } 

    k = 0;
    buffer_size = 10;
    indices = malloc(buffer_size * sizeof(int));

    while (fscanf(fp, "%d%c", &curr_value, &curr_char) == 2){
        if (buffer_size == k){
            buffer_size *= 2;
            tmp_buff = malloc(buffer_size * sizeof(int));
            if (tmp_buff == NULL) {print_err_exit();}
            for (i=0; i < k; i++){
                tmp_buff[i] = indices[i];
            }
            free(indices);
            indices = tmp_buff;
        }

        indices[k] = curr_value;
        k++;
    }

    init_centroids = alloc_matrix(k, d, 0);

    for (i=0; i < k; i++){
        row = mat_get_row(data_points, indices[i]);
        mat_set_row(init_centroids, i, row);
        free_vector(row);
    }

    free(indices);

    return init_centroids;
}

/* receives the indices array and its length-k, and outputs the corresponding centroids*/
matrix* indices_to_centroids(int *indices, int k , matrix *data_points){
    matrix* init_centroids;
    int d, i;
    vector *row;

    d = mat_get_ncolumns(data_points);

    init_centroids = alloc_matrix(k, d, 0);

    for (i=0; i < k; i++){
        row = mat_get_row(data_points, indices[i]);
        mat_set_row(init_centroids, i, row);
        
        free_vector(row);
    }

    return init_centroids;
}

/* returns: j at {0, ..., rows_offset-1} such that (*dist)(centroids[j], x) is minimal; */
int find_closest_centroid(matrix *centroids, vector *x, int rows_offset){
    int j, argmin_j;
    double min_dist, dist;
    vector *cent;

    if (rows_offset > mat_get_nrows(centroids)){
        printf("Error at find_closest_centroid(). rows_offset > centroids->d1.\n");
        exit(1);
    }

    cent = mat_get_row(centroids, 0);
    min_dist = vec_distance(cent, x); 

    free_vector(cent);

    argmin_j = 0;

    for (j=1; j < rows_offset; j++){
        cent = mat_get_row(centroids, j);
        dist = vec_distance(cent, x);
        if (dist < min_dist){
            argmin_j = j;
            min_dist = dist;
        }

        free_vector(cent);
    }

    return argmin_j;
}

/* returns: a 2 dimensional array of centroids which's size is (k, d); */
matrix* kmeans(matrix *data_points, int k, long max_iterations, int *initial_centroids_indices){

    matrix *curr_centroids, *next_centroids;
    vector *x_j, *next_sizes, *row;
    /* The last two are must when we want to free up the memory allocated for those 2 dimensional matrices; */
    int updated, i, j, argmin_j, n, d;

    n = mat_get_nrows(data_points);
    d = mat_get_ncolumns(data_points);

    if (initial_centroids_indices == NULL){
        /* k rows of vectors at R^d, all uninitialized; these are the current loop's centroids; */
        curr_centroids = alloc_matrix(k, d, 0);
        
        /* mue0 <- x0, mue1 <- x1, ......., mue_(k-1) <- x_(k-1); */
        for (i=0; i < k; i++){
            row = mat_get_row(data_points, i);
            mat_set_row(curr_centroids, i, row);

            free_vector(row);
        }
    }
    else{ 
        /* The initial centroids are given */
        curr_centroids = indices_to_centroids(initial_centroids_indices, k, data_points);
    }

    /* k rows of vectors at R^d, all initialized to zeros; these are the next loop's centroids; */
    next_centroids = alloc_matrix(k, d, 1);

    /* k counters, all unitialized; these are the next loop's groups sizes; */
    next_sizes = alloc_vector(k, VEC_INIT, VEC_IS_INT);

    /* loop intializers; */
    i = 0; updated = 1; 
    
    while (updated && (i < max_iterations)){

        for (j=0; j < n; j++){
            x_j = mat_get_row(data_points, j);
            argmin_j = find_closest_centroid(curr_centroids, x_j, k);
            add_vector_to_row(next_centroids, argmin_j, x_j);
            vector_set(next_sizes, argmin_j, vector_get(next_sizes, argmin_j) + 1);

            free_vector(x_j);
        }
        
        /* next_centroids[i] <- next_centroids[i] / next_sizes[i] (where the first is a vector and the second is a scalar); */
        divide_rows_by_vec(next_centroids, next_sizes);
    
        /* checking if mue = (mue_0, ......., mue_(k-1)) has changed; */
        updated = ~ matrix_eq(next_centroids, curr_centroids);

        /* freeing up the curr_centroids; */
        free_matrix(curr_centroids);

        /* setting up the next iteration's centroids; */
        curr_centroids = next_centroids;
        
        /* setting up the next next iteration's centroids and group sizes; */
        next_centroids = alloc_matrix(k, d, 1);

        free_vector(next_sizes);
        next_sizes = alloc_vector(k, VEC_INIT, VEC_IS_INT);

        i ++;
    }

    free_vector(next_sizes);
    free_matrix(next_centroids);

    return curr_centroids;
}
