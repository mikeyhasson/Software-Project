
#include "file_translators.h"
#include "matrix.h"

#include <stdio.h>

/* translates the filename onto a vector of integers, puts the number of the indices
   at the file in *ind_num, and returns the indices array */
int* filename_to_indices(char *filename, int *ind_num){
    FILE *fp;
    int curr_value;
    char curr_char;
    int buffer_size, *indices, *tmp_buff;
    int k, i;

    if ((fp = fopen(filename, "r+")) == NULL){
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

    *ind_num = k;

    return indices;
}

matrix* filename_to_matrix(char *filename){
    int i, n, d, buffer_size;
    char curr_char, last_char;
    double curr_value;
    vector *first_data_point;
    matrix *data_points, *tmp_mat;
    FILE *fp;

    if ((fp = fopen(filename, "r+")) == NULL){
        print_input_err_exit();
    }

    /* scanning the first data point and its dimensional - d; */
    buffer_size = INITIAL_DATA_DIM;
    /* first_data_point = (double *)malloc(buffer_size * sizeof(double)); */
    first_data_point = alloc_vector(buffer_size, ~VEC_INIT, ~VEC_IS_INT);

    d = 0;

    while (fscanf(fp, "%lf%c", &curr_value, &curr_char) == 2){
        if (d == buffer_size){
            buffer_size *= 2;
            vector_resize(first_data_point, 2 * buffer_size);
            /* first_data_point = (double *)realloc(first_data_point ,buffer_size * sizeof(double)); */
        }
        vector_set(first_data_point, d, curr_value);
        d++;
        if (curr_char == '\n'){
            break;
        }
    }

    /* scanning all of the data points and their count - n; */
    buffer_size = INITIAL_DATA_CONTAINER_SIZE;
    data_points = alloc_matrix(buffer_size, d, 0);

    mat_set_row(data_points, 0, first_data_point);
    free_vector(first_data_point);

    n = 1;
    i = 0;
    last_char = '\n';

    while (fscanf(fp, "%lf%c", &curr_value, &curr_char) == 2){
        if (n == buffer_size){
            buffer_size *= 2;
            tmp_mat = alloc_matrix(buffer_size, d, 0);
            copy_to_matrix(tmp_mat, data_points, n, d);
            free_matrix(data_points); /* freeing up previous data's memory; */
            data_points = tmp_mat;
        }

        mat_set(data_points, n, i, curr_value);
        if (curr_char == '\n'){
            n ++;
            i = -1; /* will now be increamented to 0; */
        }
        i++;
        last_char = curr_char;
    }

    /* scanning the last index at the last vector */
    fscanf(fp, "%lf", &curr_value);
    mat_set(data_points, n, d-1, curr_value);
    
    if (last_char != '\n'){
        n++;
    }
    
    if (buffer_size != n){
        /* reshape the data points matrix to the suitable dimensions */
        tmp_mat = alloc_matrix(n, d, 0);
        copy_to_matrix(tmp_mat, data_points, n, d);
        free_matrix(data_points);
        data_points = tmp_mat;
    }

    fclose(fp);

    return data_points; /* A matrix of shape (n, d) which contains the scanned data; */
}
