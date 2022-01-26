
#include "vector.h"

#include "error_msg.h"
#include "format_printer.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>


vector* alloc_vector(int size, int initialize, int is_int){
    vector *res;

    res = malloc(sizeof(vector));

    if (res == NULL){
        print_err_exit();
    }

    if (is_int == VEC_IS_INT){
        if (initialize == VEC_INIT){
            res->block1 = calloc(size, sizeof(int));
        }else{
            res->block1 = malloc(size * sizeof(int));
        }

        if (res->block1 == NULL){
            print_err_exit();
        }
    }else{
        if (initialize == VEC_INIT){
            res->block = calloc(size, sizeof(double));
        }
        else{ /* do not initialize the vector. */
            res->block = malloc(size * sizeof(double));
        }

        if (res->block == NULL){
        print_err_exit();
        }
    }
    
    res->d = size;
    res->is_int = is_int;

    return res;
}

void free_vector(vector *vec){
    if (vec->is_int == VEC_IS_INT){
        /* need to check first if this points at NULL! */
        free(vec->block1);
    }
    else{
        /* need to check first if this points at NULL! */
        free(vec->block);
    }
    free(vec);
}

int vec_get_size(vector *vec){
    return vec->d;
}

void vector_resize(vector *vec, int new_size){

    if (vec->is_int == VEC_IS_INT){
        vec->block1 = realloc(vec->block1, new_size * sizeof(int));
    }else{
        vec->block = realloc(vec->block, new_size * sizeof(double));
    }

    vec->d = new_size;
}

vector* array_to_vec(void *arr, int size, int is_int){
    vector* res;
    int i;
    int *arr1;
    double *arr2;

    if ((res = malloc(sizeof(vector))) == NULL){
        print_err_exit();
    }

    if (is_int == VEC_IS_INT){
        arr1 = (int *) arr;
        res = alloc_vector(size, ~VEC_INIT, VEC_IS_INT);
        for (i=0; i < size; i++){
            vector_set(res, i, arr1[i]);
        }
    }
    else{
        arr2 = (double *) arr;
        res = alloc_vector(size, ~VEC_INIT, ~VEC_IS_INT);
        for (i=0; i < size; i++){
            vector_set(res, i, arr2[i]);
        }
    }

    return res;
}

void *vector_to_Array(vector *vec){
    int *res, i, d;
    double *res1;

    d = vec_get_size(vec);

    if (vec->is_int == VEC_IS_INT){
        res = malloc(sizeof(int) * d);
        if (res == NULL){
            print_err_exit();
        }
        for (i=0; i < d; i++){
            res[i] = vector_get(vec, i);
        }

        return (void *) res;
    }else{
        res1 = malloc(sizeof(double) * d);
        if (res1 == NULL){
            print_err_exit();
        }
        for (i=0; i < d; i++){
            res1[i] = vector_get(vec, i);
        }

        return (void *) res1;
    }
}

double vector_get(vector *vec, int i){
    if ((i >= vec_get_size(vec)) || (i < 0)){
        printf("Error at vector_get(): index outranges vector's dimesnion.\n");
        exit(1);
    }

    if (vec->is_int == VEC_IS_INT){
        return ((double) (vec->block1)[i]);
    }else{
        return (vec->block)[i];
    }
}

void vector_set(vector *vec, int i, double value){
    if ((i >= vec_get_size(vec)) || (i < 0)){
        printf("Error at vector_set(): index %d outranges vector's dimension %d.\n", 
                    i, vec_get_size(vec));
        exit(1);
    }

    if (vec->is_int == VEC_IS_INT){
        (vec->block1)[i] = (int) value;
    }else{
        (vec->block)[i] = value;
    }   
}

double vec_distance(vector *x1, vector *x2){
    int i;
    double res; res = 0;
    
    if (vec_get_size(x1) != vec_get_size(x2)){
        printf("Error at distance()");
        exit(1);
    }
    
    for(i=0; i < x1->d; i++){
        res += pow(vector_get(x1, i) - vector_get(x2, i), 2);
    }

    return sqrt(res);
}

double vec_sum(vector *vec){
    int i;
    double res; res = 0;

    for (i=0; i < vec->d; i++){
        res += vector_get(vec, i);
    }

    return res;
}

double vec_dot(vector *x1, vector *x2){
    int i;
    double res; res = 0;

    if (vec_get_size(x1) != vec_get_size(x2)){
        printf("Error: at vec_dot(). vectors dimensions inequal.\n");
        exit(1);
    }

    for (i=0; i < x1->d; i++){
        res += vector_get(x1, i) * vector_get(x2, i);
    }

    return res;
}

void print_vector(vector *vec){
    int i, d;

    d = vec_get_size(vec);

    for (i=0; i < d-1; i++){
        print_double(vector_get(vec, i), ',', WITH_SEP);
    }

    print_double(vector_get(vec, d-1), '\n', WITH_SEP);
}

vector* int_vec(vector *vec){
    vector* res;
    int i, n;

    n = vec_get_size(vec);

    res = alloc_vector(n, ~VEC_INIT, VEC_IS_INT);

    for (i=0; i < n; i++){
        vector_set(res, i, (int) vector_get(vec, i));
    }

    return res;
}
