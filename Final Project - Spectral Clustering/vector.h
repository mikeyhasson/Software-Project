
#ifndef VECTOR_H
#define VECTOR_H

#define VEC_INIT 1
#define VEC_IS_INT 1

typedef struct{
    int d;
    double *block;
    int *block1;
    int is_int;
}vector;

vector* alloc_vector(int size, int initialize, int is_int);

void free_vector(vector *vec);

int vec_get_size(vector *vec);

void vector_resize(vector *vec, int new_size);

vector* array_to_vec(void *arr, int size, int is_int);

void *vector_to_Array(vector *vec);

double vector_get(vector *vec, int i);

void vector_set(vector *vec, int i, double value);

double vec_distance(vector *x1, vector *x2);

double vec_sum(vector *vec);

vector* int_vec(vector *vec);

double vec_dot(vector *x1, vector *x2);

void print_vector(vector *vec);

#endif
