#define PY_SSIZE_T_CLEAN
#include <Python.h>
#define PY3K

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define DEF_MAX_ITERATIONS 200
#define INITIAL_DATA_CONTAINER_SIZE 200
#define INITIAL_DATA_DIM 24

#define DEBUG
#define SHOW_ITER_NUMBER
#define SHOW_GROUPS_SIZES
#define PRINT_RESULT

#undef DEBUG
#undef SHOW_ITER_NUMBER
#undef SHOW_GROUPS_SIZES


static double **my_k_means(double **data_points, int n, int k, int d, long max_iterations, double (*dist)(double *, double *, int), double **initial_centroids);
void alloc_2d_array(int dim1, int dim2, int zeros, double ***arr, double **block);
double euclidian_dist(double *x1, double *x2, int d);
int find_closest_centroid(double **centroids, double *x, int k, int d, double (*dist)(double *, double *, int));
void add_to_vector(double *array, double *adder, int d);
void divide_vectors_by_integers(double **vectors, int *integers, int dim1, int dim2);
int cmp_2d_arrays(double **arr1, double **arr2, int dim1, int dim2);
void print_2d_array(double **arr, int dim1, int dim2);
void print_error(char *err_msg);
double **Scan_Input(int *size, int *dim);
void copy_to_array(double *buffer, double *data, int size);
void copy_to_matrice(double **buffer, double **data, int dim1, int dim2);
void print_vector(char *msg, int *vector, int size);


int main(int argc, char **argv){
    
    double **Data, **centroids;
    int n, k, d;
    long max_iterations;

    if (!((argc == 2) || (argc == 3))){
        printf("Invalid Input!\n");
        assert(0);
    }
    if (argc == 3){
        if (!(atof(argv[2]) == ((int)atof(argv[2])))){
            printf("Invalid Input!\n");
            assert(0);
        }
    }
    max_iterations = (argc == 2) ? DEF_MAX_ITERATIONS : atof(argv[2]);

    if (!(atof(argv[1]) == ((int)atof(argv[1])))){
        printf("Invalid Input!\n");
        assert(0);
    }
    k = atof(argv[1]);

    /* Scanning the Data which is a 2nd dimensional matrix of size (n, k) which will be put into those variables
    after the next command; */
    Data = Scan_Input(&n, &d);

    
    if (!((n > k) && (max_iterations > 0) && (k > 0))){
        printf("Invalid Input!\n");
        assert(0);
    }

    /* else: */

    #ifdef DEBUG
    printf("\nThe size of the data is (%d, %d).\n", n, d);
    printf("----------------------------------\nThe Data is:\n");
    print_2d_array(Data, n, d);
    printf("\n");
    #endif

    centroids = my_k_means(Data, n, k, d, max_iterations, &euclidian_dist, NULL);
    #ifdef PRINT_RESULT
    print_2d_array(centroids, k, d);
    #endif

    return 0;
}

/* 
    pre: input is valid;
    returns: an 2 dimensional array of centroids whose size is (d, k);
*/
static double **my_k_means(double **data_points, int n, int k, int d, long max_iterations, double (*dist)(double *, double *, int), double **initial_centroids){
    double **curr_centroids, **next_centroids, *x_j, *curr_centroids_allocated, *next_centroids_allocated;
    /* The last two are must when we want to free up the memory allocated for those 2 dimensional matrices; */
    int updated, i, j, argmin_j, *next_sizes;

    if (initial_centroids == NULL){
        /* k rows of vectors at R^d, all uninitialized; these are the current loop's centroids; */
        alloc_2d_array(k, d, 0, &curr_centroids, &curr_centroids_allocated); 
        
        /* mue0 <- x0, mue1 <- x1, ......., mue_(k-1) <- x_(k-1); */
        for (i=0; i < k; i++){
            for (j=0; j < d; j++){
                curr_centroids[i][j] = data_points[i][j];
            }
        }
    }  
    else{ /* The initial centroids are given */
        curr_centroids = initial_centroids;
    }

    /* k rows of vectors at R^d, all initialized to zeros; these are the next loop's centroids; */
    alloc_2d_array(k, d, 1, &next_centroids, &next_centroids_allocated);

    /* k counters, all unitialized; these are the next loop's groups sizes; */
    next_sizes = calloc(k, sizeof(int));
    if (next_sizes == NULL){
        printf("Fail when trying to allocate memory\n");
        assert(0);
    }

    #ifdef DEBUG
    printf("\n-------------------------\nMax Iterations: %ld\n-------------------------\n", max_iterations);
    #endif

    /* loop intializers; */
    i = 0; updated = 1; 
    
    while (updated && (i < max_iterations)){
        #ifdef DEBUG
        printf("Loop Number: %d\nCurrent Centroids: \n", i);
        print_2d_array(curr_centroids, k, d);
        printf("\n");
        #endif
        #ifdef SHOW_ITER_NUMBER
        printf("Loop Number: %d\n", i);
        #endif

        for (j=0; j < n; j++){
            x_j = data_points[j];
            argmin_j = find_closest_centroid(curr_centroids, x_j, k, d, dist);
            add_to_vector(next_centroids[argmin_j], x_j, d);
            next_sizes[argmin_j] ++;
        }

        #ifdef DEBUG
        printf("Next Centroids: \n");
        print_2d_array(next_centroids, k, d);
        printf("\n");
        print_vector("Groups Sizes", next_sizes, k);
        printf("\n");
        #endif
        
        /* next_centroids[i] <- next_centroids[i] / next_sizes[i] (where the first is a vector and the second is a scalar); */
        divide_vectors_by_integers(next_centroids, next_sizes, k, d);
        #ifdef DEBUG
        printf("Division's Result: \n");
        print_2d_array(next_centroids, k, d);
        printf("\n-------------------------------\n");
        #endif
        /* checking if mue = (mue_0, ......., mue_(k-1)) has changed; */
        updated = cmp_2d_arrays(next_centroids, curr_centroids, k, d);

        /* freeing up the curr_centroids; */
        free(curr_centroids);
        if (((initial_centroids == NULL) && (i == 0)) || (i > 0)){  
            free(curr_centroids_allocated);
        }

        /* setting up the next iteration's centroids; */
        curr_centroids = next_centroids;
        curr_centroids_allocated = next_centroids_allocated;
        
        /* setting up the next next iteration's centroids and group sizes; */
        alloc_2d_array(k, d, 1, &next_centroids, &next_centroids_allocated); 
        #ifdef SHOW_GROUPS_SIZES
        if (updated == 0){
            print_vector("Groups Sizes", next_sizes, k);
            printf("\n");
        }
        #endif
        next_sizes = calloc(k, sizeof(int));
        if (next_sizes == NULL){
            printf("Fail when trying to allocate memory\n");
            assert(0);
        }
        i ++;
    }

    free(next_sizes);
    free(next_centroids);
    free(next_centroids_allocated);
    #ifdef DEBUG
    printf("The Result is: \n\n");
    #endif
    return curr_centroids;
}


/* 
    allocates a contigous 2dimensional array of size (dim1, dim2), of type double;
    The array values are initialized to zeros IF AND ONLY IF (zeros != 0);
*/
void alloc_2d_array(int dim1, int dim2, int zeros, double ***alloc_array, double **alloc_block){
    double *array, **result;
    int i;

    /* A contigous array of doubles; */
    array = (zeros == 0) ? (double *)malloc(dim1 * dim2 * sizeof(double)) : (double *)calloc(dim1*dim2, sizeof(double)); 
    if (array == NULL){
        printf("Fail when trying to allocate memory\n");
        assert(0);
    }
    result = (double **)malloc(dim1 * sizeof(double *)); /* An array of (double *); */
    if (result == NULL){
        printf("Fail when trying to allocate memory\n");
        assert(0);
    }

    if ((array == NULL) || (result == NULL)) /* Memory allocation failed; */
        return;

    for (i=0; i < dim1; i++){
        result[i] = array + (i * dim2);
    }

    *alloc_array = result; /* The contigous 2 dimensional array; */
    *alloc_block = array; /* The malloc'd pointer to its pointer; */
}

double euclidian_dist(double *x1, double *x2, int d){
    int i;
    double Result;

    Result = 0;

    for (i=0; i < d; i++){
        Result += (x2[i] - x1[i]) * (x2[i] - x1[i]);
    }

    return Result;
}

/* returns: j at {0, ..., k-1} such that (*dist)(centroids[j], x) is minimal; */
int find_closest_centroid(double **centroids, double *x, int k, int d, double (*dist)(double *, double *, int)){
    int j, argmin_j;
    double min_dist, distance;

    min_dist = (*dist)(centroids[0], x, d); argmin_j = 0;
    for (j=1; j < k; j++){
        distance = (*dist)(centroids[j], x, d);
        if (distance < min_dist){
            argmin_j = j;
            min_dist = distance;
        }
    }

    return argmin_j;
}

void add_to_vector(double *array, double *adder, int d){
    int i;
    for (i=0; i < d; i++){
        array[i] += adder[i];
    }
}

void divide_vectors_by_integers(double **vectors, int *integers, int dim1, int dim2){
    int i, j;
    for (i=0; i < dim1; i++){
        for (j=0; j < dim2; j++){
            vectors[i][j] = vectors[i][j] / integers[i];
        }
    }
}

/* returns: 0 if equal and 1 if not; */
int cmp_2d_arrays(double **arr1, double **arr2, int dim1, int dim2){
    int i, j;
    for (i=0; i < dim1; i++){
        for (j=0; j < dim2; j++){
            if (arr1[i][j] != arr2[i][j]){
                return 1;
            }
        }
    }
    return 0;
}

void print_2d_array(double **arr, int dim1, int dim2){
    int i, j;
    for (i=0; i < dim1; i++){
        for (j = 0; j < dim2-1; j++){
            printf("%.4f,", arr[i][j]);
        }
        printf("%.4f\n", arr[i][dim2 - 1]); 
    }

    /* last line will be empty;  */
}

void print_error(char *err_msg){
    fprintf(stdout, "\nError: %s\n", err_msg);
}

/* returns: the scanned (double **) data matrix scanned, and puts its dimensions at the pointers "size" and "dime"; */
double **Scan_Input(int *size, int *dim){
    int i, n, d, buffer_size;
    char curr_char;
    double curr_value, *first_data_point , **arr, *block, **tmp_arr, *tmp_block;

    /* scanning the first data point and its dimensional - d; */
    buffer_size = INITIAL_DATA_DIM;
    first_data_point = (double *)malloc(buffer_size * sizeof(double));
    if (first_data_point == NULL){
        printf("Fail when trying to allocate memory\n");
        assert(0);
    }
    d = 0;

    while (scanf("%lf%c", &curr_value, &curr_char) == 2){
        if (d == buffer_size){
            buffer_size *= 2;
            first_data_point = (double *)realloc(first_data_point ,buffer_size * sizeof(double));
            if (first_data_point == NULL){
                printf("Fail when trying to allocate memory\n");
                assert(0);
            }
        }
        first_data_point[d] = curr_value;
        d++;
        if (curr_char == '\n'){
            break;
        }
    }

    /* scanning all of the data points and their count - n; */
    buffer_size = INITIAL_DATA_CONTAINER_SIZE;
    n = 1; i = 0;
    alloc_2d_array(buffer_size, d, 0, &arr, &block);
    copy_to_array(arr[0], first_data_point, d);

    while (scanf("%lf%c", &curr_value, &curr_char) == 2){
        if (n == buffer_size){
            buffer_size *= 2;
            alloc_2d_array(buffer_size, d, 0, &tmp_arr, &tmp_block);
            copy_to_matrice(tmp_arr, arr, n, d);
            /* freeing up previously used memory; */
            free(arr);
            free(block);
            arr = tmp_arr;
            block = tmp_block;
        }

        arr[n][i] = curr_value;
        if (curr_char == '\n'){
            n ++;
            i = -1; /* will now be increamented to 0; */
        }
        i++;
    }    

    *size = n; 
    *dim = d;
    return arr; /* A 2 dimensional array of size (n, k) which contains the scanned data; */
}

void copy_to_array(double *buffer, double *data, int size){
    int i;
    for (i=0; i < size; i++){
        buffer[i] = data[i];
    }
}

void copy_to_matrice(double **buffer, double **data, int dim1, int dim2){
    int i, j;
    for (i=0; i < dim1; i++){
        for (j=0; j < dim2; j++){
            buffer[i][j] = data[i][j];
        }
    }
}

void print_vector(char *msg, int *vector, int size){
    int i;
    printf("%s: ", msg);
    for (i=0; i < size - 1; i++){
        printf("%d,", vector[i]);
    }
    printf("%d\n", vector[size-1]);
}


/*--------------------------------Connecting module to python--------------------------------*/

static double extract_item(PyObject *py_list, int index, int type);
static PyObject *pack_output(double **final_centroids, int k, int d);

static PyObject* fit(PyObject *self, PyObject *args){
    
    int n, k, d, i, j;
    double **data_points, **initial_centroids, **result;
    double *data_block, *initial_centroids_block;
    long num_of_iter;

    PyObject *py_list, *output;

    /* --- Parsing the arguments ---*/
    
    if (!PyArg_ParseTuple(args, "O", &py_list)){
        return NULL;
    }

    /* py_list will contain the data points in its first n rows and  
       the initial centroids in the rows n, n+1, ..., n+k */

    py_list = PySequence_Fast(py_list, NULL); 
    /* we know that the object is iterable so we do not care of the case when 
       it doesn't and therefore we dont need to decide on an error message */

    if (!py_list){
        return NULL;
    }

    /*--- excracting n, k, d, num_of_iter ---*/
    n = (int) extract_item(py_list, 0, 1);
    k = (int) extract_item(py_list, 1, 1);
    d = (int) extract_item(py_list, 2, 1);
    num_of_iter = (long) extract_item(py_list, 3, 1);

    /* --- Converting the data to a C list (of double **) --- */
    data_points = malloc(sizeof(double *) * n);
    data_block = malloc(sizeof(double) * n * d);

    if (data_block == NULL){
        Py_DECREF(py_list);
        return PyErr_NoMemory();
    }

    for (i=0; i < n; i++){
        data_points[i] = &(data_block[i * d]);
        for (j=0; j < d; j++){
            data_points[i][j] = extract_item(py_list, 4 + i*d + j, 0);
        }
    }

    /* --- Converting the initial centroids to a C list (of double **) --- */
    initial_centroids = malloc(sizeof(double *) * k);
    initial_centroids_block = malloc(sizeof(double) * k * d);
    
    for (i=n; i < n+k; i++){
        initial_centroids[i-n] = &(initial_centroids_block[(i-n) * d]);
        for (j=0; j < d; j++){
            initial_centroids[i-n][j] = extract_item(py_list, 4 + i*d + j, 0);
        }
    }

    Py_DECREF(py_list); /* Telling python's garbage collector to collect this (dereferencing) */

    /* --- Calculating the result using the written kmeans algorithm ---*/
    result = my_k_means(data_points, n, k, d, num_of_iter, &euclidian_dist, initial_centroids);

    output = pack_output(result, k, d);
    
    if (output == NULL){ /* indicates a memory error */
        return PyErr_NoMemory();
    }

    free(result);
   
    return output;
}

static double extract_item(PyObject *py_list, int index, int type){
    PyObject *item, *f;
    double res;

    item = PySequence_Fast_GET_ITEM(py_list, index);
    if (type == 0){
        f = PyNumber_Float(item);
        res = PyFloat_AsDouble(f);
    } else{
        f = PyNumber_Long(item);
        res = PyLong_AsLong(f);
    }

    Py_DECREF(f);

    return res;
}

static PyObject *pack_output(double **final_centroids, int k, int d){
    int i, j;
    PyObject *out_list, *item;

    out_list = PyList_New(k * d);
    if (!out_list){
        Py_DECREF(out_list);
        return NULL;
    }

    for (i = 0; i < k; i++){
        for (j = 0; j < d; j++){
            item = PyFloat_FromDouble(final_centroids[i][j]);
            if (item == NULL){
                Py_DECREF(out_list);
                Py_DECREF(item);
                return NULL;
            }
            /* else */

            if (PyList_SetItem(out_list, i*d + j, item) == -1){
                return NULL;
            }
        }
    }

    return out_list;
}


static PyMethodDef capiMethods[] = 
{
    {"kmeans",
    (PyCFunction) fit,
    METH_VARARGS,
    PyDoc_STR("The Kmeans algorithm using given initial centroids")
    },
    {NULL, NULL, 0, NULL}

};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    NULL,
    -1,
    capiMethods
};

PyMODINIT_FUNC
PyInit_mykmeanssp(void){
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m){
        return NULL;
    }
    return m;
}
