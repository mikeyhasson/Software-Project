#define PY_SSIZE_T_CLEAN
#include <Python.h>

#define DEBUG

#include "spkmeans.h"
#include "nsc.h"
#include "matrix.h"
#include "kmeans.h"
#include "error_msg.h"
#include "file_translators.h"
#include "stdio.h"

/* --------------- Objects Transportation Functions --------------- */

double extract_item(PyObject *py_list, int index, int type){
    PyObject *item, *f;
    double res;

    f = NULL;

    item = PySequence_Fast_GET_ITEM(py_list, index);
    if (type == 0){
        f = PyNumber_Float(item);
        res = PyFloat_AsDouble(f);
    } else{
        f = PyNumber_Long(item);
        res = PyLong_AsLong(f);
    }

    if (f != NULL){
    	Py_DECREF(f);
    }

    return res;
}

PyObject* pack_vector(vector *vec){
    int i, n;
    PyObject *list, *item;

    item = NULL;
    n = vec_get_size(vec);    

    list = PyList_New(n);
    if (!list){
        Py_DECREF(list);
        return NULL;
    }

    for (i = 0; i < n; i++){
        item = PyFloat_FromDouble(vector_get(vec, i));
        if (item == NULL){
            Py_DECREF(list);
            Py_DECREF(item);
            return NULL;
        }
        /* else */

        if (PyList_SetItem(list, i, item) == -1){
            Py_DECREF(list);
            Py_DECREF(item);
            return NULL;
        }
    }

    return list;
}

vector* unpack_list_to_vec(PyObject *list, int size){
    int i;
    vector *vec;

    list = PySequence_Fast(list, NULL); 

    if (!list){
        return NULL;
    }

    /* --- Converting the data to a C vector --- */

    vec = alloc_vector(size, ~ VEC_INIT, ~ VEC_IS_INT);
    
    for (i=0; i < size; i++){
        vector_set(vec, i, extract_item(list, i, 0));
    }

    return vec;
}

/* packs the C matrix to a python list. */
PyObject* pack_matrix(matrix *mat){

    int i, j, n, d;
    PyObject *out_list, *item;

    item = NULL;
    n = mat_get_nrows(mat); d = mat_get_ncolumns(mat);

    out_list = PyList_New(n * d);
    if (!out_list){
        Py_DECREF(out_list);
        return NULL;
    }

    for (i = 0; i < n; i++){
        for (j = 0; j < d; j++){
            item = PyFloat_FromDouble(mat_get(mat, i, j));
            if (item == NULL){
                Py_DECREF(out_list);
                Py_DECREF(item);
                return NULL;
            }
            /* else */

            if (PyList_SetItem(out_list, i*d + j, item) == -1){
                Py_DECREF(out_list);
                Py_DECREF(item);
                return NULL;
            }
        }
    }

    return out_list;
}

/* unpacks a python list to a C matrix. */
matrix* unpack_list(PyObject *list, int n, int d){
    int i, j;
    matrix *mat;
    
    list = PySequence_Fast(list, NULL); 

    if (!list){
        return NULL;
    }

    /* --- Converting the data to a C matrix --- */

    mat = alloc_matrix(n, d, ~ INITIALIZE);
    
    for (i=0; i < n; i++){
        for (j=0; j < d; j++){
            mat_set(mat, i, j, extract_item(list, i*d + j, 0));
        }
    }

    return mat;
}


/* --------------- CAPI Interface Functions ----------------- */

static PyObject* wam_py(PyObject *self, PyObject *args){
    
    int n, d;
    matrix *data_points, *W;
    PyObject *data_list, *output;

    /* --- Parsing the arguments --- */
    if (!PyArg_ParseTuple(args, "Oii", &data_list, &n, &d)){
        return NULL;
    }

    data_points = unpack_list(data_list, n, d); /* extracting the list */

    /* --- Calculating the wam --- */
    W = wam(data_points);
    
    /* --- packing and returning the result --- */
    output = pack_matrix(W);

    free_matrix(W);
    free_matrix(data_points);

    return output;
}

static PyObject* ddg_py(PyObject *self, PyObject *args){
    int n, d;
    matrix *data_points, *DD;
    PyObject *data_list, *output;

    /* --- Parsing the arguments --- */
    if (!PyArg_ParseTuple(args, "Oii", &data_list, &n, &d)){
        return NULL;
    }

    data_points = unpack_list(data_list, n, d); /* unpacking the input to a C matrix */

    DD = ddg(data_points); /* building the ddg matrix from it */

    output = pack_matrix(DD); /* packing the C matrix to a python list */
    
    free_matrix(DD);
    free_matrix(data_points);

    return output;
}

static PyObject* lnorm_py(PyObject *self, PyObject *args){
    
    int n, d;
    matrix *data_points, *L_norm;
    PyObject *data_list, *output;

    /* --- Parsing the arguments --- */
    if (!PyArg_ParseTuple(args, "Oii", &data_list, &n, &d)){
        return NULL;
    }

    data_points = unpack_list(data_list, n, d); /* unpacking the input to a C matrix */

    L_norm = lnorm(data_points); /* building the L_norm matrix from it */

    output = pack_matrix(L_norm); /* packing the C matrix to a python list */
    
    free_matrix(L_norm);
    free_matrix(data_points);

    return output;
}

static PyObject* jacobi_py(PyObject *self, PyObject *args){

    int n;
    matrix *mat, *eig_Vectors;
    vector *eig_Values;
    PyObject *mat_list, *output, *eig_Values_list, *eig_Vectors_list;

    if (!PyArg_ParseTuple(args, "Oi", &mat_list, &n)){
        return NULL;
    }

    mat = unpack_list(mat_list, n, n);

    eig_Values = alloc_vector(n, ~ VEC_INIT, ~ VEC_IS_INT);
    eig_Vectors = alloc_matrix(n, n, ~INITIALIZE);

    jacobi(mat, &eig_Values, &eig_Vectors, EPSILON, JACOBI_MAX_ITER);

    eig_Values_list = pack_vector(eig_Values);
    eig_Vectors_list = pack_matrix(eig_Vectors);

    free_vector(eig_Values);
    free_matrix(eig_Vectors);
    free_matrix(mat);

    output = Py_BuildValue("(O,O)", eig_Values_list, eig_Vectors_list);

    return output;
}

static PyObject* sorted_jacobi_py(PyObject *self, PyObject *args){

    int n;
    matrix *mat, *eig_Vectors;
    vector *eig_Values;
    PyObject *mat_list, *output, *eig_Values_list, *eig_Vectors_list;

    if (!PyArg_ParseTuple(args, "Oi", &mat_list, &n)){
        return NULL;
    }

    mat = unpack_list(mat_list, n, n);

    eig_Values = alloc_vector(n, ~ VEC_INIT, ~ VEC_IS_INT);
    eig_Vectors = alloc_matrix(n, n, ~INITIALIZE);

    jacobi(mat, &eig_Values, &eig_Vectors, EPSILON, JACOBI_MAX_ITER);

    sort_eigen(eig_Values, eig_Vectors);

    eig_Values_list = pack_vector(eig_Values);
    eig_Vectors_list = pack_matrix(eig_Vectors);

    free_vector(eig_Values);
    free_matrix(eig_Vectors);
    free_matrix(mat);

    output = Py_BuildValue("(O,O)", eig_Values_list, eig_Vectors_list);

    return output;
}

static PyObject* kmeans_py(PyObject *self, PyObject *args){
    int n, d, k, *init_centroids_indices;
    matrix *data_points, *centroids;
    vector *temp, *temp2;
    PyObject *data_list, *out_list, *cents_list;

    if (!PyArg_ParseTuple(args, "OOiii", &data_list, &cents_list, &n, &d, &k)){
        return NULL;
    }

    data_points = unpack_list(data_list, n, d);
    temp = unpack_list_to_vec(cents_list, k);
    temp2 = int_vec(temp);
    free_vector(temp);
    init_centroids_indices = (int *) vector_to_Array(temp2);
    free_vector(temp2);

    centroids = kmeans(data_points, k, KMEANS_MAX_ITER, init_centroids_indices);

    free(init_centroids_indices);
    free_matrix(data_points);

    out_list = pack_matrix(centroids);

    free_matrix(centroids);

    return out_list;
}

/* ------------------ Setting a Connection For Python --------------------*/


static PyMethodDef capiMethods[] = {
    {"wam", (PyCFunction) wam_py, METH_VARARGS, PyDoc_STR("gets:Data points. returns: Their Weighted Adjacency matrix")},
    {"ddg", (PyCFunction) ddg_py, METH_VARARGS, PyDoc_STR("gets:Data points. returns: Their Diagonal Degree matrix.")},
    {"lnorm", (PyCFunction) lnorm_py, METH_VARARGS, PyDoc_STR("gets:Data points. returns: Their Normalized Graph Laplacian matrix.")},
    {"jacobi", (PyCFunction) jacobi_py, METH_VARARGS, PyDoc_STR("gets: A matrix. returns: its eig_Values and their eig_Vectors as a tuple of PyLists.")},
    {"sorted_jacobi", (PyCFunction) sorted_jacobi_py, METH_VARARGS, PyDoc_STR("gets: A matrix. returns: its eig_Values and their eig_Vectors as a tuple of PyLists, sorted in an accending order.")},
    {"kmeans", (PyCFunction) kmeans_py, METH_VARARGS, PyDoc_STR("gets: Data Points, num_clusters, init_cent_indices. returns: final centroids calculated by kmeans algorithm.")},
    
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "capi_spk", /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    capiMethods /* the PyMethodDef array from before containing the methods of the extension */
};

PyMODINIT_FUNC
PyInit_capi_spk(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}