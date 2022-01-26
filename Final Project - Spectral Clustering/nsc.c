
#include "nsc.h"

#include "kmeans.h"
#include "error_msg.h"
#include "matrix.h"

#include <math.h>

/* ------- Normalized Spectral Clustering Algorithm's Auxiliary Functions -------*/

matrix* Calc_W(matrix *Data){
    int i, j, n;
    matrix *W;
    vector *x_i, *x_j;

    n = mat_get_nrows(Data);

    W = alloc_matrix(n, n, ~INITIALIZE);

    for(i=0; i < n; i++){
        for(j=(i+1) ; j < n; j++){
            x_i = mat_get_row(Data, i);
            x_j = mat_get_row(Data, j);
            mat_set(W, i, j, 
                exp((-1/(double) 2) * vec_distance(x_i, x_j)));

            free_vector(x_i);
            free_vector(x_j);
        }
    }

    for (i = 0; i < n; i++){
        mat_set(W, i, i, 0);
    }

    for (i = 0; i < n; i++){
        for (j = 0; j < i; j++){
            mat_set(W, i, j, mat_get(W, j, i));
        }
    }

    return W;
}

vector *Calc_DDG_diag(matrix *W){
    vector *diag, *x_i;
    int i;

    diag = alloc_vector(mat_get_nrows(W), ~VEC_INIT, ~VEC_IS_INT);

    for(i=0; i < mat_get_nrows(W); i++){
        x_i = mat_get_row(W, i);
        vector_set(diag, i, vec_sum(x_i));

        free_vector(x_i);
    }

    return diag;
}

matrix* Calc_Diagonal_Deg_Matrix(matrix *W){
    matrix *DDG;
    vector *diag;
    
    diag = Calc_DDG_diag(W);

    DDG = vec_to_diag(diag);

    free_vector(diag);

    return DDG;
}

matrix* Calc_D_Minus_Half(matrix *D){
    int i, n;
    matrix *D_minus_half;

    n = mat_get_nrows(D);

    D_minus_half = alloc_matrix(n, n, INITIALIZE);

    for (i=0; i < n; i++){
        mat_set(D_minus_half, i, i, ((double) 1) / sqrt(mat_get(D, i, i)));
    }

    return D_minus_half;
}

matrix* Calc_L_Norm(matrix *W, vector *DDG_diag){
    int n, i, j;
    matrix *L_norm;

    n = mat_get_nrows(W);

    L_norm = alloc_matrix(n, n, ~ INITIALIZE);

    for (i=0; i < n; i++){
        for (j=0; j < i; j++){
            mat_set(L_norm, i, j, 
            (-1) * mat_get(W, i, j) / (sqrt(vector_get(DDG_diag, i)) * sqrt(vector_get(DDG_diag, j))));
        }
        mat_set(L_norm, i, i, 1);
    }

    /* L_norm is symmetric */

    for (i=0; i < n; i++){
        for (j=i+1; j < n; j++){
            mat_set(L_norm, i, j, mat_get(L_norm, j, i));
        }
    }

    return L_norm;
}

int Find_k(vector *eig_values){
    int max_val, arg_max, i, offset, curr_val, n;

    n = vec_get_size(eig_values);

    offset = n / 2; /* rounds towards |_n/2_|*/

    arg_max = 1;
    max_val = vector_get(eig_values, 1) - vector_get(eig_values, 0);

    for (i = 2; i <= offset; i ++){
        curr_val = vector_get(eig_values, i) - vector_get(eig_values, i-1);

        if (max_val < curr_val){
            max_val = curr_val;
            arg_max = i;
        }
    }

    return arg_max;
}

matrix* build_U(matrix *eig_vectors, int k){
    int n;
    matrix *U;

    n = mat_get_ncolumns(eig_vectors); 
    /* The dimension of each eigen vector (of L_norm) */

    U = alloc_matrix(n, k, ~INITIALIZE);

    copy_to_matrix(U, eig_vectors, n, k);

    return U;
}

matrix* build_T(matrix *U){
    int n, k, i, j;
    double norm;
    matrix *T;
    vector *x_i;

    n = mat_get_nrows(U); 
    k = mat_get_ncolumns(U);

    T = alloc_matrix(n, k, 0);
    
    for (i=0; i < n; i++){
        x_i = mat_get_row(U, i);
        norm = sqrt(vec_dot(x_i, x_i));
        for (j=0; j < k; j++){
            mat_set(T, i, j, mat_get(U, i, j) / norm);
        }

        free_vector(x_i);
    }

    return T;
}

vector* spk_clusters_by_centroids(matrix *data, matrix *final_centroids){
    int i, n, k;
    vector *clusters, *x_i;

    n = mat_get_nrows(data);
    k = mat_get_nrows(final_centroids);

    clusters = alloc_vector(n, ~ VEC_INIT, VEC_IS_INT);

    for (i=0; i < n; i++){
        x_i = mat_get_row(data, i);

        vector_set(clusters, 
                    i, 
                    find_closest_centroid(final_centroids, x_i, k));

        free_vector(x_i);
    }

    return clusters;
}

/* ------- NSC Algorithm's Functions Implementation ------- */

matrix* wam(matrix *data){
    return Calc_W(data);
}

matrix* ddg(matrix *data){
    matrix *W, *DD;
    
    W = Calc_W(data);
    DD = Calc_Diagonal_Deg_Matrix(W);

    free_matrix(W);
    
    return DD;
}

matrix* lnorm(matrix *data){
    matrix *W, *L_norm;
    vector *diag;
  
    W = Calc_W(data);
    diag = Calc_DDG_diag(W);

    L_norm = Calc_L_Norm(W, diag);

    free_matrix(W);
    free_vector(diag);

    return L_norm;
}

void jacobi(matrix *mat, vector **eig_Values, matrix **eig_Vectors, double epsilon,
            int jacobi_max_iter)
{   
    int n;

    n = mat_get_nrows(mat);

    if (n != mat_get_ncolumns(mat)){
        printf("Error at jacobi(). matrix must be square.\n");
        exit(1);
    }

    if ((vec_get_size(*eig_Values) != n) || (mat_get_nrows(*eig_Vectors) != n) || 
        (mat_get_ncolumns(*eig_Vectors) != n))
    {
        printf("Error at jacobi(). given parameters dimension are incompatible.\n");
        exit(1);
    }

    Find_eig(mat, eig_Values, eig_Vectors, epsilon, jacobi_max_iter);
}

matrix* spk(matrix *data, int k, int *initial_centroids_indices,
            double epsilon, int jacobi_max_iter, int kmeans_max_iter)
{
    matrix *W, *L_norm, *eig_Vectors, *U, *T, *final_centroids;
    vector *eig_Values, *DDG_diag;
    int n;

    n = mat_get_nrows(data);

    W = Calc_W(data);

    DDG_diag = Calc_DDG_diag(W);

    L_norm = Calc_L_Norm(W, DDG_diag);
    
    free_matrix(W); 
    free_vector(DDG_diag);

    eig_Vectors = alloc_matrix(n, n, ~INITIALIZE);
    eig_Values = alloc_vector(n, ~VEC_INIT, ~VEC_IS_INT);

    jacobi(L_norm, &eig_Values, &eig_Vectors, epsilon, jacobi_max_iter);
    sort_eigen(eig_Values, eig_Vectors);

    if (k == 0){
        k = Find_k(eig_Values);
    }

    free_vector(eig_Values);

    U = build_U(eig_Vectors, k);

    T = build_T(U);

    free_matrix(U);

    final_centroids = kmeans(T, k, kmeans_max_iter, initial_centroids_indices);
    /* initial_centroids_indices might be NULL */

    free_matrix(eig_Vectors);
    free_matrix(T);

    return final_centroids;
}

vector* spk_cluster(matrix *data, int k, int *initial_centroids_indices, 
                    double epsilon, int jacobi_max_iter, int kmeans_max_iter)
{
    matrix *final_centroids;
    vector *clusters;

    final_centroids = spk(data, 
                            k, 
                            initial_centroids_indices, 
                            epsilon, 
                            jacobi_max_iter,
                            kmeans_max_iter);

    clusters = spk_clusters_by_centroids(data, final_centroids);

    free_matrix(final_centroids);

    return clusters;
}
