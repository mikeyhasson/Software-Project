#include "matrix.h"

#include "vector.h"
#include "error_msg.h"
#include "format_printer.h"

#include <math.h>


matrix* alloc_matrix(int m, int n, int init){
    struct matrix *mat;

    mat = malloc(sizeof(matrix));
    if (mat == NULL){
        print_err_exit();
    }

    if (init == INITIALIZE){
        mat -> block = calloc(m * n, sizeof(double));
    }else{ 
        /* dont intialize matrix's data to zeros */
        mat -> block = malloc(m * n * sizeof(double));
    }
    
    if (mat -> block == NULL){
        print_err_exit();
    }

    mat -> d1 = m;
    mat -> d2 = n;

    return mat;
}

void free_matrix(matrix *mat){
    free(mat -> block);
    free(mat);
}

int mat_get_nrows(matrix *mat){
    return mat->d1;
}

int mat_get_ncolumns(matrix *mat){
    return mat->d2;
}

vector* mat_get_row(matrix *mat, int i){

    if ((i < 0) || (i >= mat->d1)){
        printf("mat_get_row() Error: Invalid row index. index: %d, nrows: %d.\n", i, mat->d1);
        exit(1);
    }

    return array_to_vec((void *)(&((mat->block)[i*(mat->d2)])), mat->d2, ~ VEC_IS_INT);
}

double mat_get(matrix *mat, int i, int j){

    if ((i < 0) || (i >= mat->d1)){
        printf("mat_get() Error: Invalid row index. index: %d, nrows: %d.\n", i, mat->d1);
        exit(1);
    }
    if ((j < 0) || (j >= mat->d2)){
        printf("mat_get() Error: Invalid column index. index: %d, ncolumns: %d.\n", j, mat->d2);
        exit(1);
    }

    return (mat -> block)[i * (mat->d2) + j];
}

matrix* vec_to_diag(vector *vals){
    int n, i;
    matrix *res;

    n = vec_get_size(vals);
    res = alloc_matrix(n, n, INITIALIZE);

    for(i=0; i < n; i++){
        mat_set(res, i, i, vector_get(vals, i));
    }

    return res;
}

void mat_set(matrix *mat, int i, int j, double val){

    if ((i < 0) || (i >= mat->d1)){
        printf("Invalid row index.\n");
        exit(1);
    }
    if ((j < 0) || (j >= mat->d2)){
        printf("Invalid column index.\n");
        exit(1);
    }

    (mat -> block)[i * (mat->d2) + j] = val;
}

void mat_set_row(matrix *mat, int i, vector *vec){
    int j;

    if ((i < 0) || (i >= mat->d1)){
        printf("Error at mat_set_row(): row index out of range.\n");
        exit(1);
    }

    for (j=0; j < mat->d2; j++){
        mat_set(mat, i, j, vector_get(vec, j));
    }
}

matrix* diag(int n, double val){
    int i;
    matrix *I;

    if (n < 0){
        printf("Invalid dimension.\n");
        exit(1);
    }

    I = alloc_matrix(n, n, INITIALIZE);

    for (i = 0; i < n; i++){
        mat_set(I, i, i, val);
    }

    return I;
}

matrix* mat_mul(matrix *A, matrix *B){
    int i, j, k;
    matrix *C;
    double dot;

    if (A->d2 != B->d1){
        printf("Matrix multiplaction dimension error. \n");
        exit(1);
    }

    C = alloc_matrix(A->d1, B->d2, ~INITIALIZE);

    for (i=0; i < (A -> d1); i++){
        for (j = 0; j < (B -> d2); j++){
            dot = 0;
            for (k = 0; k < (A -> d2); k++){
                dot += mat_get(A, i, k) * mat_get(B, k, j);
            }
            mat_set(C, i, j, dot);
        }
    }

    return C;
}

matrix* mat_add(matrix *A, matrix *B){
    matrix *C;
    int i, j;

    if (A->d1 != B->d1){
        printf("Error at mat_mul: A->d1 != B->d1. /n");
        exit(1);
    }
    if (A->d2 != B->d2){
        printf("Error at mat_mul: A->d2 != B->d2. /n");
        exit(1);
    }

    C = alloc_matrix(A->d1, A->d2, ~INITIALIZE);

    for (i = 0; i < A->d1; i++){
        for (j = 0; j < A->d2; j++){
            mat_set(C, i, j, mat_get(A, i, j) + mat_get(B, i, j));
        }
    }

    return C;
}

matrix* mat_sub(matrix *A, matrix *B){
    matrix *C;
    int i, j;

    if (A->d1 != B->d1){
        printf("Error at mat_mul: A->d1 != B->d1. /n");
        exit(1);
    }
    if (A->d2 != B->d2){
        printf("Error at mat_mul: A->d2 != B->d2. /n");
        exit(1);
    }

    C = alloc_matrix(A->d1, A->d2, ~INITIALIZE);

    for (i = 0; i < A->d1; i++){
        for (j = 0; j < A->d2; j++){
            mat_set(C, i, j, mat_get(A, i, j) - mat_get(B, i, j));
        }
    }

    return C;
}

void add_vector_to_row(matrix *A, int i, vector *vec){
    int j;

    if (vec->d != A->d2){
       printf("Error at: add_Vector_to_Row(): d!=d2");
       exit(1);
    }
    
    for (j=0; j < A->d2; j++){
        mat_set(A, i, j, mat_get(A, i, j) + vector_get(vec, j));
    }
}

void divide_rows_by_vec(matrix *A, vector *vec){
    int i, j;
    
    if (vec->d != A->d1){
        printf("Error at divide_rows_by_Vec(): vec->d != A->d1");
        exit(1);
    }

    for (i=0; i < A->d1; i++){
        for (j=0; j < A->d2; j++){
            mat_set(A, i, j, mat_get(A, i, j) / vector_get(vec, i));
        }
    }
}

int matrix_eq(matrix *A, matrix *B){
    int i, j;

    if ((A -> d1) != (B -> d1)){
        printf("Error at are_matrix_eq(): A -> d1 != B -> d1");
        exit(1);
    }

    if ((A -> d2) != (B -> d2)){
        printf("Error at are_matrix_eq(): A -> d2 != B -> d2");
        exit(1);
    }

    for (i=0; i < A->d1; i++){
        for (j=0; j < A->d2; j++){
            if (mat_get(A, i, j) != mat_get(B, i, j)){
                return 0;
            }
        }
    }

    return 1;
}

matrix *transpose(matrix *A){
    matrix *res;

    int m, n, i, j;

    m = mat_get_nrows(A);
    n = mat_get_ncolumns(A);

    res = alloc_matrix(n, m, ~INITIALIZE);

    for (i=0; i < n; i++){
        for (j=0; j < m; j++){
            mat_set(res, i, j, mat_get(A, j, i));
        }
    }

    return res;
}

void set_row_to_vec(matrix *A, int i, vector *vec){
    int j;

    if (vec->d != A ->d2){
        printf("Error: dimension error at set_row_to_vec(): vec->d != A->d2.\n");
        exit(1);
    }
    
    for (j=0; j < A->d2; j++){
        mat_set(A, i, j, vector_get(vec, j));
    }
}

void copy_to_matrix(matrix *to, matrix *from, int rows, int cols){
    int i, j;

    if ((rows > from->d1) || (rows > to->d1) || (rows < 0)){
        printf("Error at copy_to_matrix: more rows than possible");
        exit(1);
    }
    if ((cols > from->d2) || (cols > to->d2) || (cols < 0)){
        printf("Error at copy_to_matrix: more cols than possible");
        exit(1);
    }

    for (i=0; i < rows; i++){
        for (j=0; j < cols; j++){
            mat_set(to, i, j, mat_get(from, i, j));
        }
    }
}

void mat_print(matrix *A){
    int i, j, n, d;

    n = mat_get_nrows(A); d = mat_get_ncolumns(A);

    for (i=0; i < n-1; i++){
        for (j=0; j < d-1; j++){
            print_double(mat_get(A, i, j), ',', WITH_SEP);
        }
        print_double(mat_get(A, i, d-1), '\n', WITH_SEP);
    }
    for (j=0; j < d-1; j++){
        print_double(mat_get(A, n-1, j), ',', WITH_SEP);
    }
    
    print_double(mat_get(A, n-1, d-1), 0, ~ WITH_SEP);
}

/* -- jacobi's algorithm implementation -- */

/* given matrix A, return 1 if and only if it's not diagonal. otherwise - 0 */
int not_diagonal(matrix *A){
    int i,j,n=A->d1;

    for (i=0; i < n; i++){
        for (j=0; j < n;j++){
            if ((i != j) && mat_get(A, i, j) != 0){
                return 1;
            }
        }
    }
    return 0;
}

/* calculates the index of the off-diagonal element with the largest absolute value */
void get_max_off_diagonal(matrix *A, int *i_max, int *j_max){
    int i,j,n=A -> d1;
    double A_max = fabs(mat_get(A,0,1));
    *i_max=0;
    *j_max=1;

    for (i=0;i< n;i++){
        for (j=i+1;j<n;j++){
            if ((i!=j) && fabs(mat_get(A,i,j))>A_max){
                A_max=fabs(mat_get(A,i,j));
                *i_max=i;
                *j_max=j;
            }   
        }
    }
}

/* calculates the rotation matrix P */
void build_rotation_matrix(matrix *P, matrix *A,int *i_max_ptr,int *j_max_ptr, double *c_ptr,double *s_ptr){
    int i_max,j_max;
    double theta,t,c,s;

    get_max_off_diagonal(A, &i_max, &j_max);
    theta = (mat_get(A,j_max,j_max)-mat_get(A,i_max,i_max))/(2*mat_get(A,i_max,j_max));
    t= 1 /(fabs(theta) + sqrt(pow(theta, 2) + 1));
    t= (theta < 0) ? ((-1)*t) : t;
    c= 1 /(sqrt(pow(t, 2)+1));
    s= t * c;
    if (*i_max_ptr != -1){ /* resetting P to eye. */
        mat_set(P,*i_max_ptr,*i_max_ptr,1);
        mat_set(P,*j_max_ptr,*j_max_ptr,1);
        mat_set(P,*i_max_ptr,*j_max_ptr,0);
        mat_set(P,*j_max_ptr,*i_max_ptr,0);
    }
    /* new i,j indices. */
    mat_set(P,i_max,i_max,c);
    mat_set(P,j_max,j_max,c);
    mat_set(P,i_max,j_max,s);
    mat_set(P,j_max,i_max,(-1)*s);

    *i_max_ptr=i_max;
    *j_max_ptr=j_max;
    *c_ptr=c;
    *s_ptr=s;
}

/* applies the transformation of A=P^T*A*P as described in jacobi algorithm */
double get_similar_matrix(matrix *A,int i,int j,double c,double s){
    int r, n;
    double Ari_old,Arj_old,Ari_new,Arj_new,Aii,Ajj,dif;

    dif = 0;
    n = mat_get_nrows(A);

    for (r=0; r < n; r++){
        if (!(r==i || r==j)){
            Ari_old= mat_get(A,r,i);
            Arj_old = mat_get(A,r,j);
            Ari_new= c*Ari_old-s*Arj_old;
            Arj_new = c*Arj_old+s*Ari_old;

            dif+=pow(Ari_old,2)+pow(Arj_old,2)-pow(Ari_new,2)-pow(Arj_new,2);
            mat_set(A,i,r,Ari_new);
            mat_set(A,r,i,Ari_new);
            mat_set(A,r,j,Arj_new);
            mat_set(A,j,r,Arj_new);
        }
    }

    Aii=mat_get(A,i,i);
    Ajj=mat_get(A,j,j);
    mat_set(A,i,i, pow(c,2)*Aii+pow(s,2)*Ajj-2*s*c*mat_get(A,i,j));
    mat_set(A,j,j, pow(s,2)*Aii+pow(c,2)*Ajj+2*s*c*mat_get(A,i,j));
    dif += pow(mat_get(A,i,j), 2);
    mat_set(A,i,j,0);
    mat_set(A,j,i,0);
    dif*=2;
    return dif;
}

typedef struct eigen {
    vector* vector;
    double val;
    int original_index;
} eigen;

/* comparator function for struct eigen */
int compare_eigens (const void * a, const void * b) {
    double c;

    if ((c = ( (((eigen*)a) ->val) - (((eigen*)b) ->val))) > 0){
        return 1;
    }else{
        if (c < 0){
            return -1;
        }else{
            /* c == 0. */

            return (((eigen*)a)->original_index) - (((eigen*)b)->original_index);
            
            /* we want the sort to be stable; */
        }
    }
}

/* sorts eigen values and eigen vectors */
void sort_eigen(vector *eigen_values, matrix *eigen_vectors){
    int i,j,n;
    eigen *eigen_lst;

    n = mat_get_nrows(eigen_vectors);

    eigen_lst = malloc(n * sizeof(eigen));

    if (eigen_lst == NULL){
        print_err_exit();
    }

    for (i = 0; i < n; i++){
        eigen_lst[i].val = vector_get(eigen_values, i);
        eigen_lst[i].vector = alloc_vector(n, ~ VEC_INIT, ~ VEC_IS_INT);
        eigen_lst[i].original_index = i;
    
        for(j=0; j<n; j++){
            vector_set(eigen_lst[i].vector, j, mat_get(eigen_vectors, j, i));
        }
    }

    qsort(eigen_lst, n, sizeof(eigen), compare_eigens);

    for (i=0; i<n; i++){
        vector_set(eigen_values, i, eigen_lst[i].val);

        for(j=0; j < n; j++){
            mat_set(eigen_vectors, j, i, vector_get(eigen_lst[i].vector, j));
        }

        free_vector(eigen_lst[i].vector);
    }

    free(eigen_lst);
}

/* description in the header file. */
void Find_eig(matrix *A, vector **eig_values, matrix **eig_vectors, 
                double epsilon, int jacobi_max_iter)
{
    int n, i, i_piv, j_piv,cnt;
    double c, s, diff;
    matrix *P, *B, *temp;

    n = mat_get_nrows(A);

    P = diag(n, 1); /* eye matrix */
    B = alloc_matrix(n, n, ~INITIALIZE);

    copy_to_matrix(B, A, n, n);

    *eig_vectors = diag(n, 1); /* eye matrix */

    i_piv = -1;
    cnt = 0;
    diff = epsilon + 1111;

    while ((diff > epsilon) && (cnt < jacobi_max_iter) && not_diagonal(B)){
        build_rotation_matrix(P, B, &i_piv, &j_piv, &c, &s);
        temp = mat_mul(*eig_vectors, P);
        free_matrix(*eig_vectors);
        *eig_vectors = temp;
        diff = get_similar_matrix(B, i_piv, j_piv, c, s);
        cnt++;
    }
    
    for (i=0; i<n; i++){
        vector_set(*eig_values, i, mat_get(B, i, i));
    }

    free_matrix(P);
    free_matrix(B);
}

void Find_eig_Sorted(matrix *A, vector **eig_values, matrix **eig_vectors,
             double epsilon, int jacobi_max_iter)
    {
        Find_eig(A, eig_values, eig_vectors, epsilon, jacobi_max_iter);

        sort_eigen(*eig_values, *eig_vectors);
    }
