#include "spkmeans.h"

#include "nsc.h"
#include "matrix.h"
#include "kmeans.h"
#include "error_msg.h"
#include "file_translators.h"

#include <stdlib.h>
#include <string.h>

#include "format_printer.h"

/* ------------ Main Interface ------------*/

/* Note: It's the only function exceeds 80 lines of code. */
int main(int argc, char **argv){
    int k, n;
    char *goal;
    matrix *data_points, *W, *DD, *L_norm, *eig_Vectors, *final_centroids, *tmp;
    vector *eig_Values;

    if (argc != 4){
        print_input_err_exit();
    }

    goal = argv[2];
    
    data_points = filename_to_matrix(argv[3]);

    n = mat_get_nrows(data_points);

    if (strcmp(goal, "wam") == 0){
        W = wam(data_points);

        /* outputing the results */
        mat_print(W);

        /* freeing memory */
        free_matrix(W);

        return 0;
    }

    if (strcmp(goal, "ddg") == 0){
        DD = ddg(data_points);

        /* outputing the results */
        mat_print(DD);

        /* freeing memory */
        free_matrix(DD);

        return 0;
    }

    if (strcmp(goal, "lnorm") == 0){
        L_norm = lnorm(data_points);

        /* outputing the results */
        mat_print(L_norm);

        /* freeing memory */
        free_matrix(L_norm);

        return 0;
    }

    if (strcmp(goal, "jacobi") == 0){
        
        eig_Vectors = alloc_matrix(n, n, ~INITIALIZE);
        eig_Values = alloc_vector(n, ~VEC_INIT, ~VEC_IS_INT);

        jacobi(data_points, &eig_Values, &eig_Vectors, EPSILON, JACOBI_MAX_ITER);

        tmp = eig_Vectors;
        eig_Vectors = transpose(eig_Vectors);
        free_matrix(tmp);

        /* outputing the results */
        print_vector(eig_Values);
        mat_print(eig_Vectors);

        /* freeing memory */
        free_vector(eig_Values);
        free_matrix(eig_Vectors);
        
        return 0;
    }

    if (strcmp(goal, "spk") == 0){

        if (atof(argv[1]) != (int)(atof(argv[1]))){
            print_input_err_exit();
        }

        k = atof(argv[1]); /* auto-casting occures */

        if (k < 0){
            print_input_err_exit();
        }

        final_centroids = spk(data_points, 
                                k, 
                                NULL, 
                                EPSILON, 
                                JACOBI_MAX_ITER,
                                KMEANS_MAX_ITER);
        
        /* outputing the results */
        mat_print(final_centroids);

        /* freeing memory */
        
        free_matrix(final_centroids);

        return 0;
    }
    
    /* else: input's invalid */
    print_input_err_exit();
    
    return 1;
}










