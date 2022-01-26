#ifndef FILETRANS_H
#define FILETRANS_H

#include <stdlib.h>
#include "matrix.h"
#include "error_msg.h"

#define INITIAL_DATA_DIM 10
#define INITIAL_DATA_CONTAINER_SIZE 50

int* filename_to_indices(char *filename, int *ind_num);

matrix* filename_to_matrix(char *filename);

#endif
