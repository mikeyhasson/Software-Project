
#include "error_msg.h"

#include <stdio.h>
#include <stdlib.h>

void print_input_err_exit(){
    printf("Invalid Input!");
    exit(1);
}

void print_err_exit(){
    printf("An Error Has Occured");
    exit(1);
}
