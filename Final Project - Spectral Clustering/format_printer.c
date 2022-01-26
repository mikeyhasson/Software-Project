
#include "format_printer.h"

#include <stdio.h>
#include <math.h>

void print_double(double c, char sep, int with_sep){
    if (with_sep == WITH_SEP){
        if ((c < 0) && (c > -0.00005)){
            printf("%.4f%c", 0.0, sep);
        }
        else{
            printf("%.4f%c", c, sep);
        }
    }else{
        if ((c < 0) && (c > -0.00005)){
            printf("%.4f", 0.0);
        }
        else{
            printf("%.4f", c);
        }
    }
}
