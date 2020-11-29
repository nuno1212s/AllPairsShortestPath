#include "matrixutil.h"


#include <memory.h>
#include <stdlib.h>

int matrix_alloc(unsigned int size, unsigned int **matrix) {
    if (*matrix != NULL) {
        return 0;
    }

    *matrix = (unsigned int *) malloc(size * size * sizeof(unsigned int));
    return 1;
}

int matrix_copy(unsigned int size, unsigned int *dest, unsigned int *source) {
    if (dest == NULL || source == NULL)
        return 0 ;

    memcpy(dest, source, size * size * sizeof(unsigned int));

    return 1;
}

int matrix_free(unsigned int **mx) {
    if (*mx != NULL)
        free(*mx);
    else
        return 0;

    *mx = NULL;
    return 1;
}

int matrix_fill(unsigned int matrixSize, unsigned int *matrix, unsigned int value) {
    if (matrix == NULL)
        return 0 ;

    for (int i = 0; i != matrixSize * matrixSize; ++i)
        matrix[i] = value;

    return 1;
}