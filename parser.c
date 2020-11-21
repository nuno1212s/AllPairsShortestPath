#include "parser.h"
#include <stdlib.h>

int parseMatrix(FILE *fp, unsigned int matrixSize, unsigned int (*matrix)[matrixSize]) {

    for (int row = 0; row < matrixSize; row++) {

        for (int column = 0; column < matrixSize; column++) {

            fscanf(fp, "%d", &matrix[row][column]);

        }
    }

    return 1;
}

int printMatrix(FILE *fp, unsigned int matrixSize,unsigned int (*matrix)[matrixSize]) {

    for (int row = 0; row < matrixSize; row++) {
        for (int column = 0; column < matrixSize; column++) {

            fprintf(fp, "%d ", matrix[row][column]);

        }

        fprintf(fp, "\n");
    }

    return 1;
}