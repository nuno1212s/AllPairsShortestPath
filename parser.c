#include "parser.h"
#include <stdlib.h>
#include <limits.h>
#include "matrixutil.h"

int parseMatrix(FILE *fp, unsigned int matrixSize, unsigned int *matrix) {

    for (int row = 0; row < matrixSize; row++) {

        for (int column = 0; column < matrixSize; column++) {

            fscanf(fp, "%d", &matrix[PROJECT(matrixSize, row, column)]);

        }
    }

    return 1;
}

int printMatrix(FILE *fp, unsigned int matrixSize, unsigned int *matrix) {

    for (int row = 0; row < matrixSize; row++) {
        for (int column = 0; column < matrixSize; column++) {

            int pos = PROJECT(matrixSize, row, column);

            if (matrix[pos] >= INT_MAX - 2) {
                fprintf(fp, "IM ");
            } else {
                fprintf(fp, "%d ", matrix[pos]);
            }

        }

        fprintf(fp, "\n");
    }

    return 1;
}