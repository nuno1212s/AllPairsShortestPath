#include "matrix_multiplication.h"
#include "matrixutil.h"
#include <limits.h>
#include <memory.h>
#include <stdio.h>

#define min(a, b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

int prepareMatrixForAllPairs(unsigned int matrixSize, unsigned int *matrix) {

    for (int i = 0; i < matrixSize; i++) {
        for (int j = 0; j < matrixSize; j++) {

            if (i != j && matrix[PROJECT(matrixSize, i, j)] == 0) {
                //Use only INT_MAX so that dik + dkj doesn't surpass the limit of unsigned ints
                matrix[PROJECT(matrixSize, i, j)] = INT_MAX - 1;
            }

        }
    }

    return 1;
}

int doSpecialMatrixMultiply(unsigned int matrixSize, unsigned int *matrixA,
                            unsigned int *matrixB,
                            unsigned int *matrixC) {

    for (int i = 0; i < matrixSize; i++) {

        for (int j = 0; j < matrixSize; j++) {

            for (int k = 0; k < matrixSize; k++) {

                matrixC[PROJECT(matrixSize, i, j)] = min(matrixC[PROJECT(matrixSize, i, j)]
                        , matrixA[PROJECT(matrixSize, i, k)] + matrixB[PROJECT(matrixSize, k, j)]);

            }
        }

    }

    return 1;
}

int doAllPairsShortestPaths(unsigned int matrixSize, unsigned int *matrixW,
                            unsigned int *matrixD) {

    int m = 1;

    matrix_fill(matrixSize, matrixD, INT_MAX - 1);

    while (m < matrixSize - 1) {

        doSpecialMatrixMultiply(matrixSize, matrixW, matrixW, matrixD);

        m *= 2;

        //Move the matrixD (D^2m) back to matrixW
        matrix_copy(matrixSize, matrixW, matrixD);
    }

    return 1;
}