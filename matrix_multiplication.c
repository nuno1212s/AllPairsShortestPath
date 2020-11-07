#include "matrix_multiplication.h"
#include <limits.h>
#include <memory.h>
#include <stdio.h>

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

int prepareMatrixForAllPairs(unsigned int matrixSize, unsigned int (*matrix)[matrixSize]) {

    printf("prepping...\n");

    for (int i = 0; i < matrixSize; i++) {
        for (int j = 0; j < matrixSize; j++) {

            if (i != j && matrix[i][j] == 0) {
                matrix[i][j] = INT_MAX;
            }

        }
    }

    return 1;
}

int doSpecialMatrixMultiply(unsigned int matrixSize, unsigned int (*matrixA)[matrixSize],
                        unsigned int (*matrixB)[matrixSize],
                        unsigned int (*matrixC)[matrixSize]) {

    for (int i = 0; i < matrixSize; i++) {

        for (int j = 0; j < matrixSize; j++) {

            matrixC[i][j] = UINT_MAX;

            for (int k = 0; k < matrixSize; k++) {

                matrixC[i][j] = min(matrixC[i][j], matrixA[i][k] + matrixB[k][j]);

            }

        }

    }

    return 1;
}

int doAllPairsShortestPaths(unsigned int matrixSize, unsigned int (*matrixW)[matrixSize],
                            unsigned int (*matrixD)[matrixSize]) {

    int m = 1;

    while (m < matrixSize - 1) {

        doSpecialMatrixMultiply(matrixSize, matrixW, matrixW, matrixD);

        m *= 2;

        //Move the matrixD (D^2m) back to matrixW
        memcpy(matrixW, matrixD, sizeof(int) * matrixSize * matrixSize);
    }

    return 1;
}