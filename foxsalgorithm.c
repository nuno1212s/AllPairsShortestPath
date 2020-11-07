#include "foxsalgorithm.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

int subdivideMatrix(unsigned int matrixSize, unsigned int (*matrix)[matrixSize],
                    unsigned int processCount, unsigned int Q,
                unsigned int (*matrixDestinations)[matrixSize / Q][matrixSize / Q]) {

    int perMatrixSize = matrixSize / Q;

    int currentRowIndex = 0;

    for (int matrixCounter = 0; matrixCounter < processCount; matrixCounter++) {

        int startRow = currentRowIndex * perMatrixSize, endRow = startRow + perMatrixSize,
                startColumn = (matrixCounter % Q) * perMatrixSize, endColumn = startColumn + perMatrixSize;

        for (int row = 0; row < perMatrixSize; row++) {
            for (int column = 0; column < perMatrixSize; column++) {
                matrixDestinations[matrixCounter][row][column] = matrix[row + startRow][column + startColumn];
            }
        }

        if ((matrixCounter + 1) % Q == 0) {
            currentRowIndex++;
        }
    }

    return 1;
}