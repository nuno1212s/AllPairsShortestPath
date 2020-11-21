#ifndef TRABALHO_1_MATRIX_MULTIPLICATION_H
#define TRABALHO_1_MATRIX_MULTIPLICATION_H

/**
 * Prepare a matrix for the all pairs algorithm
 * @param matrixSize The size of the matrix to prepare
 * @param matrix The matrix to prepare
 * @return Whether the preparation was successful
 */
int prepareMatrixForAllPairs(unsigned int matrixSize, unsigned int (*matrix)[matrixSize]);

/**
 * Perform the special matrix multiply
 *
 * @param matrixSize The size of the matrix
 * @param matrixA The matrix A
 * @param matrixB The matrix B
 * @param matrixC The destination matrix C
 * @return Whether the multiplication was successful
 */
int doSpecialMatrixMultiply(unsigned int matrixSize, unsigned int (*matrixA)[matrixSize],
                            unsigned int(*matrixB)[matrixSize],
                            unsigned int (*matrixC)[matrixSize]);

/**
 * Perform the all pairs shortest paths algorithm
 *
 * @param matrixSize The size of the matrix
 * @param matrixW The weight matrix, after being passed through prepareMatrix (This matrix will be modified)
 * @param matrixD The result matrix.
 * @return Whether the algorithm was successful
 */
int doAllPairsShortestPaths(unsigned int matrixSize, unsigned int (*matrixW)[matrixSize],
                            unsigned int (*matrixD)[matrixSize]);

#endif //TRABALHO_1_MATRIX_MULTIPLICATION_H
